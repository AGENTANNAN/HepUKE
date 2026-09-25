# ============================================================
# Dataset preparation
# ============================================================
# Multi-energy scan: 56 CMS points, sqrt(s) = 3.51 - 4.95 GeV, 44.55 fb^-1.
# Use ConExc-style scan with representative BESIII datasets.
data_3773  = DatasetManager.real_data.find("712_3773")
data_4180  = DatasetManager.real_data.find("703_4180")
data_4230  = DatasetManager.real_data.find("703_4230")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# --------------------------------------------------------------------
# Signal decay card: e+e- -> K+ Xi0 anti-Sigma-, PHSP,
# Xi0 -> Lambda pi0, Lambda -> p pi-, pi0 -> gamma gamma,
# anti-Sigma- -> anti-n pi- (inferred as missing recoil).
# Continuum production; use psi(4260) as the top KKMC mother per BESIII convention.
# --------------------------------------------------------------------
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000  K+  Xi0  anti-Sigma-        PHSP;
  Enddecay

  Decay Xi0
  1.0000  Lambda0  pi0                PHSP;
  Enddecay

  Decay anti-Xi0
  1.0000  anti-Lambda0  pi0           PHSP;
  Enddecay

  Decay Lambda0
  1.0000  p+  pi-                     PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000  anti-p-  pi+                PHSP;
  Enddecay

  Decay anti-Sigma-
  1.0000  anti-n0  pi+                PHSP;
  Enddecay

  Decay Sigma+
  1.0000  n0  pi+                     PHSP;
  Enddecay

  Decay pi0
  1.0000  gamma  gamma                PHSP;
  Enddecay

  End
DECAYCARD

# Signal exclusive MC generated at 10^5 events per c.m. point in the scan
exMC_signal_scan = DatasetManager.create_exclusive_mc_for([data_3773, data_4180, data_4230]) do |config|
  config.sample_name    = "ee_KpXi0SigmaBarM"
  config.events         = 100_000
  config.decay_card     = decay_card_signal
  config.cross_section  = :default
end
exMC_signal_scan.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

# ============================================================
# Event selection: partial reconstruction technique.
# Reconstruct K+ and Xi0 (-> Lambda pi0 -> p pi- gamma gamma);
# infer anti-Sigma- from recoil of K+ Xi0 system.
# ============================================================
alg = Algorithm.new("EEtoKpXi0SigmaBarM")
alg.set_header(["EEtoKpXi0SigmaBarMAlg/EEtoKpXi0SigmaBarM.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })

event_selection = Selection.new
event_selection.select_track {
                 cos_theta 0.93
                 # No IP requirement (Lambda gives a secondary vertex)
                 nChrp ">=2"      # at least K+ and p
                 nChrn ">=1"      # at least pi-
               }
               .select_photon {
                 energyThreshold_b 0.025
                 energyThreshold_e 0.050
                 tdc_emc_start 0
                 tdc_emc_end   14
                 angle_to_track 10.0
                 nGam ">=2"
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :proton, against: [:kaon, :pion]
                 identify :kaon,   against: [:pion, :proton]
                 identify :pion,   against: [:kaon, :proton]
                 nprp ">=1"
                 nkp  ">=1"
                 npim ">=1"
               }
               # Reconstruct pi0 from gamma gamma via 1C mass constraint
               .kalman_kinematic_fit([:gamma, :gamma]) {
                 invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                 chi2_cut 500
                 npi0 ">=1"
               }
               # Reconstruct Lambda -> p pi- via vertex + secondary vertex fit
               .secondary_vertex_fit([:prp, :pim]) {
                 build_virtual_particle(:Lambda).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               }
               # Partial reconstruction: reconstruct K+ side (kp) + Xi0 (-> Lambda pi0);
               # anti-Sigma- inferred from recoil of the K+ Xi0 system.
               # DecayCardResolver.rec_id_list (assumed layout):
               #   0 => psi(4260)  (top mother, skip)
               #   1 => K+                                     (tag)
               #   2 => Xi0        -> Lambda pi0               (tag)
               #   3 => anti-Sigma-                            (missing recoil)
               #   ...
               .partial_miss([3]) {
                 best_combination_by_mass :Xi0, 1.31486
                 # Recoil against K+ Xi0 should sit near the anti-Sigma- mass (~1.197 GeV)
                 require_recoil_mass 1.10, 1.30
               }

alg.note(:lambda_mass_window,
         "|M(p pi-) - m_Lambda| <= 5 MeV/c^2, applied in ROOT; L/deltaL > 2 to suppress backgrounds.")
   .note(:xi0_mass_window,
         "|M(Lambda pi0) - m_Xi0| < 10 MeV/c^2; best Xi0 candidate chosen by min |Delta M|.")
   .note(:isr_vp_correction,
         "Cross section obtained iteratively with ISR (1+delta) and VP correction 1/|1-Pi|^2 " \
         "following Refs. [51-53]; final Born cross section extracted from fit to M_recoil(K+ Xi0).")

alg.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on all datasets: real data at multiple energy points, inclusive MC, and the scanned signal MC
alg.execute_on([data_3773, data_4180, data_4230, incMC_3773] + exMC_signal_scan)
