# BESIII DSL: 2405.12809v1 — Precision measurement of BF(J/psi -> K+K-) via psi(2S) -> pi+pi- J/psi
# psi(2S) data, 448.1e6 events, BOSS 709
# Relative to J/psi -> mu+mu-; no PID used, momentum-based track assignment

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: psi(2S) -> pi+ pi- J/psi, J/psi -> K+ K-
decay_card_sig = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi+ pi- J/psi JPIPI;
    Enddecay

    Decay J/psi
    1.000 K+ K- VSS;
    Enddecay
    End
DECAYCARD

exMC_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_kk_exclusive_mc"
  config.related_dataset = psip_data
  config.events = 750000
  config.decay_card = decay_card_sig
  config.cross_section = :default
end

# Decay card: psi(2S) -> pi+ pi- J/psi, J/psi -> mu+ mu- (reference channel)
decay_card_ref = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi+ pi- J/psi JPIPI;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- VLL;
    Enddecay
    End
DECAYCARD

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_mumu_reference_mc"
  config.related_dataset = psip_data
  config.events = 750000
  config.decay_card = decay_card_ref
  config.cross_section = :default
end

# =============================================
# Algorithm for J/psi -> K+K-  (same selection as mu+mu-, diff only in X_vis window)
# Steps expressible in DSL: track selection + common vertex fit
# Post-vertex-fit cuts (RM, cos_theta, momentum, E_dep/p, X_vis) are ROOT-level
# =============================================
alg_kk = Algorithm.new("JpsiKKViaPsi2S")
alg_kk.set_header(["JpsiKKViaPsi2SAlg/JpsiKKViaPsi2S.h"])
       .set_constant({ "ECMS" => [:double, 3.686] })
       .note(:no_pid, "No particle identification is applied; lower-momentum tracks (<1.0 GeV/c) are assigned as pions, higher-momentum tracks (>1.2 GeV/c) as kaons at ROOT level.")
       .note(:common_vertex_fit, "All 4 charged tracks are forced to a common vertex via vertex_fit(0,1,2,3) in the kinematic fit; fit must be successful.")
       .note(:electron_veto, "E_dep/p < 0.8 applied at ROOT level to veto J/psi -> e+e- background; no DSL method for per-track EMC energy deposit ratio.")
       .note(:emc_energy_window, "Total E_EMC required in [0.3, 2.5] GeV at ROOT level to suppress backgrounds.")
       .note(:cos_theta_cuts, "cos_theta(pipi) < 0.5 (photon conversion veto) and cos_theta(KK) < -0.95 applied at ROOT level.")
       .note(:recoil_mass_window, "RM(pi+pi-) in [3.087, 3.107] GeV/c^2 (J/psi tag) applied at ROOT level.")
       .note(:momentum_sum_veto, "|sum p_i|/E_Jpsi < 0.025 applied at ROOT level to reject >2-body decays.")
       .note(:xvis_signal_region, "X_vis = (E_K+ + E_K-)/E_Jpsi in [0.98, 1.01] (KKR) for J/psi->K+K- and [1.02, 1.07] (MMR) for J/psi->mu+mu-; applied at ROOT level.")
       .note(:xvis_min, "X_vis > 0.90 applied at ROOT level to reject low-mass tail (ppbar background).")

sel_kk = Selection.new
sel_kk.select_track {
         cos_theta 0.80
         Vz 10.0
         Vr 1.0
         nChrp "==2"
         nChrn "==2"
         nNet "==0"
       }
       # No photon selection needed (all-neutral final state reconstruction in this decay)
       # No PID: assign all positive to pip, negative to pim (kaon vs pion resolved by momentum at ROOT level)
       .assign({:chrgp => :pip, :chrgn => :pim})
       # Common vertex fit for all 4 tracks (the only BOSS-level kinematic constraint)
       .kinematic_fit([:pip, :pim, :pip, :pim]) {
         nominal
         vertex_fit([0, 1, 2, 3])
         chi2_cut 200
       }

alg_kk.with_decay_card(decay_card_sig).apply(sel_kk)
alg_kk.execute_on([psip_data, psip_incMC, exMC_kk, exMC_mumu])