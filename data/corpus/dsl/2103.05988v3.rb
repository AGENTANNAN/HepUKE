# Measurement of D -> K-pi+pi+pi- and D -> K-pi+pi0 coherence factors
# and average strong-phase differences in quantum-correlated DDbar decays
# BESIII: 2.93 fb-1 at psi(3770)

# --- Datasets ---
psip_data = DatasetManager.real_data.find("712_3773")
psip_incMC = DatasetManager.inclusive_mc.find("712_3773")

# --- Decay cards for two signal modes ---
# Mode I: D0 -> K-pi+pi+pi-
decay_card_K3pi = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay
    Decay D0
    1.000 K- pi+ pi+ pi- PHSP;
    Enddecay
    Decay anti-D0
    1.000 K+ pi- pi- pi+ PHSP;
    Enddecay
End
DECAYCARD

# Mode II: D0 -> K-pi+pi0
decay_card_Kpipi0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay
    Decay D0
    1.000 K- pi+ pi0 PHSP;
    Enddecay
    Decay anti-D0
    1.000 K+ pi- pi0 PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
End
DECAYCARD

exMC_K3pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0_K3pi_signal"
  config.related_dataset = psip_data
  config.events = 200000
  config.decay_card = decay_card_K3pi
  config.cross_section = :default
end

exMC_Kpipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0_Kpipi0_signal"
  config.related_dataset = psip_data
  config.events = 200000
  config.decay_card = decay_card_Kpipi0
  config.cross_section = :default
end

# === Mode I Algorithm: D -> K-pi+pi+pi- ===
alg_K3pi = Algorithm.new("D0toK3pi")
alg_K3pi.set_header(["D0toK3piAlg/D0toK3pi.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })

sel_K3pi = Selection.new
sel_K3pi.select_track {
           cos_theta 0.93
           Vz 100.0
           Vr 10.0
           nChrp ">=2"
           nChrn ">=2"
           nNet "==0"
         }
         .select_photon {
           tdc_emc_start 0
           tdc_emc_end 14
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           angle_to_track 20.0
         }
         .pid(method: :probability) {
           prob_cut 0.001
           identify :kaon, against: [:pion]
           identify :pion, against: [:kaon]
           nkm ">=1"
           nkp ">=1"
           npip ">=2"
           npim ">=2"
         }
         # D0 mass-constrained kinematic fit for signal D candidate
         .kinematic_fit([:km, :pip, :pip, :pim]) {
           nominal
           invariant_mass_of(:km, :pip, :pip, :pim).constrain_to_nominal_mass_of(:D0)
           chi2_cut 200
         }

alg_K3pi.with_decay_card(decay_card_K3pi).apply(sel_K3pi)
alg_K3pi.note(:ddbar_correlation, "Quantum-correlated DDbar pairs at psi(3770).
  Double-tag technique with three classes of tags:
  - CP tags (even: K+K-, pi+pi-, pi+pi-pi0, KsPi0, KsOmega, KsEta';
     odd: KsPi, KsEta, KsOmega_alt, KsEta'_alt, KsPhi, Kspipi0);
     KL0 modes use missing-mass technique Mmiss^2 peaking at KL0 mass.
  - Like-sign and opposite-sign flavour tags (K-pi+, K-pi+pi+pi-, K-pi+pi0).
  - Self-conjugate KsPiPi tags in 16 Dalitz bins.
  Both D mesons fully reconstructed; yields extracted from MBC fits.")
  .note(:ks_veto, "Ks veto via flight significance L/sigma_L > 2 applied to
  pi+pi- pairs in D->K-pi+pi+pi- signal to suppress D->KsKpi background.
  Ks mass region [0.487,0.511] excluded from binned phase space.")
  .note(:deltaE_cut, "DeltaE = E_D - sqrt(s)/2 within +/- 3 sigma of DeltaE peak;
  resolution varies by mode: from [-18,+17] MeV (K3pi) to [-69,+44] MeV (KsPi0).")
  .note(:best_candidate, "If multiple DT candidates in one event, combination with
  average reconstructed D mass closest to nominal D0 mass is chosen.")
  .note(:binned_phase_space, "D->K-pi+pi+pi- phase space divided into 4 bins
  based on LHCb amplitude model; bin boundaries chosen to equalise product of
  integrated CF and DCS amplitudes across bins.")
  .note(:peaking_backgrounds, "Doubly misidentified opposite-sign events calibrated
  via data-driven misID rates in momentum bins. D->KsKpi background constrained
  using known branching fractions and quantum-correlation corrections.")
  .note(:cp_tag_fractions, "CP-even fraction F+ used for quasi-CP-eigenstate tags:
  pi+pi-pi0 (F+ > 95%). Non-resonant KsPiPiPi0 backgrounds in KsOmega/KsEta tags
  at 20%/10% level, mildly CP-odd / CP-neutral.")

# === Mode II Algorithm: D -> K-pi+pi0 ===
alg_Kpipi0 = Algorithm.new("D0toKpipi0")
alg_Kpipi0.set_header(["D0toKpipi0Alg/D0toKpipi0.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })

sel_Kpipi0 = Selection.new
sel_Kpipi0.select_track {
            cos_theta 0.93
            Vz 100.0
            Vr 10.0
            nChrp ">=1"
            nChrn ">=1"
            nNet "==0"
          }
          .select_photon {
            tdc_emc_start 0
            tdc_emc_end 14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            angle_to_track 20.0
            nGam ">=2"
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion]
            identify :pion, against: [:kaon]
            nkm ">=1"
            nkp ">=1"
            npip ">=1"
            npim ">=1"
          }
          # pi0 reconstruction via Kalman mass-constrained fit
          .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"
          }
          # D0 mass-constrained kinematic fit
          .kinematic_fit([:km, :pip, :pi0]) {
            nominal
            invariant_mass_of(:km, :pip, :pi0).constrain_to_nominal_mass_of(:D0)
            chi2_cut 200
          }

alg_Kpipi0.with_decay_card(decay_card_Kpipi0).apply(sel_Kpipi0)
alg_Kpipi0.note(:ddbar_correlation, "Same double-tag strategy as Mode I with CP tags,
  flavour tags, and KsPiPi tags. Quantum-correlated DDbar pairs at psi(3770).")
  .note(:deltaE_cut, "DeltaE within +/- 3 sigma of peak; resolution mode-dependent.")
  .note(:kl_missing_mass, "CP-tag modes containing KL0: Mmiss^2 calculated from
  signal D momentum and remaining event content; peaking at KL0 mass squared.")
  .note(:pi0_mass_window, "Di-photon invariant mass in [115,150] MeV/c^2 with
  at least one barrel photon; 1C kinematic fit to nominal pi0 mass.")
  .note(:cosmic_veto, "For D->K+K- and D->pipi- tag modes: TOF time difference
  between two tracks < 5 ns, neither track identified as e/mu, to suppress
  cosmic and Bhabha backgrounds.")

alg_K3pi.execute_on([psip_data, psip_incMC, exMC_K3pi])
alg_Kpipi0.execute_on([psip_data, psip_incMC, exMC_Kpipi0])