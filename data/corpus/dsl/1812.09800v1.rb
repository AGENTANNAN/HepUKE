# arXiv: 1812.09800v1
# Observation of e+e- -> Ds+ Dbar(*0) K- and study of the P-wave Ds mesons
# BESIII Collaboration

# ============================================================================
# Dataset preparation
# ============================================================================

# Primary data: sqrt(s) = 4.600 GeV, 567 pb^-1
ds_4600 = DatasetManager.real_data.find("703_4600")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

# Cross-section scan datasets
ds_4575 = DatasetManager.real_data.find("703_4575")  # 4.575 GeV, 48 pb^-1
ds_4530 = DatasetManager.real_data.find("703_4530")  # 4.527 GeV, 110 pb^-1
ds_4470 = DatasetManager.real_data.find("703_4470")  # 4.467 GeV, 110 pb^-1
ds_4420 = DatasetManager.real_data.find("703_4420")  # 4.416 GeV, 1029 pb^-1

# ============================================================================
# Channel A: e+e- -> Ds+ Dbar*0 K-
# Dbar*0 identified via recoil mass (partial reconstruction)
# Dominated by Ds+ Ds1(2536)-, Ds1(2536)- -> Dbar*0 K-
# ============================================================================

decay_card_dstar = <<~DECAYCARD
  Decay psi(4260)
  1.000  D_s+ anti-D*0 K-   PHSP;
  Enddecay

  Decay D_s+
  1.000  K+ K- pi+          PHSP;
  Enddecay

  Decay anti-D*0
  1.000  anti-D0 pi0        VSS;
  Enddecay

  Decay pi0
  1.000  gamma gamma        PHSP;
  Enddecay

  End
DECAYCARD

exMC_dstar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_DsDstar0K"
  config.related_dataset = ds_4600
  config.events          = 500_000
  config.decay_card      = decay_card_dstar
  config.cross_section   = :default
end

algorithm_dstar = Algorithm.new("DsDstar0K")
algorithm_dstar.set_header(["DsDstar0KAlg/DsDstar0K.h"])
               .set_constant("ECMS" => [:double, 4.600])

selection_dstar = Selection.new
selection_dstar
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=2"
    nNet "==0"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon]
    identify :kaon, against: [:pion]
  }
  .partial_miss([2]) {
    best_combination_by_mass :D_s_p, 1.96834
  }

algorithm_dstar
  .note(:ds_mass_window, "Ds+ mass window (1.955, 1.980) GeV/c^2 applied to M(K+K-pi+)")
  .note(:ds_subresonance, "Ds+ sub-resonance selection: region A M(K+K-) < 1.05 GeV/c^2 (phi pi+), region B 0.863 < M(K-pi+) < 0.930 GeV/c^2 (K+ K*0); both A+B combined for cross-section measurement")
  .note(:lambda_c_veto, "Require RQ(Ds+) < 2.59 GeV/c^2 to suppress e+e- -> Lambda_c+ anti-Lambda_c- background")
  .note(:rq_dstar_window, "RQ(K-Ds+) signal region for anti-D*0: (1.993, 2.024) GeV/c^2")
  .note(:resonance_fit_method, "Ds1(2536)- resonance parameters extracted via unbinned ML fit to RQ(Ds+) spectrum; S-wave + D-wave BW convolved with resolution; cross section uses ISR correction (KKMC) and vacuum polarization factor")
  .with_decay_card(decay_card_dstar)
  .apply(selection_dstar)

# ============================================================================
# Channel B: e+e- -> Ds+ Dbar0 K-
# Dbar0 identified via recoil mass (partial reconstruction)
# Dominated by Ds+ Ds2*(2573)-, Ds2*(2573)- -> Dbar0 K-
# ============================================================================

decay_card_dbar0 = <<~DECAYCARD
  Decay psi(4260)
  1.000  D_s+ anti-D0 K-    PHSP;
  Enddecay

  Decay D_s+
  1.000  K+ K- pi+          PHSP;
  Enddecay

  End
DECAYCARD

exMC_dbar0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_DsDbar0K"
  config.related_dataset = ds_4600
  config.events          = 500_000
  config.decay_card      = decay_card_dbar0
  config.cross_section   = :default
end

algorithm_dbar0 = Algorithm.new("DsDbar0K")
algorithm_dbar0.set_header(["DsDbar0KAlg/DsDbar0K.h"])
               .set_constant("ECMS" => [:double, 4.600])

selection_dbar0 = Selection.new
selection_dbar0
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=2"
    nNet "==0"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon]
    identify :kaon, against: [:pion]
  }
  .partial_miss([2]) {
    best_combination_by_mass :D_s_p, 1.96834
  }

algorithm_dbar0
  .note(:ds_mass_window, "Ds+ mass window (1.955, 1.980) GeV/c^2 applied")
  .note(:ds2_phi_region, "For Ds2*(2573)- study: only Ds+ candidates from phi region A (M(K+K-) < 1.05 GeV/c^2) to suppress Lambda_c+ background")
  .note(:rq_dbar0_window, "RQ(K-Ds+) signal region for anti-D0: (1.850, 1.880) GeV/c^2")
  .note(:resonance_fit_method, "Ds2*(2573)- resonance parameters extracted via unbinned ML fit to RQ(Ds+) spectrum; D-wave BW convolved with resolution; J^P determined to be 2+ via helicity angle distribution chi^2 test")
  .note(:spin_parity, "Ds2*(2573)- spin-parity J^P = 2+ determined; J^P = 1- disfavored (chi^2 = 278.67 vs 7.85)")
  .with_decay_card(decay_card_dbar0)
  .apply(selection_dbar0)

# ============================================================================
# Execute
# ============================================================================

root_dstar = algorithm_dstar.execute_on([ds_4600, incMC_4600, exMC_dstar])
root_dbar0 = algorithm_dbar0.execute_on([ds_4600, incMC_4600, exMC_dbar0])