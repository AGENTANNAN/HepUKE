# ==========================================================================
# Dataset preparation
# ==========================================================================
# Real data and inclusive MC at sqrt(s) = 4.600 GeV (the energy on which the
# event selection below is executed).
data_4600  = DatasetManager.real_data.find("703_4600")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

# Additional energy points prepared for the cross-section scan.
# They are loaded here but NOT executed: the selection runs on the 4.600 GeV data.
data_4575 = DatasetManager.real_data.find("703_4575")   # 4.575 GeV
data_4527 = DatasetManager.real_data.find("703_4530")   # 4.527 GeV
data_4467 = DatasetManager.real_data.find("703_4470")   # 4.467 GeV
data_4416 = DatasetManager.real_data.find("703_4420")   # 4.416 GeV

# ==========================================================================
# Decay cards (EvtGen format)
# ==========================================================================
# Channel I (Dbar*0): e+e- -> Ds+ Dbar*0 K-, dominated by Ds1(2536)- -> Dbar*0 K-,
# with Ds+ -> K+ K- pi+,  Dbar*0 -> Dbar0 pi0,  pi0 -> gamma gamma
decay_card_DsDstar0K = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Ds+ anti-D*0 K-      PHSP;
    Enddecay

    Decay Ds+
    1.0000 K+ K- pi+            PHSP;
    Enddecay

    Decay anti-D*0
    1.0000 anti-D0 pi0          PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma          PHSP;
    Enddecay

    End
DECAYCARD

# Channel II (Dbar0): e+e- -> Ds+ Dbar0 K-, dominated by Ds2*(2573)- -> Dbar0 K-,
# with Ds+ -> K+ K- pi+
decay_card_DsD0K = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Ds+ anti-D0 K-       PHSP;
    Enddecay

    Decay Ds+
    1.0000 K+ K- pi+            PHSP;
    Enddecay

    End
DECAYCARD

# ==========================================================================
# Exclusive signal MC: 500k events for each of the two modes
# ==========================================================================
exMC_DsDstar0K = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_DsDstar0K_4600"
  config.related_dataset = data_4600
  config.events          = 500000
  config.decay_card      = decay_card_DsDstar0K
  config.cross_section   = :default
end

exMC_DsD0K = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_DsD0K_4600"
  config.related_dataset = data_4600
  config.events          = 500000
  config.decay_card      = decay_card_DsD0K
  config.cross_section   = :default
end

# ==========================================================================
# Channel I: e+e- -> Ds+ Dbar*0 K-
# ==========================================================================
alg_I = Algorithm.new("DsDstar0K")
alg_I.set_header(["DsDstar0KAlg/DsDstar0K.h"])
     .set_constant({"ECMS" => [:double, 4.600]})
     .note(:ds_mass_window, "Ds+ -> K+K-pi+ candidates required in [1.955, 1.980] GeV/c2")
     .note(:background_veto, "Lambda_c+ background suppressed by requiring the Ds+ recoil mass RQ(Ds+) < 2.59 GeV/c2")
     .note(:sub_resonance_regions, "Ds+ sub-resonance regions A (phi pi+: M(K+K-) < 1.05 GeV/c2) and B (K+ K*0: 0.863 < M(K-pi+) < 0.930 GeV/c2) combined for the cross-section measurement")
     .note(:fit_model, "unbinned maximum-likelihood fit to RQ(Ds+) for the Ds1(2536)- and Ds2*(2573)- parameters (S/D-wave Breit-Wigner convolved with resolution); ISR and vacuum-polarization corrections applied to the Born cross section")
     .note(:spin_parity, "Ds2*(2573)- J^P = 2+ is favored")

sel_I = Selection.new
sel_I.select_track {
        cos_theta  0.93
        Vz         10.0
        Vr         1.0
        nChrp      ">=2"
        nChrn      ">=2"
        nNet       "==0"
      }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon]
        identify :kaon, against: [:pion]
      }
     .partial_rec([1, 3, 4, 5, 6]) {
        best_combination_by_mass :Dsp, 1.96834
        require_recoil_mass 1.993, 2.024
      }

alg_I.with_decay_card(decay_card_DsDstar0K).apply(sel_I)

# ==========================================================================
# Channel II: e+e- -> Ds+ Dbar0 K-
# ==========================================================================
alg_II = Algorithm.new("DsD0K")
alg_II.set_header(["DsD0KAlg/DsD0K.h"])
      .set_constant({"ECMS" => [:double, 4.600]})
      .note(:ds_mass_window, "Ds+ -> K+K-pi+ candidates required in [1.955, 1.980] GeV/c2")
      .note(:background_veto, "Lambda_c+ suppressed by using only Ds+ sub-resonance region A (phi pi+: M(K+K-) < 1.05 GeV/c2)")
      .note(:fit_model, "unbinned maximum-likelihood fit to RQ(Ds+) for the Ds1(2536)- and Ds2*(2573)- parameters (S/D-wave Breit-Wigner convolved with resolution); ISR and vacuum-polarization corrections applied to the Born cross section")
      .note(:spin_parity, "Ds2*(2573)- J^P = 2+ is favored")

sel_II = Selection.new
sel_II.select_track {
         cos_theta  0.93
         Vz         10.0
         Vr         1.0
         nChrp      ">=2"
         nChrn      ">=2"
         nNet       "==0"
       }
      .pid(method: :probability) {
         prob_cut 0.001
         identify :pion, against: [:kaon]
         identify :kaon, against: [:pion]
       }
      .partial_rec([1, 3, 4, 5, 6]) {
         best_combination_by_mass :Dsp, 1.96834
         require_recoil_mass 1.850, 1.880
       }

alg_II.with_decay_card(decay_card_DsD0K).apply(sel_II)

# ==========================================================================
# Execute on the 4.600 GeV data, inclusive MC and the two signal MC samples
# ==========================================================================
root_files_I  = alg_I.execute_on([data_4600, incMC_4600, exMC_DsDstar0K])
root_files_II = alg_II.execute_on([data_4600, incMC_4600, exMC_DsD0K])