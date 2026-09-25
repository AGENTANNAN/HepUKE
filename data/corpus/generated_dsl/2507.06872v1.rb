# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# J/psi (3.097 GeV) real data and inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card: lepton-number-violating signal J/psi -> K+ K+ e- e-
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 K+ K+ e- e- PHSP;
    Enddecay
    End
DECAYCARD

# Decay card: charge-conjugate mode J/psi -> K- K- e+ e+
decay_card_signal_cc = <<~DECAYCARD
    Decay J/psi
    1.000 K- K- e+ e+ PHSP;
    Enddecay
    End
DECAYCARD

# 100k exclusive phase-space MC events for the signal
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_jpsi_to_kkee"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Analogous 100k exclusive MC for the charge-conjugate mode
exMC_signal_cc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_jpsi_to_kkee_cc"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_signal_cc
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ---------- Signal mode: J/psi -> K+ K+ e- e- ----------
alg_name_signal = "LNVKKee"
alg_signal = Algorithm.new(alg_name_signal)
alg_signal.set_header(["#{alg_name_signal}Alg/#{alg_name_signal}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})

selection_signal = Selection.new
selection_signal
  .select_track {                                   # charged track selection
     cos_theta 0.93                                 # |cos(theta)| < 0.93
     Vz        10.0                                 # |Vz| < 10 cm
     Vr        1.0                                  # Vr < 1 cm
     nChrp     "==2"                                # two positive tracks
     nChrn     "==2"                                # two negative tracks -> four tracks total
     nNet      "==0"                                # net charge zero
  }
  .pid(method: :probability) {                      # probability PID, 0.001 cut
     prob_cut 0.001
     identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.8
     identify :kaon, against: [:pion]               # pi/K separation
     nkp "==2"                                      # 2 K+
     nkm "==0"                                      # 0 K-
     nlp "==0"                                      # 0 l+
     nlm "==2"                                      # 2 l-
  }
  # nominal 4C fit to K+ K+ e- e-; J/psi mass window applied as a pre-fit selection
  .kinematic_fit([:kp, :kp, :lm, :lm]) {
     nominal
     constrain_four_momentum
     invariant_mass_of(:kp, :kp, :lm, :lm).within(3.07, 3.12)
     chi2_cut 20
  }
  # competing 4C hypotheses (chi2 stored for the ROOT-level "smallest-chi2" selection)
  .assign({:chrgp => :kp, :chrgn => :km})
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:kp, :km, :pip, :pim]) {          # K+ K- pi+ pi-
     constrain_four_momentum
  }
  .kinematic_fit([:kp, :km, :kp, :km]) {            # K+ K- K+ K-
     constrain_four_momentum
  }
  .kinematic_fit([:pip, :pim, :pip, :pim]) {        # pi+ pi- pi+ pi-
     constrain_four_momentum
  }

alg_signal
  .note(:background_veto,
        "e+e- pairs originating from gamma conversions are vetoed before the 4C kinematic fit to suppress gamma -> e+e- contamination of the K+K+e-e- final state")
  .note(:pid_criteria,
        "on top of the probability PID: electron requires CL_e>0.001, CL_e/(CL_e+CL_K+CL_pi)>0.8 and 0.8<E/p<1.2; kaon requires CL_K>0, CL_K>CL_pi and E/p<0.8 (electron veto)")
  .with_decay_card(decay_card_signal)
  .apply(selection_signal)

alg_signal.execute_on([jpsi_data, jpsi_incMC, exMC_signal])

# ---------- Charge-conjugate mode: J/psi -> K- K- e+ e+ (counters swapped) ----------
alg_name_cc = "LNVKKeeCC"
alg_cc = Algorithm.new(alg_name_cc)
alg_cc.set_header(["#{alg_name_cc}Alg/#{alg_name_cc}.h"])
      .set_constant({"ECMS" => [:double, 3.097]})
      .set_alias({"std::vector<double>" => "Vdouble"})

selection_cc = Selection.new
selection_cc
  .select_track {                                   # same track quality cuts
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
     nChrp     "==2"
     nChrn     "==2"
     nNet      "==0"
  }
  .pid(method: :probability) {                      # counters swapped wrt the signal mode
     prob_cut 0.001
     identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.8
     identify :kaon, against: [:pion]
     nkp "==0"
     nkm "==2"
     nlp "==2"
     nlm "==0"
  }
  # nominal 4C fit to K- K- e+ e-; J/psi mass window as a pre-fit selection
  .kinematic_fit([:km, :km, :lp, :lp]) {
     nominal
     constrain_four_momentum
     invariant_mass_of(:km, :km, :lp, :lp).within(3.07, 3.12)
     chi2_cut 20
  }
  # competing 4C hypotheses (chi2 stored for the ROOT-level "smallest-chi2" selection)
  .assign({:chrgp => :kp, :chrgn => :km})
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:kp, :km, :pip, :pim]) {          # K+ K- pi+ pi-
     constrain_four_momentum
  }
  .kinematic_fit([:kp, :km, :kp, :km]) {            # K+ K- K+ K-
     constrain_four_momentum
  }
  .kinematic_fit([:pip, :pim, :pip, :pim]) {        # pi+ pi- pi+ pi-
     constrain_four_momentum
  }

alg_cc
  .note(:background_veto,
        "e+e- pairs originating from gamma conversions are vetoed before the 4C kinematic fit to suppress gamma -> e+e- contamination of the K-K-e+e+ final state")
  .note(:pid_criteria,
        "on top of the probability PID: electron requires CL_e>0.001, CL_e/(CL_e+CL_K+CL_pi)>0.8 and 0.8<E/p<1.2; kaon requires CL_K>0, CL_K>CL_pi and E/p<0.8 (electron veto)")
  .with_decay_card(decay_card_signal_cc)
  .apply(selection_cc)

alg_cc.execute_on([jpsi_data, jpsi_incMC, exMC_signal_cc])