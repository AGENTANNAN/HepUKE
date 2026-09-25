# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")        # J/psi(3097) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")     # corresponding inclusive MC sample
qed_data   = DatasetManager.real_data.find("708_3080")        # 3.080 GeV off-resonance data used for the QED continuum background

# Single decay card covering all three radiative modes J/psi -> gamma P with P -> gamma gamma.
# The three channels share the identical all-photon final state (3 gamma), so one algorithm/card is used.
decay_card_for_signal = <<~DECAYCARD
    Decay J/psi
    0.3333 gamma pi0  HELAMP 1 0 1 0;
    0.3333 gamma eta  HELAMP 1 0 1 0;
    0.3334 gamma eta' HELAMP 1 0 1 0;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Signal exclusive MC: 2.3 million events, generated with the three-mode decay card
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gammaP_3gamma"
  config.related_dataset = jpsi_data            # associated real data set for the simulation
  config.events          = 2_300_000            # 2.3 M signal events
  config.decay_card      = decay_card_for_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiGammaP"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})          # CMS energy 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# Full event selection chain
event_selection = Selection.new
event_selection
  .select_track {                       # Charged-track quality cuts + full charged-track veto
     cos_theta 0.93                     # |cos(theta)| < 0.93 for any track
     Vz        10.0                     # |Vz| < 10 cm
     Vr        1.0                      # Vr < 1 cm in the transverse plane
     nChrp     "==0"                    # no positively charged track survives -> veto
     nChrn     "==0"                    # no negatively charged track survives -> veto
     nTot      "==0"                    # exactly zero charged tracks
  }
  .select_photon {                      # Photon selection (purely neutral final state)
     tdc_emc_start     -500.0           # EMC time window: -500 ns
     tdc_emc_end        500.0           # EMC time window: +500 ns
     energyThreshold_b  0.080           # 80 MeV threshold in the barrel region
     energyThreshold_e  0.080           # 80 MeV threshold in the endcap region
     nGam               ">=3"           # at least three photons
  }
  # No charged-particle PID: the final state is purely neutral.
  # Nominal 4C kinematic fit to gamma gamma gamma
  .kinematic_fit([:gamma, :gamma, :gamma]) {
     nominal                            # nominal fit: corrected four-momenta are saved
     constrain_four_momentum            # 4C energy-momentum constraint
     chi2_cut 50                        # chi2(3gamma) < 50
  }
  # Competing hypothesis 1: e+e- -> gamma gamma (2-gamma 4C fit, chi2 stored for the veto chi2(3gamma) < chi2(2gamma))
  .kinematic_fit([:gamma, :gamma]) {
     constrain_four_momentum
  }
  # Competing hypothesis 2: pi0 pi0 (4-gamma 4C fit, chi2 stored for the veto chi2(3gamma) < chi2(4gamma))
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma]) {
     constrain_four_momentum
  }

# Generate the algorithm for the three-mode decay card and execute on all datasets
my_algorithm.with_decay_card(decay_card_for_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal, qed_data])