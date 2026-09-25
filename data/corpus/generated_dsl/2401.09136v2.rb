# Core DSL classes and dependencies are loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

# Decay card for J/psi -> gamma eta, eta -> gamma e+ e-   (final state: gamma gamma e+ e-)
decay_card_eta = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for J/psi -> gamma eta', eta' -> gamma e+ e-  (final state: gamma gamma e+ e-)
decay_card_etap = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 100k events for each of the two modes
exMC_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_eta_gammaee"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_eta
  config.cross_section   = :default
end

exMC_etap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etaprime_gammaee"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_etap
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name_eta = "JpsiGammaEtaGammaEE"
alg_eta = Algorithm.new(alg_name_eta)
alg_eta.set_header(["#{alg_name_eta}Alg/#{alg_name_eta}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})   # ECMS = 3.097 GeV

alg_name_etap = "JpsiGammaEtaPrimeGammaEE"
alg_etap = Algorithm.new(alg_name_etap)
alg_etap.set_header(["#{alg_name_etap}Alg/#{alg_name_etap}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})

# Common selection chain shared by both the eta and eta' channels
event_selection_common = Selection.new
  .select_track {
    cos_theta 0.93      # |cos(theta)| < 0.93
    Vz        10.0      # |Vz| < 10 cm
    Vr        1.0       # Vr < 1 cm in transverse plane
    nChrp     "==1"     # exactly one positive track
    nChrn     "==1"     # exactly one negative track
    nNet      "==0"     # net charge zero
  }
  .select_photon {
    tdc_emc_start     0      # TDC window 0-14
    tdc_emc_end       14
    energyThreshold_b 0.025  # barrel energy > 25 MeV
    energyThreshold_e 0.050  # endcap energy > 50 MeV
    angle_to_track    10.0   # angle to nearest charged track > 10 deg
    nGam              ">=2"  # at least two photons -> gamma gamma
  }

# eta channel: gamma gamma e+e-, electron/pion likelihood cut > 0.5
eta_selection = event_selection_common.dup
  .pid(method: :probability) {
    prob_cut 0.5
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
    nlp "==1"   # one positive lepton
    nlm "==1"   # one negative lepton
  }
  # Nominal 4C kinematic fit to the gamma gamma e+ e- hypothesis
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 100
  }

# eta' channel: same chain, electron/pion likelihood cut > 0.95
etap_selection = event_selection_common.dup
  .pid(method: :probability) {
    prob_cut 0.95
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  }
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 100
  }

# Generate the algorithms for each decay mode
alg_eta.with_decay_card(decay_card_eta).apply(eta_selection)
alg_etap.with_decay_card(decay_card_etap).apply(etap_selection)

# Execute on the datasets, producing ROOT files
root_files_eta  = alg_eta.execute_on([jpsi_data, jpsi_incMC, exMC_eta])
root_files_etap = alg_etap.execute_on([jpsi_data, jpsi_incMC, exMC_etap])