# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # 1.31e9 J/psi events at sqrt(s) = 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive MC sample

### Decay cards (EvtGen format) ###
# Mode I: J/psi -> gamma eta, eta -> pi+ pi- pi0 (eta Dalitz), pi0 -> gamma gamma
decay_card_eta_pipipim0 = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta PHSP;
    Enddecay

    Decay eta
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: J/psi -> gamma eta, eta -> pi0 pi0 pi0, pi0 -> gamma gamma
decay_card_eta_3pi0 = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta PHSP;
    Enddecay

    Decay eta
    1.0000 pi0 pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode III: J/psi -> gamma eta', eta' -> pi0 pi0 pi0, pi0 -> gamma gamma
decay_card_etap_3pi0 = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 pi0 pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples (200k events each) ###
exMC_eta_pipipim0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_eta_pipipim0"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_eta_pipipim0
  config.cross_section   = :default
end

exMC_eta_3pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_eta_3pi0"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_eta_3pi0
  config.cross_section   = :default
end

exMC_etap_3pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etap_3pi0"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_etap_3pi0
  config.cross_section   = :default
end

### Event selection (BOSS) — Mode I: eta -> pi+ pi- pi0 ###
alg_name_eta_pipipim0 = "JpsiGammaEtaToPiPiPi0"
alg_eta_pipipim0 = Algorithm.new(alg_name_eta_pipipim0)
alg_eta_pipipim0.set_header(["#{alg_name_eta_pipipim0}Alg/#{alg_name_eta_pipipim0}.h"])
                 .set_constant({"ECMS" => [:double, 3.097]})
                 .set_alias({"std::vector<double>" => "Vdouble"})

sel_eta_pipipim0 = Selection.new
  # Charged track selection: exactly one pi+ and one pi-, net charge zero
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  # Photon selection: at least three good photons
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"
  }
  # Nominal 6C fit: 4C + m(gamma gamma)=m(pi0) + m(gamma gamma pi+ pi-)=m(eta)
  .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma, :pip, :pim).constrain_to_nominal_mass_of(:eta)
    chi2_cut 80
  }
  # Competing 4C hypotheses: store chi2 (the signal-vs-competitor veto is applied later in ROOT)
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }
  .kinematic_fit([:pip, :pim, :gamma, :gamma]) {
    constrain_four_momentum
  }

alg_eta_pipipim0
  .note(:radiative_photon_assignment, "the most energetic photon is taken as the radiative photon; the remaining two photons are combined into the pi0 candidate")
  .with_decay_card(decay_card_eta_pipipim0)
  .apply(sel_eta_pipipim0)

alg_eta_pipipim0.execute_on([jpsi_data, jpsi_incMC, exMC_eta_pipipim0])

### Event selection (BOSS) — Mode II: eta -> 3 pi0 ###
alg_name_eta_3pi0 = "JpsiGammaEtaTo3Pi0"
alg_eta_3pi0 = Algorithm.new(alg_name_eta_3pi0)
alg_eta_3pi0.set_header(["#{alg_name_eta_3pi0}Alg/#{alg_name_eta_3pi0}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

sel_eta_3pi0 = Selection.new
  # Require no charged tracks
  .select_track {
    nChrp "==0"
    nChrn "==0"
  }
  # At least seven good photons
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=7"
  }
  # Reconstruct pi0 from photon pairs (1C mass constraint); require at least three pi0
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=3"
  }
  # Nominal 7C fit to gamma pi0 pi0 pi0: 4C + m(3 pi0) = m(eta)
  .kinematic_fit([:gamma, :pi0, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:pi0, :pi0, :pi0).constrain_to_nominal_mass_of(:eta)
    chi2_cut 70
  }

alg_eta_3pi0
  .note(:pi0_decay_angle_cut, "each pi0 candidate is required to have |cos(theta_decay)| < 0.95 in its gamma-gamma rest frame")
  .with_decay_card(decay_card_eta_3pi0)
  .apply(sel_eta_3pi0)

alg_eta_3pi0.execute_on([jpsi_data, jpsi_incMC, exMC_eta_3pi0])

### Event selection (BOSS) — Mode III: eta' -> 3 pi0 ###
alg_name_etap_3pi0 = "JpsiGammaEtaPrimeTo3Pi0"
alg_etap_3pi0 = Algorithm.new(alg_name_etap_3pi0)
alg_etap_3pi0.set_header(["#{alg_name_etap_3pi0}Alg/#{alg_name_etap_3pi0}.h"])
             .set_constant({"ECMS" => [:double, 3.097]})
             .set_alias({"std::vector<double>" => "Vdouble"})

sel_etap_3pi0 = Selection.new
  # Require no charged tracks
  .select_track {
    nChrp "==0"
    nChrn "==0"
  }
  # At least seven good photons
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=7"
  }
  # Reconstruct pi0 from photon pairs (1C mass constraint); require at least three pi0
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=3"
  }
  # Nominal 7C fit to gamma pi0 pi0 pi0: 4C + m(3 pi0) = m(eta'); veto |m(gamma pi0) - m(omega)| < 0.05 GeV
  .kinematic_fit([:gamma, :pi0, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:pi0, :pi0, :pi0).constrain_to_nominal_mass_of(:etap)
    invariant_mass_of(:gamma, :pi0).out_of(0.7327, 0.8327)
    chi2_cut 70
  }
  # Competing 7C hypothesis J/psi -> gamma eta pi0 pi0 (eta -> gamma gamma): store chi2 for the ROOT-level veto
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  }

alg_etap_3pi0
  .note(:pi0_decay_angle_cut, "each pi0 candidate is required to have |cos(theta_decay)| < 0.95 in its gamma-gamma rest frame")
  .note(:background_veto, "reject events in which any gamma-gamma pair has invariant mass within 0.03 GeV of m(eta)")
  .with_decay_card(decay_card_etap_3pi0)
  .apply(sel_etap_3pi0)

alg_etap_3pi0.execute_on([jpsi_data, jpsi_incMC, exMC_etap_3pi0])