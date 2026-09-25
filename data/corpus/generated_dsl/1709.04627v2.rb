# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi(3097) real data sample (1.31e9 events)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC sample

# Decay card: charged mode J/psi -> gamma eta', eta' -> eta pi+ pi-, eta -> gamma gamma
decay_card_charged = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'  PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi-  PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: neutral mode J/psi -> gamma eta', eta' -> eta pi0 pi0, eta -> gamma gamma, pi0 -> gamma gamma
decay_card_neutral = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'  PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi0 pi0  PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 500k events for each eta' decay mode
exMC_charged = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etap_etapipi"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_charged
  config.cross_section   = :default
end

exMC_neutral = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etap_etapi0pi0"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_neutral
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# ================= Mode I (charged): eta' -> eta pi+ pi-, eta -> gamma gamma =================
alg_name_charged = "JpsiGammaEtapCharged"
alg_charged = Algorithm.new(alg_name_charged)
alg_charged.set_header(["#{alg_name_charged}Alg/#{alg_name_charged}.h"])
           .set_constant({"ECMS" => [:double, 3.097]})       # J/psi(3097) CMS energy in GeV
           .set_alias({"std::vector<double>" => "Vdouble"})

sel_charged = Selection.new
  .select_track {                       # exactly two oppositely charged tracks
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     "==1"                     # exactly 1 positive track
    nChrn     "==1"                     # exactly 1 negative track
    nNet      "==0"                     # net charge zero
  }
  .pid(method: :probability) {          # probability PID
    prob_cut 0.001                      # prob > 0.001
    identify :pion, against: [:kaon, :proton]  # pi+ and pi- vs K and p
    npip "==1"                          # one pi+
    npim "==1"                          # one pi-
  }
  .select_photon {                      # at least three photons
    tdc_emc_start     0                 # EMC timing 0 ... 700 ns
    tdc_emc_end       14
    angle_to_track    10.0              # > 10 deg photon-track opening angle
    energyThreshold_b 0.025             # 25 MeV barrel
    energyThreshold_e 0.050             # 50 MeV endcap
    nGam ">=3"                          # at least 3 photons
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C Kalman fit: eta -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"                          # at least one eta candidate
  }
  .kinematic_fit([:gamma, :eta, :pip, :pim]) {  # 6C fit to gamma eta pi+ pi-
    nominal
    constrain_four_momentum             # 4C
    invariant_mass_of(:eta).constrain_to_nominal_mass_of(:eta)             # + m(eta)
    invariant_mass_of(:eta, :pip, :pim).constrain_to_nominal_mass_of(:etap) # + m(eta')
    chi2_cut 200                        # loose chi2, tight cut applied in ROOT
  }

# The most energetic photon is assigned as the radiative photon from J/psi -> gamma eta'
alg_charged.note(:radiative_photon_assignment,
  "the most energetic photon in the event is assigned as the radiative photon from " \
  "J/psi -> gamma eta'; the remaining photon clusters form the eta (-> gamma gamma) candidate")

alg_charged.with_decay_card(decay_card_charged).apply(sel_charged)

# ================= Mode II (neutral): eta' -> eta pi0 pi0, eta -> gamma gamma, pi0 -> gamma gamma =================
alg_name_neutral = "JpsiGammaEtapNeutral"
alg_neutral = Algorithm.new(alg_name_neutral)
alg_neutral.set_header(["#{alg_name_neutral}Alg/#{alg_name_neutral}.h"])
           .set_constant({"ECMS" => [:double, 3.097]})       # J/psi(3097) CMS energy in GeV
           .set_alias({"std::vector<double>" => "Vdouble"})

sel_neutral = Selection.new
  .select_track {                       # no charged tracks
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==0"
    nChrn     "==0"
    nNet      "==0"
  }
  .select_photon {                      # at least seven photons
    tdc_emc_start     0                 # EMC timing 0 ... 700 ns
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025             # 25 MeV barrel
    energyThreshold_e 0.050             # 50 MeV endcap
    nGam ">=7"                          # at least 7 photons
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C Kalman fit: pi0 -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 "==2"                          # exactly two pi0 candidates
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C Kalman fit: eta -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"                          # at least one eta candidate
  }
  .kinematic_fit([:gamma, :eta, :pi0, :pi0]) {  # 8C fit to gamma eta pi0 pi0
    nominal
    constrain_four_momentum             # 4C
    invariant_mass_of(:eta).constrain_to_nominal_mass_of(:eta)              # + m(eta)
    invariant_mass_of(:pi0).constrain_to_nominal_mass_of(:pi0)              # + m(pi0) (x2)
    invariant_mass_of(:eta, :pi0, :pi0).constrain_to_nominal_mass_of(:etap) # + m(eta')
    chi2_cut 200                        # loose chi2, tight cut applied in ROOT
  }

alg_neutral.note(:radiative_photon_assignment,
    "the most energetic photon in the event is assigned as the radiative photon from " \
    "J/psi -> gamma eta'; the remaining photon clusters form the eta (-> gamma gamma) and " \
    "the two pi0 (-> gamma gamma) candidates")
  .note(:background_veto,
    "pi0 candidates required to satisfy |cos(theta_decay)| < 0.95 (photon decay angle in the " \
    "pi0 rest frame) to suppress photon-pair miscombination; applied after the Kalman pi0 reconstruction")

alg_neutral.with_decay_card(decay_card_neutral).apply(sel_neutral)

### Execute on datasets and produce ROOT files ###
root_files_charged = alg_charged.execute_on([jpsi_data, jpsi_incMC, exMC_charged])
root_files_neutral = alg_neutral.execute_on([jpsi_data, jpsi_incMC, exMC_neutral])