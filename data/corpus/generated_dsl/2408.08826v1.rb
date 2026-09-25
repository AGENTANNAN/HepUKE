# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi real data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC sample

# Decay cards (EvtGen format) for J/psi -> gamma D0, one card per D0 decay mode
decay_card_mode1 = <<~DECAYCARD
    Decay J/psi
    1.000 gamma D0 PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay

    End
DECAYCARD

decay_card_mode2 = <<~DECAYCARD
    Decay J/psi
    1.000 gamma D0 PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_mode3 = <<~DECAYCARD
    Decay J/psi
    1.000 gamma D0 PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples: 500k events for each D0 decay mode
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gammaD0_KPi"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_mode1
  config.cross_section   = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gammaD0_KPiPi0"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_mode2
  config.cross_section   = :default
end

exMC_mode3 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gammaD0_KPiPiPi"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_mode3
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Each D0 decay mode has a different final state -> its own Algorithm object
alg_name_mode1 = "JpsiGammaD0KPi"
alg_mode1 = Algorithm.new(alg_name_mode1)
alg_mode1.set_header(["#{alg_name_mode1}Alg/#{alg_name_mode1}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

alg_name_mode2 = "JpsiGammaD0KPiPi0"
alg_mode2 = Algorithm.new(alg_name_mode2)
alg_mode2.set_header(["#{alg_name_mode2}Alg/#{alg_name_mode2}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

alg_name_mode3 = "JpsiGammaD0KPiPiPi"
alg_mode3 = Algorithm.new(alg_name_mode3)
alg_mode3.set_header(["#{alg_name_mode3}Alg/#{alg_name_mode3}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

# Common selection shared by all three modes
event_selection_common = Selection.new
  .select_track {
    cos_theta 0.93          # |cos(theta)| < 0.93
    Vz        10.0          # |Vz| < 10 cm
    Vr        1.0           # Vr < 1 cm
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025 # 25 MeV in the barrel
    energyThreshold_e 0.050 # 50 MeV in the endcap
  }

# --- Mode I: D0 -> K- pi+ ---
mode1_selection = event_selection_common.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]     # K+/K- separated from pions
    identify :pion, against: [:kaon]     # pi+/pi- separated from kaons
    nkm  "==1"    # one K-
    nkp  "==0"    # no K+
    npip "==1"    # one pi+
    npim "==0"    # no pi-
  }
  .kinematic_fit([:km, :pip, :gamma]) {  # nominal 4C fit to K- pi+ gamma
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  .assign({:km => :pim})                 # pi+pi-gamma hypothesis: treat K- as pi-
  .kinematic_fit([:pip, :pim, :gamma]) {
    constrain_four_momentum
  }
  .assign({:pip => :kp})                 # K+K-gamma hypothesis: treat pi+ as K+
  .kinematic_fit([:kp, :km, :gamma]) {
    constrain_four_momentum
  }
  .kinematic_fit([:pip, :pim, :gamma, :gamma]) {  # pi+pi-gammagamma hypothesis
    constrain_four_momentum
  }

# --- Mode II: D0 -> K- pi+ pi0 (pi0 -> gamma gamma) ---
mode2_selection = event_selection_common.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkm  "==1"    # one K-
    nkp  "==0"    # no K+
    npip "==1"    # one pi+
    npim "==0"    # no pi-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {  # 1C fit: pi0 from two photons
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  }
  .kinematic_fit([:km, :pip, :gamma, :gamma, :gamma]) {  # nominal 5C fit to K- pi+ gamma gamma gamma
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  }
  .assign({:km => :pim})                 # pi+pi-pi0 gamma hypothesis
  .kinematic_fit([:pip, :pim, :pi0, :gamma]) {
    constrain_four_momentum
  }
  .assign({:pip => :kp})                 # K+K-pi0 gamma hypothesis
  .kinematic_fit([:kp, :km, :pi0, :gamma]) {
    constrain_four_momentum
  }
  .kinematic_fit([:pip, :pim, :pi0, :gamma, :gamma]) {  # pi+pi-pi0 gamma gamma hypothesis
    constrain_four_momentum
  }

# --- Mode III: D0 -> K- pi+ pi+ pi- ---
mode3_selection = event_selection_common.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkm  "==1"    # one K-
    nkp  "==0"    # no K+
    npip "==2"    # two pi+
    npim "==1"    # one pi-
  }
  .kinematic_fit([:km, :pip, :pip, :pim, :gamma]) {  # nominal 4C fit to K- pi+ pi+ pi- gamma
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  .assign({:km => :pim})                 # pi+pi-pi+pi-gamma hypothesis: treat K- as pi-
  .kinematic_fit([:pip, :pip, :pim, :pim, :gamma]) {
    constrain_four_momentum
  }

# Generate the algorithm code for each mode/decay card
alg_mode1.with_decay_card(decay_card_mode1).apply(mode1_selection)
alg_mode2.with_decay_card(decay_card_mode2).apply(mode2_selection)
alg_mode3.with_decay_card(decay_card_mode3).apply(mode3_selection)

# Execute on real data, inclusive MC, and the corresponding exclusive MC sample
alg_mode1.execute_on([jpsi_data, jpsi_incMC, exMC_mode1])
alg_mode2.execute_on([jpsi_data, jpsi_incMC, exMC_mode2])
alg_mode3.execute_on([jpsi_data, jpsi_incMC, exMC_mode3])