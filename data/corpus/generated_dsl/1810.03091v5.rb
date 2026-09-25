### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding inclusive MC

# Decay card — Mode I: J/psi -> e+ e- eta, eta -> gamma gamma
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000 e+ e- eta PHOTOS VLL;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card — Mode II: J/psi -> e+ e- eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.0000 e+ e- eta PHOTOS VLL;
    Enddecay

    Decay eta
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC — 100k events for each mode
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_ee_eta_gammagamma"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_ee_eta_pipimpi0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# ================= Mode I: J/psi -> e+ e- eta, eta -> gamma gamma =================
alg_name_I = "JpsiEtaGG"
alg_modeI = Algorithm.new(alg_name_I)
alg_modeI.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})           # J/psi c.m. energy
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
  .select_track {                                              # exactly 1 positive + 1 negative track (l+ l-)
    cos_theta 0.93                                             # |cos(theta)| < 0.93
    Vz        10.0                                             # |Vz| < 10 cm
    Vr        1.0                                              # Vr < 1 cm
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {                                             # at least 2 photons
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0                                     # >= 10 deg from any charged track
    energyThreshold_b 0.025                                    # E > 25 MeV (barrel)
    energyThreshold_e 0.050                                    # E > 50 MeV (endcap)
    nGam              ">=2"
  }
  .pid(method: :probability) {                                 # hybrid: high-momentum leptons
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6  # p>1.0 -> lepton; EMC E>0.6 -> e, else mu
    nlp "==1"
    nlm "==1"
  }
  .kinematic_fit([:lp, :lm, :gamma, :gamma]) {                 # nominal 4C fit of l+ l- gamma gamma
    nominal
    vertex_fit([0, 1])                                         # vertex fit on the two lepton tracks
    constrain_four_momentum
    chi2_cut 100
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

# ================= Mode II: J/psi -> e+ e- eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma =================
alg_name_II = "JpsiEtaPiPiPi0"
alg_modeII = Algorithm.new(alg_name_II)
alg_modeII.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})          # J/psi c.m. energy
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
  .select_track {                                              # exactly 2 positive + 2 negative tracks
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {                                             # at least 2 photons (from pi0 -> gamma gamma)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"
  }
  .pid(method: :probability) {                                 # high-momentum leptons + non-lepton tracks as pions
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon, :proton]
    nlp   "==1"
    nlm   "==1"
    npip  "==1"
    npim  "==1"
  }
  .kinematic_fit([:lp, :lm, :pip, :pim, :gamma, :gamma]) {     # nominal 4C fit of l+ l- pi+ pi- gamma gamma
    nominal
    vertex_fit([0, 1, 2, 3])                                   # vertex fit on the four charged tracks
    constrain_four_momentum
    chi2_cut 100
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

### Execution ###
root_files_modeI  = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])
root_files_modeII = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])