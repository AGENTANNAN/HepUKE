### Dataset preparation ###
# J/psi (3.097 GeV) real data and its inclusive MC
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card Mode I: J/psi -> eta' e+ e-, eta' -> gamma pi+ pi-
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000 eta' e+ e- PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card Mode II: J/psi -> eta' e+ e-, eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.0000 eta' e+ e- PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for Mode I (200k events)
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_etaprime_ee_gammapipi"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

# Exclusive MC for Mode II (100k events)
exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_etaprime_ee_pipieta"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# ===================== Mode I: eta' -> gamma pi+ pi- =====================
alg_name_modeI = "JPsiEtaPrimeEEGammaPiPi"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
  .select_track {          # charged track selection: exactly 2 positive, 2 negative, net charge 0
    cos_theta 0.93         # |cos(theta)| < 0.93
    Vz        100.0        # |Vz| < 100 cm
    Vr        10.0         # Vr < 10 cm
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {         # photon selection: at least one photon (Mode I)
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025   # 25 MeV in the barrel
    energyThreshold_e 0.050   # 50 MeV in the endcap
    angle_to_track    10.0    # > 10 deg from any charged track
    nGam              ">=1"
  }
  .pid(method: :probability) {   # PID: probability method, prob cut 0.001
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6  # p>1.0 -> lepton; EMC E>0.6 -> e, else mu
    identify :pion, against: [:kaon]  # pi/K separation; remaining tracks are the pi+ pi-
    nlp   "==1"
    nlm   "==1"
    npip  "==1"
    npim  "==1"
  }
  .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {   # nominal 4C fit to gamma pi+ pi- e+ e-
    nominal
    vertex_fit([1, 2, 3, 4])                          # vertex constraint on the four charged tracks
    constrain_four_momentum                           # 4C energy-momentum constraint
    invariant_mass_of(:gamma, :lp, :lm).out_of(0.10, 0.16)  # pi0 Dalitz veto M(gamma e+ e-) not in [0.10,0.16]
    chi2_cut 100
  }

alg_modeI
  .note(:background_veto, "gamma -> e+ e- conversion veto (transverse distance delta_xy < 2 cm) applied to all photon candidates; no dedicated DSL primitive exists for the conversion-vertex distance")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)
alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])

# ===================== Mode II: eta' -> pi+ pi- eta, eta -> gamma gamma =====================
alg_name_modeII = "JPsiEtaPrimeEEPipiEta"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
  .select_track {          # same charged track selection
    cos_theta 0.93
    Vz        100.0
    Vr        10.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {         # photon selection: at least two photons (Mode II)
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  }
  .pid(method: :probability) {   # same PID
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]
    nlp   "==1"
    nlm   "==1"
    npip  "==1"
    npim  "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # Kalman fit to reconstruct eta from gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)  # eta mass constraint
    invariant_mass_of(:gamma, :gamma).between(0.48, 0.60)                 # M(gamma gamma) in [0.48,0.60]
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:eta, :pip, :pim, :lp, :lm]) {   # nominal 4C fit to eta pi+ pi- e+ e-
    nominal
    vertex_fit([1, 2, 3, 4])                       # vertex constraint on the four charged tracks
    constrain_four_momentum
    chi2_cut 100
  }

alg_modeII
  .note(:background_veto, "gamma -> e+ e- conversion veto (transverse distance delta_xy < 2 cm) applied to all photon candidates; no dedicated DSL primitive exists for the conversion-vertex distance")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)
alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])