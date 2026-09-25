### Datasets ###
# psi(3686) real data at 3.686 GeV (2712 x 10^6 psi(3686) events) and matching inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

### Decay cards (EvtGen format) ###
# Mode I: psi(3686) -> pi0 h_c, h_c -> pi+ pi- pi0   (bachelor pi0 and signal pi0 -> gamma gamma)
decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: psi(3686) -> pi0 h_c, h_c -> K+ K- pi0
decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 K+ K- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode III: psi(3686) -> pi0 h_c, h_c -> K+ K- eta
decay_card_modeIII = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 K+ K- eta PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode IV: psi(3686) -> pi0 h_c, h_c -> pi+ pi- eta
decay_card_modeIV = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 pi+ pi- eta PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples (500k events each, one per h_c decay mode) ###
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_pip_pim_pi0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_kp_km_pi0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

exMC_modeIII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_kp_km_eta"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_modeIII
  config.cross_section   = :default
end

exMC_modeIV = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_pip_pim_eta"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_modeIV
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Common charged-track and photon pre-selection shared by all four modes
base_selection = Selection.new
  .select_track {                 # charged track selection
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        10.0                # |Vz| < 10 cm
    Vr        1.0                 # Vr < 1 cm
    nChrp     "==1"               # exactly one positive track
    nChrn     "==1"               # exactly one negative track
    nNet      "==0"               # net charge zero
  }
  .select_photon {                # photon selection
    tdc_emc_start     0          # TDC start
    tdc_emc_end       14         # TDC end
    energyThreshold_b 0.025      # 25 MeV (barrel)
    energyThreshold_e 0.050      # 50 MeV (endcap)
    angle_to_track    10.0       # at least 10 deg from any charged track
    nGam              ">=4"      # bachelor pi0/eta and signal pi0/eta -> 4 photons
  }

# ------------------------------------------------------------------
# Mode I : psi(3686) -> pi0 h_c, h_c -> pi+ pi- pi0
# ------------------------------------------------------------------
alg_name_modeI = "PsiPpi0hcPipPimPi0"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

selection_modeI = base_selection.dup
  .pid(method: :probability) {                        # probability PID
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]         # pi+/pi- vs K and p
    npip ">=1"                                        # at least one pi+
    npim ">=1"                                        # at least one pi-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {           # reconstruct pi0 -> gamma gamma (1C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=2"                                        # bachelor pi0 + signal pi0 candidates
  }
  .kinematic_fit([:pip, :pim, :pi0, :pi0]) {          # nominal 6C: 4C + bachelor/signal pi0 masses
    nominal
    constrain_four_momentum
    chi2_cut 200                                      # loose BOSS cut (tight cuts applied in ROOT)
  }
  .kinematic_fit([:pip, :pim, :pi0]) {                # competing one-pi0 (fewer-photon) hypothesis
    constrain_four_momentum                           # no chi2_cut, no nominal -> store competing chi2
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(selection_modeI)

# ------------------------------------------------------------------
# Mode II : psi(3686) -> pi0 h_c, h_c -> K+ K- pi0
# ------------------------------------------------------------------
alg_name_modeII = "PsiPpi0hcKpKmPi0"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .set_alias({"std::vector<double>" => "Vdouble"})

selection_modeII = base_selection.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]         # K+/K- vs pi and p
    nkp ">=1"                                         # at least one K+
    nkm ">=1"                                         # at least one K-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=2"
  }
  .kinematic_fit([:kp, :km, :pi0, :pi0]) {            # nominal 6C
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  .kinematic_fit([:kp, :km, :pi0]) {                  # competing one-pi0 hypothesis
    constrain_four_momentum
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(selection_modeII)

# ------------------------------------------------------------------
# Mode III : psi(3686) -> pi0 h_c, h_c -> K+ K- eta
# ------------------------------------------------------------------
alg_name_modeIII = "PsiPpi0hcKpKmEta"
alg_modeIII = Algorithm.new(alg_name_modeIII)
alg_modeIII.set_header(["#{alg_name_modeIII}Alg/#{alg_name_modeIII}.h"])
           .set_constant({"ECMS" => [:double, 3.686]})
           .set_alias({"std::vector<double>" => "Vdouble"})

selection_modeIII = base_selection.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp ">=1"
    nkm ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {           # bachelor pi0 -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {           # signal eta -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta ">=1"
  }
  .kinematic_fit([:kp, :km, :pi0, :eta]) {            # nominal 6C: 4C + bachelor pi0 and signal eta masses
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  .kinematic_fit([:kp, :km, :pi0]) {                  # competing one-pi0 (fewer-photon) hypothesis
    constrain_four_momentum
  }

alg_modeIII.with_decay_card(decay_card_modeIII).apply(selection_modeIII)

# ------------------------------------------------------------------
# Mode IV : psi(3686) -> pi0 h_c, h_c -> pi+ pi- eta
# ------------------------------------------------------------------
alg_name_modeIV = "PsiPpi0hcPipPimEta"
alg_modeIV = Algorithm.new(alg_name_modeIV)
alg_modeIV.set_header(["#{alg_name_modeIV}Alg/#{alg_name_modeIV}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .set_alias({"std::vector<double>" => "Vdouble"})

selection_modeIV = base_selection.dup
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=1"
    npim ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {           # bachelor pi0 -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {           # signal eta -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta ">=1"
  }
  .kinematic_fit([:pip, :pim, :pi0, :eta]) {          # nominal 6C: 4C + bachelor pi0 and signal eta masses
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  .kinematic_fit([:pip, :pim, :pi0]) {                # competing one-pi0 (fewer-photon) hypothesis
    constrain_four_momentum
  }

alg_modeIV.with_decay_card(decay_card_modeIV).apply(selection_modeIV)

### Execute on real data, inclusive MC and the mode-specific exclusive MC ###
root_files_modeI   = alg_modeI.execute_on([psip_data, psip_incMC, exMC_modeI])
root_files_modeII  = alg_modeII.execute_on([psip_data, psip_incMC, exMC_modeII])
root_files_modeIII = alg_modeIII.execute_on([psip_data, psip_incMC, exMC_modeIII])
root_files_modeIV  = alg_modeIV.execute_on([psip_data, psip_incMC, exMC_modeIV])