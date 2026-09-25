# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # corresponding inclusive MC sample

# --- Decay cards (EvtGen format) for the four signal modes ---

# Mode I: psi(2S) -> pi0 hc, hc -> gamma pi+ pi-
decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 hc PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay hc
    1.000 gamma pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: psi(2S) -> pi0 hc, hc -> gamma pi+ pi- eta, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 hc PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay hc
    1.000 gamma pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode III: psi(2S) -> pi0 hc, hc -> gamma pi+ pi- pi+ pi-
decay_card_modeIII = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 hc PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay hc
    1.000 gamma pi+ pi- pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode IV: psi(2S) -> pi0 hc, hc -> gamma p+ anti-p-
decay_card_modeIV = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 hc PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay hc
    1.000 gamma p+ anti-p- PHSP;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC samples (100k events per mode) ---
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exmc_3686_pi0hc_gamPipPim"
    config.related_dataset = psip_data
    config.events         = 100000
    config.decay_card     = decay_card_modeI
    config.cross_section  = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exmc_3686_pi0hc_gamPipPimEta"
    config.related_dataset = psip_data
    config.events         = 100000
    config.decay_card     = decay_card_modeII
    config.cross_section  = :default
end

exMC_modeIII = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exmc_3686_pi0hc_gam4Pi"
    config.related_dataset = psip_data
    config.events         = 100000
    config.decay_card     = decay_card_modeIII
    config.cross_section  = :default
end

exMC_modeIV = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exmc_3686_pi0hc_gamPpbar"
    config.related_dataset = psip_data
    config.events         = 100000
    config.decay_card     = decay_card_modeIV
    config.cross_section  = :default
end

### Event selection (BOSS) — one Algorithm + Selection per decay mode ###

# ===================== Mode I: hc -> gamma pi+ pi- =====================
alg_name_I = "HcToGamPipPim"
alg_modeI = Algorithm.new(alg_name_I)
alg_modeI.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
sel_modeI.select_track {              # charged-track selection: pi+ pi-
    cos_theta 0.93
    Vz        100.0
    Vr        10.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {                    # photon selection: radiative gamma + pi0 -> gamma gamma
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"           # at least three photons
  }
  .pid(method: :probability) {        # pions identified against kaons
    prob_cut 0.001
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C Kalman fit: gamma gamma -> pi0 mass
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pi0]) {  # 5C fit: 4C + pi0 mass constraint
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

# ===================== Mode II: hc -> gamma pi+ pi- eta (eta -> gamma gamma) =====================
alg_name_II = "HcToGamPipPimEta"
alg_modeII = Algorithm.new(alg_name_II)
alg_modeII.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
sel_modeII.select_track {             # charged-track selection: pi+ pi-
    cos_theta 0.93
    Vz        100.0
    Vr        10.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {                    # photon selection: radiative gamma + pi0 gamma gamma + eta gamma gamma
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=5"           # at least five photons
  }
  .pid(method: :probability) {        # pions identified against kaons
    prob_cut 0.001
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C Kalman fit: gamma gamma -> pi0 mass
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C Kalman fit: gamma gamma -> eta mass
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pi0, :eta]) {  # 5C fit: 4C + pi0 mass (eta constrained)
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

# ===================== Mode III: hc -> gamma 2(pi+ pi-) =====================
alg_name_III = "HcToGam4Pi"
alg_modeIII = Algorithm.new(alg_name_III)
alg_modeIII.set_header(["#{alg_name_III}Alg/#{alg_name_III}.h"])
           .set_constant({"ECMS" => [:double, 3.686]})
           .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeIII = Selection.new
sel_modeIII.select_track {            # charged-track selection: 2(pi+ pi-)
    cos_theta 0.93
    Vz        100.0
    Vr        10.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {                    # photon selection: radiative gamma + pi0 -> gamma gamma
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"           # at least three photons
  }
  .pid(method: :probability) {        # pions identified against kaons
    prob_cut 0.001
    identify :pion, against: [:kaon]
    npip "==2"
    npim "==2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C Kalman fit: gamma gamma -> pi0 mass
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :pip, :pip, :pim, :pim, :pi0]) {  # 5C fit: 4C + pi0 mass constraint
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_modeIII.with_decay_card(decay_card_modeIII).apply(sel_modeIII)

# ===================== Mode IV: hc -> gamma p+ anti-p- =====================
alg_name_IV = "HcToGamPpbar"
alg_modeIV = Algorithm.new(alg_name_IV)
alg_modeIV.set_header(["#{alg_name_IV}Alg/#{alg_name_IV}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeIV = Selection.new
sel_modeIV.select_track {             # charged-track selection: p+ anti-p-
    cos_theta 0.93
    Vz        100.0
    Vr        10.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {                    # photon selection: radiative gamma + pi0 -> gamma gamma
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"           # at least three photons
  }
  .pid(method: :probability) {        # protons identified against kaons and pions
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"
    nprm "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C Kalman fit: gamma gamma -> pi0 mass
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :prp, :prm, :pi0]) {  # 5C fit: 4C + pi0 mass constraint
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_modeIV.with_decay_card(decay_card_modeIV).apply(sel_modeIV)

### Execute each algorithm on data, inclusive MC, and its exclusive MC ###
root_files_modeI   = alg_modeI.execute_on([psip_data, psip_incMC, exMC_modeI])
root_files_modeII  = alg_modeII.execute_on([psip_data, psip_incMC, exMC_modeII])
root_files_modeIII = alg_modeIII.execute_on([psip_data, psip_incMC, exMC_modeIII])
root_files_modeIV  = alg_modeIV.execute_on([psip_data, psip_incMC, exMC_modeIV])