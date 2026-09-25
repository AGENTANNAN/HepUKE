# frozen_string_literal: true

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC

# ---- Decay cards (EvtGen) for the five hadronic h_c modes ----
# Mode 1: psi(2S) -> pi0 h_c, h_c -> 3(pi+pi-)pi0
decay_card_mode1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 pi+ pi- pi+ pi- pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 2: psi(2S) -> pi0 h_c, h_c -> 2(pi+pi-)pi0 eta
decay_card_mode2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 pi+ pi- pi+ pi- pi0 eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 3: psi(2S) -> pi0 h_c, h_c -> 2(pi+pi-)eta
decay_card_mode3 = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 pi+ pi- pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 4: psi(2S) -> pi0 h_c, h_c -> p pbar
decay_card_mode4 = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 p+ anti-p- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 5: psi(2S) -> pi0 h_c, h_c -> 2(pi+pi-)omega, omega -> pi+pi-pi0
decay_card_mode5 = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 pi+ pi- pi+ pi- omega PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---- 500k-event exclusive MC samples, one per h_c mode ----
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_pi0hc_3pipi0"
  config.related_dataset = psip_data
  config.events         = 500000
  config.decay_card     = decay_card_mode1
  config.cross_section  = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_pi0hc_2pipi_pi0eta"
  config.related_dataset = psip_data
  config.events         = 500000
  config.decay_card     = decay_card_mode2
  config.cross_section  = :default
end

exMC_mode3 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_pi0hc_2pipi_eta"
  config.related_dataset = psip_data
  config.events         = 500000
  config.decay_card     = decay_card_mode3
  config.cross_section  = :default
end

exMC_mode4 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_pi0hc_ppbar"
  config.related_dataset = psip_data
  config.events         = 500000
  config.decay_card     = decay_card_mode4
  config.cross_section  = :default
end

exMC_mode5 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_pi0hc_2pipi_omega"
  config.related_dataset = psip_data
  config.events         = 500000
  config.decay_card     = decay_card_mode5
  config.cross_section  = :default
end

### Event selection (BOSS) ###

# ---------------------------------------------------------------
# Mode 1: 3(pi+pi-)pi0  (also shared by Mode 5: 2(pi+pi-)omega)
# ---------------------------------------------------------------
alg1_name = "Pi0HcCharged6Pi2Pi0"
alg_mode1 = Algorithm.new(alg1_name)
alg_mode1.set_header(["#{alg1_name}Alg/#{alg1_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

selection_mode1 = Selection.new
  .select_track {                # charged track selection
    cos_theta 0.93               # |cos(theta)| < 0.93
    Vz        10.0               # |Vz| < 10 cm
    Vr        1.0                # Vr < 1 cm
    nChrp     "==3"              # exactly 3 positive tracks (3 pi+)
    nChrn     "==3"              # exactly 3 negative tracks (3 pi-)
    nNet      "==0"              # net charge zero
  }
  .select_photon {               # photon selection
    tdc_emc_start     0          # EMC timing 0-700 ns
    tdc_emc_end       14
    angle_to_track    10.0       # > 10 deg from nearest charged track
    energyThreshold_b 0.025      # > 25 MeV in barrel
    energyThreshold_e 0.050      # > 50 MeV in endcap
    nGam              ">=4"      # 2 pi0 -> >=4 photons
  }
  .pid(method: :probability) {   # PID: pions against kaons/protons
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==3"
    npim "==3"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct pi0 from gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=2"                                 # two pi0 in the final state
  }
  .kinematic_fit([:pip, :pip, :pip, :pim, :pim, :pim, :pi0, :pi0]) {  # nominal 4C fit
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode1.note(:background_veto, "J/psi-related backgrounds suppressed by recoil-mass windows and photon-count hypothesis tests (specific windows not given in the description; final windows applied at ROOT level)")

alg_mode1.with_decay_card(decay_card_mode1).apply(selection_mode1)

# ---------------------------------------------------------------
# Mode 2: 2(pi+pi-)pi0 eta
# ---------------------------------------------------------------
alg2_name = "Pi0Hc4PiPi0Eta"
alg_mode2 = Algorithm.new(alg2_name)
alg_mode2.set_header(["#{alg2_name}Alg/#{alg2_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

selection_mode2 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"              # exactly 2 pi+
    nChrn     "==2"              # exactly 2 pi-
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=6"      # 2 pi0 + eta -> >=6 photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==2"
    npim "==2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # eta -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta ">=1"
  }
  .kinematic_fit([:pip, :pip, :pim, :pim, :pi0, :pi0, :eta]) {  # nominal 4C fit
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode2.note(:background_veto, "J/psi-related backgrounds suppressed by recoil-mass windows and photon-count hypothesis tests (specific windows not given in the description; final windows applied at ROOT level)")

alg_mode2.with_decay_card(decay_card_mode2).apply(selection_mode2)

# ---------------------------------------------------------------
# Mode 3: 2(pi+pi-)eta
# ---------------------------------------------------------------
alg3_name = "Pi0Hc4PiEta"
alg_mode3 = Algorithm.new(alg3_name)
alg_mode3.set_header(["#{alg3_name}Alg/#{alg3_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

selection_mode3 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"      # pi0 + eta -> >=4 photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==2"
    npim "==2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # eta -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta ">=1"
  }
  .kinematic_fit([:pip, :pip, :pim, :pim, :pi0, :eta]) {  # nominal 4C fit
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode3.note(:background_veto, "J/psi-related backgrounds suppressed by recoil-mass windows and photon-count hypothesis tests (specific windows not given in the description; final windows applied at ROOT level)")

alg_mode3.with_decay_card(decay_card_mode3).apply(selection_mode3)

# ---------------------------------------------------------------
# Mode 4: p pbar
# ---------------------------------------------------------------
alg4_name = "Pi0HcPPbar"
alg_mode4 = Algorithm.new(alg4_name)
alg_mode4.set_header(["#{alg4_name}Alg/#{alg4_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

selection_mode4 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"              # exactly 1 proton
    nChrn     "==1"              # exactly 1 anti-proton
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"      # pi0 -> >=2 photons
  }
  .pid(method: :probability) {   # PID: protons against kaons/pions
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"
    nprm "==1"
  }
  .remove(:prp) { condition "fabs(cos_theta_of(:prp)) > 0.8" }   # |cos(theta_p)| < 0.8
  .remove(:prm) { condition "fabs(cos_theta_of(:prm)) > 0.8" }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=1"
  }
  .kinematic_fit([:prp, :prm, :pi0]) {        # nominal 4C fit
    nominal
    vertex_fit([0, 1])                        # common charged-track vertex for p and pbar
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode4.note(:background_veto, "J/psi-related backgrounds suppressed by recoil-mass windows and photon-count hypothesis tests (specific windows not given in the description; final windows applied at ROOT level)")

alg_mode4.with_decay_card(decay_card_mode4).apply(selection_mode4)

# ---------------------------------------------------------------
# Mode 5: 2(pi+pi-)omega (omega -> pi+pi-pi0) — shares Mode 1 selection,
# extracted by a simultaneous fit to M(pi+pi-pi0) in omega signal/sideband
# regions (performed at ROOT level).
# ---------------------------------------------------------------
alg5_name = "Pi0Hc4PiOmega"
alg_mode5 = Algorithm.new(alg5_name)
alg_mode5.set_header(["#{alg5_name}Alg/#{alg5_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

alg_mode5.note(:background_veto, "J/psi-related backgrounds suppressed by recoil-mass windows and photon-count hypothesis tests (specific windows not given in the description; final windows applied at ROOT level)")
alg_mode5.note(:signal_extraction, "Mode 5 shares the 3(pi+pi-)pi0 (Mode 1) selection; the omega is extracted by a simultaneous fit to the pi+pi-pi0 invariant mass in omega signal and sideband regions, performed at the ROOT analysis level")

alg_mode5.with_decay_card(decay_card_mode5).apply(selection_mode1.dup)

### Generate job scripts and run on data / inclusive MC / signal MC ###
root_files_mode1 = alg_mode1.execute_on([psip_data, psip_incMC, exMC_mode1])
root_files_mode5 = alg_mode5.execute_on([psip_data, psip_incMC, exMC_mode5])
root_files_mode2 = alg_mode2.execute_on([psip_data, psip_incMC, exMC_mode2])
root_files_mode3 = alg_mode3.execute_on([psip_data, psip_incMC, exMC_mode3])
root_files_mode4 = alg_mode4.execute_on([psip_data, psip_incMC, exMC_mode4])