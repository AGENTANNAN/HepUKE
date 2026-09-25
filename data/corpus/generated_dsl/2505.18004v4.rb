# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Seven e+e- data points, 4.600 - 4.699 GeV (~4.5 fb^-1).  Sample name = [BOSS_version]_[Ecm in MeV]
data_points = [
  DatasetManager.real_data.find("703_4600"),   # Ecm = 4599.53 MeV
  DatasetManager.real_data.find("706_4610"),   # Ecm = 4611.86 MeV
  DatasetManager.real_data.find("706_4620"),   # Ecm = 4628.00 MeV
  DatasetManager.real_data.find("706_4640"),   # Ecm = 4640.91 MeV
  DatasetManager.real_data.find("706_4660"),   # Ecm = 4661.24 MeV
  DatasetManager.real_data.find("706_4680"),   # Ecm = 4681.92 MeV
  DatasetManager.real_data.find("706_4700"),   # Ecm = 4698.82 MeV
]

# Matching inclusive MC samples (one per energy point)
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
]

### Decay cards (EvtGen format) ###
# The e+e- system is generated with psi(4260) as the KKMC top mother (BESIII convention);
# the Lambda_c+ carries the signal decay, the anti-Lambda_c- side decays inclusively (dominant mode).

# --- Signal: Lambda_c+ -> Sigma+ eta, eta -> gamma gamma ---
decay_card_lc_sigma_eta_gg = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Sigma+ eta PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K+ pi- PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- Signal: Lambda_c+ -> Sigma+ eta, eta -> pi+ pi- pi0 ---
decay_card_lc_sigma_eta_3pi = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Sigma+ eta PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K+ pi- PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay

  Decay eta
  1.0000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- Signal: Lambda_c+ -> Sigma+ eta', eta' -> pi+ pi- eta, eta -> gamma gamma ---
decay_card_lc_sigma_etap = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Sigma+ eta' PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K+ pi- PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- Reference: Lambda_c+ -> Sigma+ pi0 ---
decay_card_lc_sigma_pi0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Sigma+ pi0 PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K+ pi- PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- Reference: Lambda_c+ -> Sigma+ omega, omega -> pi+ pi- pi0 ---
decay_card_lc_sigma_omega = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Sigma+ omega PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K+ pi- PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay

  Decay omega
  1.0000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC (200k events per mode, produced at each of the seven energy points) ###
exMC_lc_sigma_eta_gg = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_Lc2SigmaEta_etagg"
  config.events      = 200000
  config.decay_card  = decay_card_lc_sigma_eta_gg
  config.cross_section = :default
end

exMC_lc_sigma_eta_3pi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_Lc2SigmaEta_eta3pi"
  config.events      = 200000
  config.decay_card  = decay_card_lc_sigma_eta_3pi
  config.cross_section = :default
end

exMC_lc_sigma_etap = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_Lc2SigmaEtaP_etap2pipieta"
  config.events      = 200000
  config.decay_card  = decay_card_lc_sigma_etap
  config.cross_section = :default
end

exMC_lc_sigma_pi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_Lc2SigmaPi0"
  config.events      = 200000
  config.decay_card  = decay_card_lc_sigma_pi0
  config.cross_section = :default
end

exMC_lc_sigma_omega = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_Lc2SigmaOmega"
  config.events      = 200000
  config.decay_card  = decay_card_lc_sigma_omega
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ==========================================================================
# Mode 1: Lambda_c+ -> Sigma+ eta, eta -> gamma gamma, Sigma+ -> p pi0
# ==========================================================================
alg_name_1 = "Lc2SigmaEtaEta2gg"
alg_1 = Algorithm.new(alg_name_1)
alg_1.set_header(["#{alg_name_1}Alg/#{alg_name_1}.h"])
     .set_constant({"ECMS" => [:double, 4.6]})   # placeholder; see note below
     .note(:ebeam_per_run, "Data span seven energy points 4.600-4.699 GeV; the single ECMS "
           "constant is a placeholder - the 4C fit must use the per-run measured beam energy "
           "for each of the seven datasets.")

sel_1 = Selection.new
  .select_track {
    cos_theta 0.93     # |cos(theta)| < 0.93
    Vz        10.0     # |Vz| < 10 cm
    Vr        1.0      # Vr < 1 cm
    nChrp    ">=1"     # eta -> gamma gamma final state: at least one positive track
    nChrn    ">=1"     # ... and at least one negative track
  }
  .select_photon {
    energyThreshold_b 0.025   # EMC barrel E > 25 MeV
    energyThreshold_e 0.050   # EMC endcap E > 50 MeV
    tdc_emc_start     0       # TDC window 0-14
    tdc_emc_end       14
    angle_to_track    10.0    # > 10 deg from any charged track
    nGam             ">=4"    # at least four photons
  }
  .pid(method: :probability) {
    prob_cut 0.0
    identify :proton, against: [:kaon, :pion]   # p vs K/pi
    identify :kaon,   against: [:pion]          # K vs pi
    identify :pion,   against: [:kaon]          # pi vs K
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # reconstruct pi0 from a gamma gamma pair
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 100
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # reconstruct eta from a gamma gamma pair
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 100
    neta ">=1"
  }
  .kinematic_fit([:prp, :pi0, :eta]) {          # 4C fit on p pi0 + eta
    nominal
    constrain_four_momentum
    invariant_mass_of(:prp, :pi0).constrain_to_nominal_mass_of(:Sigma_p)  # M(p pi0) -> Sigma+
    invariant_mass_of(:prp, :pi0).within(1.174, 1.200)   # Sigma+ signal window (GeV/c^2)
    invariant_mass_of(:gamma, :gamma).within(0.500, 0.560)  # eta -> gamma gamma window
    chi2_cut 17
  }

alg_1.with_decay_card(decay_card_lc_sigma_eta_gg).apply(sel_1)

# ==========================================================================
# Mode 2: Lambda_c+ -> Sigma+ eta, eta -> pi+ pi- pi0, Sigma+ -> p pi0
# ==========================================================================
alg_name_2 = "Lc2SigmaEtaEta23pi"
alg_2 = Algorithm.new(alg_name_2)
alg_2.set_header(["#{alg_name_2}Alg/#{alg_name_2}.h"])
     .set_constant({"ECMS" => [:double, 4.6]})
     .note(:ebeam_per_run, "Per-run beam energy required for the 4C fit (seven energy points).")

sel_2 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp    ">=2"     # p and pi+ (three-pion final state)
    nChrn    ">=2"     # pi- (and room for charge-conjugate background tracks)
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    nGam             ">=4"
  }
  .pid(method: :probability) {
    prob_cut 0.0
    identify :proton, against: [:kaon, :pion]
    identify :kaon,   against: [:pion]
    identify :pion,   against: [:kaon]
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # reconstruct the two pi0 (Sigma+ and eta daughters)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 100
    npi0 ">=2"
  }
  .kinematic_fit([:prp, :pi0, :pip, :pim, :pi0]) {  # 4C fit on p pi0 + eta(-> pi+ pi- pi0)
    nominal
    constrain_four_momentum
    invariant_mass_of(:prp, :pi0).constrain_to_nominal_mass_of(:Sigma_p)
    invariant_mass_of(:prp, :pi0).within(1.174, 1.200)          # Sigma+ window
    invariant_mass_of(:pip, :pim, :pi0).within(0.535, 0.560)    # eta -> 3pi window
    invariant_mass_of(:pi0, :pi0).out_of(0.440, 0.520)          # p K_S0 veto
    chi2_cut 17
  }

alg_2.with_decay_card(decay_card_lc_sigma_eta_3pi).apply(sel_2)

# ==========================================================================
# Mode 3: Lambda_c+ -> Sigma+ eta', eta' -> pi+ pi- eta, eta -> gamma gamma
# ==========================================================================
alg_name_3 = "Lc2SigmaEtaP"
alg_3 = Algorithm.new(alg_name_3)
alg_3.set_header(["#{alg_name_3}Alg/#{alg_name_3}.h"])
     .set_constant({"ECMS" => [:double, 4.6]})
     .note(:ebeam_per_run, "Per-run beam energy required for the 4C fit (seven energy points).")

sel_3 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp    ">=2"     # p and pi+ (three-pion final state)
    nChrn    ">=2"     # pi-
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    nGam             ">=4"
  }
  .pid(method: :probability) {
    prob_cut 0.0
    identify :proton, against: [:kaon, :pion]
    identify :kaon,   against: [:pion]
    identify :pion,   against: [:kaon]
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # reconstruct pi0 (Sigma+ daughter)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 100
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # reconstruct eta (eta' -> pi+ pi- eta)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 100
    neta ">=1"
  }
  .kinematic_fit([:prp, :pi0, :pip, :pim, :eta]) {  # 4C fit on p pi0 + eta'(-> pi+ pi- eta)
    nominal
    constrain_four_momentum
    invariant_mass_of(:prp, :pi0).constrain_to_nominal_mass_of(:Sigma_p)
    invariant_mass_of(:prp, :pi0).within(1.174, 1.200)          # Sigma+ window
    invariant_mass_of(:pip, :pim, :eta).within(0.946, 0.968)    # eta' window
    chi2_cut 30
  }

alg_3.with_decay_card(decay_card_lc_sigma_etap).apply(sel_3)

# ==========================================================================
# Reference mode 4: Lambda_c+ -> Sigma+ pi0, Sigma+ -> p pi0, pi0 -> gamma gamma
# ==========================================================================
alg_name_4 = "Lc2SigmaPi0"
alg_4 = Algorithm.new(alg_name_4)
alg_4.set_header(["#{alg_name_4}Alg/#{alg_name_4}.h"])
     .set_constant({"ECMS" => [:double, 4.6]})
     .note(:ebeam_per_run, "Per-run beam energy required for the 4C fit (seven energy points).")

sel_4 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp    ">=1"     # p only
    nChrn    ">=1"
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    nGam             ">=4"
  }
  .pid(method: :probability) {
    prob_cut 0.0
    identify :proton, against: [:kaon, :pion]
    identify :kaon,   against: [:pion]
    identify :pion,   against: [:kaon]
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # reconstruct the two pi0 from gamma gamma pairs
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 100
    npi0 ">=2"
  }
  .kinematic_fit([:prp, :pi0, :pi0]) {          # 4C fit on p pi0 + pi0
    nominal
    constrain_four_momentum
    invariant_mass_of(:prp, :pi0).constrain_to_nominal_mass_of(:Sigma_p)
    invariant_mass_of(:prp, :pi0).within(1.174, 1.200)   # Sigma+ window
    invariant_mass_of(:pi0, :pi0).out_of(0.440, 0.520)   # p K_S0 veto
    chi2_cut 17
  }

alg_4.with_decay_card(decay_card_lc_sigma_pi0).apply(sel_4)

# ==========================================================================
# Reference mode 5: Lambda_c+ -> Sigma+ omega, omega -> pi+ pi- pi0
# ==========================================================================
alg_name_5 = "Lc2SigmaOmega"
alg_5 = Algorithm.new(alg_name_5)
alg_5.set_header(["#{alg_name_5}Alg/#{alg_name_5}.h"])
     .set_constant({"ECMS" => [:double, 4.6]})
     .note(:ebeam_per_run, "Per-run beam energy required for the 4C fit (seven energy points).")

sel_5 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp    ">=2"     # p and pi+ (omega -> pi+ pi- pi0)
    nChrn    ">=2"     # pi-
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    nGam             ">=4"
  }
  .pid(method: :probability) {
    prob_cut 0.0
    identify :proton, against: [:kaon, :pion]
    identify :kaon,   against: [:pion]
    identify :pion,   against: [:kaon]
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # reconstruct the two pi0 (Sigma+ and omega daughters)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 100
    npi0 ">=2"
  }
  .kinematic_fit([:prp, :pi0, :pip, :pim, :pi0]) {  # 4C fit on p pi0 + omega(-> pi+ pi- pi0)
    nominal
    constrain_four_momentum
    invariant_mass_of(:prp, :pi0).constrain_to_nominal_mass_of(:Sigma_p)
    invariant_mass_of(:prp, :pi0).within(1.174, 1.200)   # Sigma+ window
    invariant_mass_of(:pi0, :pi0).out_of(0.440, 0.520)   # p K_S0 veto
    chi2_cut 17
  }

alg_5.with_decay_card(decay_card_lc_sigma_omega).apply(sel_5)

### Execute on the seven data points, matching inclusive MC and the signal/reference exclusive MC ###
root_files_1 = alg_1.execute_on(data_points + incMC_points + exMC_lc_sigma_eta_gg)
root_files_2 = alg_2.execute_on(data_points + incMC_points + exMC_lc_sigma_eta_3pi)
root_files_3 = alg_3.execute_on(data_points + incMC_points + exMC_lc_sigma_etap)
root_files_4 = alg_4.execute_on(data_points + incMC_points + exMC_lc_sigma_pi0)
root_files_5 = alg_5.execute_on(data_points + incMC_points + exMC_lc_sigma_omega)