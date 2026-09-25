# =============================================================================
#  e+e- -> gamma chi_c1(3872), chi_c1(3872) -> pi+ pi- eta
#    Mode I : eta -> gamma gamma
#    Mode II: eta -> pi+ pi- pi0, pi0 -> gamma gamma
#  Data: 4.13 - 4.34 GeV scan (705 real-data points) + corresponding inclusive MC
# =============================================================================

### ------------------------- Dataset description ------------------------- ###
# Real-data points of the 4.13-4.34 GeV scan (BOSS 7.0.5, sample name = boss_energy)
data_points = [
  DatasetManager.real_data.find("705_4130"),   # 4128.48 MeV
  DatasetManager.real_data.find("705_4160"),   # 4157.44 MeV
  DatasetManager.real_data.find("705_4290"),   # 4287.88 MeV
  DatasetManager.real_data.find("705_4315"),   # 4312.05 MeV
  DatasetManager.real_data.find("705_4340"),   # 4337.39 MeV
]

# Corresponding inclusive MC samples (one per energy point)
incMC_points = [
  DatasetManager.inclusive_mc.find("705_4130"),
  DatasetManager.inclusive_mc.find("705_4160"),
  DatasetManager.inclusive_mc.find("705_4290"),
  DatasetManager.inclusive_mc.find("705_4315"),
  DatasetManager.inclusive_mc.find("705_4340"),
]

### ---------------------------- Decay cards ----------------------------- ###
# Mode I: chi_c1(3872) -> pi+ pi- eta, eta -> gamma gamma
decay_card_mode_I = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma chi_c1(3872) PHSP;
  Enddecay

  Decay chi_c1(3872)
  1.0000 pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: chi_c1(3872) -> pi+ pi- eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_mode_II = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma chi_c1(3872) PHSP;
  Enddecay

  Decay chi_c1(3872)
  1.0000 pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.0000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### --------------------------- Exclusive MC ----------------------------- ###
# 100k events of signal MC generated for each mode at every energy point.
exMC_mode_I = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_chic1_eta_gammagamma"  # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_mode_I
  config.cross_section = :default
end

exMC_mode_II = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_chic1_eta_3pi"         # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_mode_II
  config.cross_section = :default
end

# =============================================================================
#  Mode I :  e+e- -> gamma chi_c1(3872), chi_c1(3872) -> pi+ pi- eta, eta -> gamma gamma
# =============================================================================
alg_modeI = Algorithm.new("Chic1EtaGamGam")
alg_modeI.set_header(["Chic1EtaGamGamAlg/Chic1EtaGamGam.h"])
         .set_constant({"ECMS" => [:double, 4.26]})   # nominal generator energy

selection_mode_I = Selection.new
selection_mode_I
  .select_track {                # charged-track quality cuts + exactly 1 pi+ / 1 pi-
    cos_theta 0.93               # |cos(theta)| < 0.93
    Vz        10.0               # |Vz| < 10 cm
    Vr        1.0                # Vr < 1 cm
    nChrp     "==1"              # exactly one positive track
    nChrn     "==1"              # exactly one negative track
    nNet      "==0"              # net charge zero
  }
  .select_photon {               # good photon selection, at least 3 photons
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0       # >= 10 deg from any charged track
    energyThreshold_b 0.025      # 25 MeV in the barrel
    energyThreshold_e 0.050      # 50 MeV in the endcap
    nGam              ">=3"      # gamma_rad + gamma_1 + gamma_2
  }
  .assign({:chrgp => :pip, :chrgn => :pim})   # the single (pi+, pi-) pair is taken as pi+pi-
  # 5C fit: gamma_rad pi+ pi- gamma gamma, with M(gamma gamma) constrained to m(eta)
  # and four-momentum conservation (4C + 1 mass = 5C)
  .kinematic_fit([:gamma, :pip, :pim, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  }

# Background vetoes that are applied on fitted quantities / explicit photon-pair
# combinations and therefore cannot be expressed by the type-based participant lists.
alg_modeI
  .note(:background_veto,
    "reject pi0 -> gamma gamma mis-combinations when |M(gamma_rad gamma_1(2)) - m(pi0)| " \
    "< 20 MeV/c^2; each radiative photon combined with either eta photon is checked " \
    "after the 5C fit (evaluated in ROOT on kinematically fitted momenta)")
  .note(:radiative_dimuon_veto,
    "suppress the radiative dimuon background by requiring cos(theta_pipi) < -0.96, " \
    "where theta_pipi is the opening angle between the two pions")
  .with_decay_card(decay_card_mode_I)
  .apply(selection_mode_I)

# =============================================================================
#  Mode II : e+e- -> gamma chi_c1(3872), chi_c1(3872) -> pi+ pi- eta,
#            eta -> pi+ pi- pi0, pi0 -> gamma gamma
# =============================================================================
alg_modeII = Algorithm.new("Chic1EtaThreePi")
alg_modeII.set_header(["Chic1EtaThreePiAlg/Chic1EtaThreePi.h"])
          .set_constant({"ECMS" => [:double, 4.26]})   # nominal generator energy

selection_mode_II = Selection.new
selection_mode_II
  .select_track {                # charged-track quality cuts + exactly 2 pi+ / 2 pi-
    cos_theta 0.93               # |cos(theta)| < 0.93
    Vz        10.0               # |Vz| < 10 cm
    Vr        1.0                # Vr < 1 cm
    nChrp     "==2"              # exactly two positive tracks
    nChrn     "==2"              # exactly two negative tracks
    nNet      "==0"              # net charge zero
  }
  .select_photon {               # good photon selection, at least 3 photons
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0       # >= 10 deg from any charged track
    energyThreshold_b 0.025      # 25 MeV in the barrel
    energyThreshold_e 0.050      # 50 MeV in the endcap
    nGam              ">=3"      # gamma_rad + gamma_1 + gamma_2
  }
  .pid(method: :probability) {   # probability PID: pions against kaons
    prob_cut 0.001
    identify :pion, against: [:kaon]   # pi+ and pi- (charge-conjugation shorthand)
    npip ">=1"                         # at least one pi+
    npim ">=1"                         # at least one pi-
  }
  # intermediate pi0 -> gamma gamma with gamma-gamma mass constraint (1C Kalman fit)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"                         # at least one pi0 candidate
  }
  # 5C fit: pi+ pi- pi+ pi- pi0 gamma, with M(pi+pi-pi+pi-pi0) constrained to m(eta)
  # and four-momentum conservation (4C + 1 mass = 5C)
  .kinematic_fit([:pip, :pim, :pip, :pim, :pi0, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:pip, :pim, :pip, :pim, :pi0).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  }

# Veto against e+e- -> pi+ pi- pi+ pi- pi0, applied on explicit gamma_rad / eta-photon
# combinations after the fit - not expressible through the type-based participant lists.
alg_modeII
  .note(:background_veto,
    "veto e+e- -> pi+ pi- pi+ pi- pi0 by rejecting events with " \
    "M(gamma_rad gamma_1) in [125.6, 150.0] MeV/c^2 or " \
    "M(gamma_rad gamma_2) in [115.7, 160.0] MeV/c^2; " \
    "gamma_1/gamma_2 are the two eta photons, evaluated on fitted momenta in ROOT")
  .with_decay_card(decay_card_mode_II)
  .apply(selection_mode_II)

### ------------------------------ Execute --------------------------------- ###
datasets = data_points + incMC_points

root_files_mode_I  = alg_modeI.execute_on(datasets + exMC_mode_I)    # data + incMC + mode-I signal MC
root_files_mode_II = alg_modeII.execute_on(datasets + exMC_mode_II)  # data + incMC + mode-II signal MC