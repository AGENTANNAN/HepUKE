# ============================================================================
# BESIII R-scan analysis of e+e- -> omega eta and e+e- -> omega pi0
# (omega -> pi+ pi- pi0, pi0 / eta -> gamma gamma)
# over the 22 scan points 2.000 - 3.080 GeV (BOSS 713), pi+pi-4gamma final state
# ============================================================================

### Datasets ###
# 22 R-scan energy points, 2.000 ... 3.080 GeV (2015 R scan, BOSS 713)
scan_sample_names = %w[
  713_Rscan_2000 713_Rscan_2050 713_Rscan_2100 713_Rscan_2125
  713_Rscan_2150 713_Rscan_2175 713_Rscan_2200 713_Rscan_2232
  713_Rscan_2309 713_Rscan_2386 713_Rscan_2396 713_Rscan_2500
  713_Rscan_2644 713_Rscan_2646 713_Rscan_2700 713_Rscan_2800
  713_Rscan_2900 713_Rscan_2950 713_Rscan_2981 713_Rscan_3000
  713_Rscan_3020 713_Rscan_3080
]
scan_data  = scan_sample_names.map { |n| DatasetManager.real_data.find(n) }
scan_incMC = scan_sample_names.map { |n| DatasetManager.inclusive_mc.find(n) }
data_2125  = DatasetManager.real_data.find("713_Rscan_2125")  # reference point for the signal MC

### Decay cards ###
# Mode I : e+e- -> omega eta ; omega -> pi+ pi- pi0 ; eta -> gamma gamma ; pi0 -> gamma gamma
decay_card_omega_eta = <<~DECAYCARD
  Decay psi(4260)
  1.000 omega eta PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: e+e- -> omega pi0 ; omega -> pi+ pi- pi0 ; both pi0 -> gamma gamma
decay_card_omega_pi0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 omega pi0 PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC (200k events per mode, both matched to the 2.125 GeV data set) ###
exMC_omega_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_omega_eta"
  config.related_dataset = data_2125
  config.events          = 200_000
  config.decay_card      = decay_card_omega_eta
  config.cross_section   = :default
end

exMC_omega_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_omega_pi0"
  config.related_dataset = data_2125
  config.events          = 200_000
  config.decay_card      = decay_card_omega_pi0
  config.cross_section   = :default
end

### Algorithms (one per signal mode - Rule T1) ###
alg_name_omega_eta = "OmegaEta"
alg_omega_eta = Algorithm.new(alg_name_omega_eta)
alg_omega_eta.set_header(["#{alg_name_omega_eta}Alg/#{alg_name_omega_eta}.h"])
             .set_constant({"ECMS" => [:double, 2.125]})

alg_name_omega_pi0 = "OmegaPi0"
alg_omega_pi0 = Algorithm.new(alg_name_omega_pi0)
alg_omega_pi0.set_header(["#{alg_name_omega_pi0}Alg/#{alg_name_omega_pi0}.h"])
             .set_constant({"ECMS" => [:double, 2.125]})

### Common preselection (charged tracks, photons, PID - shared by both modes) ###
event_selection_common = Selection.new
  .select_track {
    cos_theta 0.93      # |cos(theta)| < 0.93
    Vr        1.0       # Vr < 1 cm
    Vz        10.0      # |Vz| < 10 cm
    nChrp     "==1"     # exactly one positive charged track
    nChrn     "==1"     # exactly one negative charged track
    nNet      "==0"     # net charge zero
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14        # EMC time window [0, 700] ns
    angle_to_track    10.0      # photon-track opening angle > 10 degrees
    energyThreshold_b 0.025     # E > 25 MeV in the barrel (|cos(theta)| < 0.80)
    energyThreshold_e 0.050     # E > 50 MeV in the endcap (0.86 < |cos(theta)| < 0.92)
    nGam              ">=4"     # at least four photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]   # pi+ and pi- (charge-conjugation shorthand)
    npip "==1"
    npim "==1"
  }

### Mode I : omega eta ###
omega_eta_selection = event_selection_common.dup
  # pi0 from a gamma-gamma pair, Kalman mass-constrained fit (chi2 < 25)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # eta from a gamma-gamma pair, Kalman mass-constrained fit (chi2 < 25)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  # Nominal 4C fit to pi+ pi- pi0 eta (loose chi2; the optimal tight cut is applied in ROOT)
  .kinematic_fit([:pip, :pim, :pi0, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 70
  }
  # Competing hypothesis pi+ pi- 5gamma (one extra photon); chi2 stored for the
  # ROOT-level veto chi2(pi+pi-4gamma) < chi2(pi+pi-5gamma)
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }

alg_omega_eta
  .note(:combinatorial_selection, "the pi0-eta pairing is chosen to minimise the combined pi0/eta mass chi2 and is required to be smaller than the pi0pi0 and etaeta alternative pairings")
  .note(:background_veto, "veto |E(gamma3)-E(gamma4)|/p(eta) < 0.9 against e+e- -> omega gamma_ISR, applied together with the chi2(pi+pi-4gamma) < chi2(pi+pi-5gamma) veto in ROOT")
  .with_decay_card(decay_card_omega_eta)
  .apply(omega_eta_selection)

alg_omega_eta.execute_on(scan_data + scan_incMC + [exMC_omega_eta])

### Mode II : omega pi0 ###
omega_pi0_selection = event_selection_common.dup
  # two pi0 candidates, each from a gamma-gamma pair, Kalman mass-constrained fit (chi2 < 25)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  # Nominal 4C fit to pi+ pi- pi0 (loose chi2; the optimal tight cut is applied in ROOT)
  .kinematic_fit([:pip, :pim, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 70
  }

alg_omega_pi0
  .note(:combinatorial_selection, "the pi0pi0 pairing is chosen by smallest Kalman chi2 and the pi0 whose pi+pi-pi0 mass is closest to m_omega is taken as the omega daughter")
  .with_decay_card(decay_card_omega_pi0)
  .apply(omega_pi0_selection)

alg_omega_pi0.execute_on(scan_data + scan_incMC + [exMC_omega_pi0])