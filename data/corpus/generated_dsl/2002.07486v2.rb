# ============================================================
# Dataset description — J/psi @ 3.097 GeV (BOSS 708)
# ============================================================
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi real data (BOSS 708)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC

# ------------------------------------------------------------
# Decay cards (EvtGen syntax) — one per signal mode
# ------------------------------------------------------------
# Mode 1: J/psi -> gamma pi0 eta', eta' -> pi+ pi- eta, eta -> gamma gamma, pi0 -> gamma gamma
decay_card_gamma_pi0_etap = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma pi0 eta' PHSP;
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

# Mode 2: J/psi -> omega eta', omega -> gamma pi0, eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_omega_etap = <<~DECAYCARD
    Decay J/psi
    1.0000 omega eta' PHSP;
    Enddecay

    Decay omega
    1.0000 gamma pi0 PHSP;
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

# Mode 3: J/psi -> gamma eta_c, eta_c -> pi0 eta', eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_gamma_etac = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0000 pi0 eta' PHSP;
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

# ------------------------------------------------------------
# Exclusive MC samples — 200k events each
# ------------------------------------------------------------
exMC_gamma_pi0_etap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_pi0_etap"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_gamma_pi0_etap
  config.cross_section   = :default
end

exMC_omega_etap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_omega_etap"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_omega_etap
  config.cross_section   = :default
end

exMC_gamma_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etac"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_gamma_etac
  config.cross_section   = :default
end

# ============================================================
# Common event selection (shared by all three signal modes)
# ============================================================
common_selection = Selection.new
  .select_track {            # Charged track selection
    cos_theta 0.93           # |cos(theta)| < 0.93
    Vz        10.0           # |Vz| < 10 cm
    Vr        1.0            # Vr < 1 cm
    nChrp     ">=1"          # At least one positive track
    nChrn     ">=1"          # At least one negative track
    nNet      "==0"          # Net charge zero
  }
  .select_photon {           # Photon selection
    tdc_emc_start     0      # EMC TDC window start
    tdc_emc_end       14     # EMC TDC window end
    energyThreshold_b 0.025  # Barrel energy threshold 25 MeV
    energyThreshold_e 0.050  # Endcap energy threshold 50 MeV
    angle_to_track    10.0   # More than 10 degrees from any charged track
    nGam              ">=5"  # At least five photons
  }
  .pid(method: :probability) {   # PID by the probability method
    prob_cut 0.001
    identify :pion, against: [:kaon]   # pi/K separation (pi+ and pi-)
    npip ">=1"                          # At least one pi+
    npim ">=1"                          # At least one pi-
  }
  # Nominal 4C kinematic fit to pi+ pi- 5gamma (four-momentum constrained, chi2 < 200)
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  # Competing pi+ pi- 6gamma hypothesis (no chi2_cut, no nominal):
  # its chi2 is stored and the comparison with the nominal chi2 is applied later in ROOT
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }
  # Mass-constrained Kalman fit: rebuild pi0 from gamma gamma
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"                          # At least one pi0 candidate
  }
  # Mass-constrained Kalman fit: rebuild eta from gamma gamma
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta ">=1"                          # At least one eta candidate
  }
  # 5C fit to pi+ pi- 3gamma eta, eta mass constrained (four-momentum + eta mass)
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :eta]) {
    constrain_four_momentum
    chi2_cut 30
  }

# ============================================================
# Algorithm 1 — J/psi -> gamma pi0 eta'
# ============================================================
alg_name_1 = "JpsiGammaPi0EtaP"
alg_mode1 = Algorithm.new(alg_name_1)
alg_mode1.set_header(["#{alg_name_1}Alg/#{alg_name_1}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})
alg_mode1.with_decay_card(decay_card_gamma_pi0_etap).apply(common_selection.dup)
alg_mode1.execute_on([jpsi_data, jpsi_incMC, exMC_gamma_pi0_etap])

# ============================================================
# Algorithm 2 — J/psi -> omega eta' (omega -> gamma pi0)
# ============================================================
alg_name_2 = "JpsiOmegaEtaP"
alg_mode2 = Algorithm.new(alg_name_2)
alg_mode2.set_header(["#{alg_name_2}Alg/#{alg_name_2}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})
alg_mode2.with_decay_card(decay_card_omega_etap).apply(common_selection.dup)
alg_mode2.execute_on([jpsi_data, jpsi_incMC, exMC_omega_etap])

# ============================================================
# Algorithm 3 — J/psi -> gamma eta_c (eta_c -> pi0 eta')
# ============================================================
alg_name_3 = "JpsiGammaEtaC"
alg_mode3 = Algorithm.new(alg_name_3)
alg_mode3.set_header(["#{alg_name_3}Alg/#{alg_name_3}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})
alg_mode3.with_decay_card(decay_card_gamma_etac).apply(common_selection.dup)
alg_mode3.execute_on([jpsi_data, jpsi_incMC, exMC_gamma_etac])