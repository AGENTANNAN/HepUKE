# ================================================================
# BESIII : search for the lepton-number-violating decay
#          omega -> pi+ pi- e+ e-  in  J/psi -> omega eta, eta -> gamma gamma
#          at the J/psi peak (sqrt(s) = 3.097 GeV)
#   * reference mode : omega -> pi+ pi- pi0 (pi0 -> gamma gamma)
#   * signal mode    : omega -> pi+ pi- e+ e-
#   * both channels  : eta -> gamma gamma
# ================================================================

### ----------------------------- Datasets ----------------------------- ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi-peak real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # matching inclusive MC

### --------------------------- Decay cards --------------------------- ###
# Reference mode : J/psi -> omega eta, omega -> pi+ pi- pi0, pi0 -> g g, eta -> g g
decay_card_ref = <<~DECAYCARD
    Decay J/psi
    1.0000 omega eta PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Signal mode : J/psi -> omega eta, omega -> pi+ pi- e+ e-, eta -> g g
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 omega eta PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- e+ e- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### --------------------- Exclusive MC (100k events) --------------------- ###
exMC_ref = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_708_3097_omega_eta_ref"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_ref
  config.cross_section   = :default
end

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_708_3097_omega_eta_lnv"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### ================= Reference channel : pi+ pi- 4 gamma ================= ###
alg_ref_name = "JpsiOmegaEtaRef"
alg_ref = Algorithm.new(alg_ref_name)
alg_ref.set_header(["#{alg_ref_name}Alg/#{alg_ref_name}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})

ref_selection = Selection.new
  # exactly two charged tracks, net charge 0
  .select_track {
    cos_theta 0.93   # |cos(theta)| < 0.93
    Vz        10.0   # |Vz| < 10 cm
    Vr        1.0    # Vr < 1 cm
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  # at least four photons, E > 25 MeV, angle to nearest charged track > 10 deg
  .select_photon {
    tdc_emc_start      0
    tdc_emc_end        14
    angle_to_track     10.0
    energyThreshold_b  0.025
    energyThreshold_e  0.050
    nGam               ">=4"
  }
  # pion identification (CL > 0.001 against kaons and protons)
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  }
  # pi0 -> gamma gamma, 1C Kalman fit, chi2 < 200
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  }
  # eta -> gamma gamma, 1C Kalman fit, chi2 < 200
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta ">=1"
  }
  # nominal 5C fit : four-momentum + pi0 mass constraint, chi2 < 20
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
  }
  # competing 4C hypotheses (chi2 stored only; discrimination performed in ROOT)
  .kinematic_fit([:pip, :pim, :pip, :pim, :gamma, :gamma]) {
    constrain_four_momentum
  }
  .assign({:pip => :kp, :pim => :km})
  .kinematic_fit([:kp, :km, :kp, :km, :gamma, :gamma]) {
    constrain_four_momentum
  }

alg_ref.with_decay_card(decay_card_ref).apply(ref_selection)

### ================= Signal channel : e+e- pi+ pi- 2 gamma ================= ###
alg_sig_name = "JpsiOmegaEtaLNV"
alg_sig = Algorithm.new(alg_sig_name)
alg_sig.set_header(["#{alg_sig_name}Alg/#{alg_sig_name}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})
       # The description uses a CL-based electron selection
       # (CL_e > 0.001 and CL_e/(CL_e+CL_K+CL_pi) > 0.8). The DSL exposes only
       # identify_high_momentum_leptons (momentum / EMC-energy recipe), so the
       # residual difference is carried as a PID systematic.
       .note(:pid_correction_method,
             "signal channel electron identification uses confidence-level cuts " \
             "CL_e > 0.001 and CL_e/(CL_e+CL_K+CL_pi) > 0.8; replaced here by " \
             "identify_high_momentum_leptons (momentum / EMC-energy recipe) - " \
             "difference treated as a PID systematic")

sig_selection = Selection.new
  # four charged tracks (e+ e- pi+ pi-), net charge 0
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  # at least two photons above 25 MeV
  .select_photon {
    tdc_emc_start      0
    tdc_emc_end        14
    angle_to_track     10.0
    energyThreshold_b  0.025
    energyThreshold_e  0.050
    nGam               ">=2"
  }
  # electron identification together with pion identification
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon, :proton]
    nlp  "==1"
    nlm  "==1"
    npip "==1"
    npim "==1"
  }
  # eta -> gamma gamma, 1C Kalman fit, chi2 < 200
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta ">=1"
  }
  # nominal 4C fit (no eta mass constraint in the fit), chi2 < 10
  .kinematic_fit([:lp, :lm, :pip, :pim, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 10
  }
  # competing 4C hypotheses (chi2 stored only; discrimination performed in ROOT)
  #   pi+ pi- pi+ pi- gamma gamma  (e+e- reinterpreted as pions)
  .assign({:lp => :pip, :lm => :pim})
  .kinematic_fit([:pip, :pim, :pip, :pim, :gamma, :gamma]) {
    constrain_four_momentum
  }
  #   K+ K- K+ K- gamma gamma      (all four charged tracks as kaons)
  .assign({:pip => :kp, :pim => :km, :lp => :kp, :lm => :km})
  .kinematic_fit([:kp, :km, :kp, :km, :gamma, :gamma]) {
    constrain_four_momentum
  }
  #   pi+ pi- K+ K- gamma gamma    (e+e- reinterpreted as kaons)
  .assign({:lp => :kp, :lm => :km})
  .kinematic_fit([:pip, :pim, :kp, :km, :gamma, :gamma]) {
    constrain_four_momentum
  }
  #   pi+ pi- p pbar gamma gamma   (e+e- reinterpreted as protons)
  .assign({:lp => :prp, :lm => :prm})
  .kinematic_fit([:pip, :pim, :prp, :prm, :gamma, :gamma]) {
    constrain_four_momentum
  }

alg_sig.with_decay_card(decay_card_signal).apply(sig_selection)

### ----------------------------- Execute ----------------------------- ###
alg_ref.execute_on([jpsi_data, jpsi_incMC, exMC_ref, exMC_signal])
alg_sig.execute_on([jpsi_data, jpsi_incMC, exMC_ref, exMC_signal])