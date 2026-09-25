# ============================================================
# Dataset preparation
# ============================================================
jpsi_data   = DatasetManager.real_data.find("708_3097")      # J/psi(3097) real data
jpsi_incMC  = DatasetManager.inclusive_mc.find("708_3097")   # J/psi(3097) inclusive MC
cont_3080   = DatasetManager.real_data.find("708_3080")      # 3.080 GeV continuum data (background)

# Decay card: signal  J/psi -> gamma pi0 pi0  (phase space)
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Decay card: background  J/psi -> gamma eta, eta -> pi0 pi0 pi0
decay_card_eta = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta PHSP;
  Enddecay

  Decay eta
  1.0000 pi0 pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Decay card: background  J/psi -> gamma eta', eta' -> eta pi0 pi0, eta -> gamma gamma
decay_card_etap = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.0000 eta pi0 pi0 PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# One 100k-event exclusive MC per mode
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammapi0pi0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

exMC_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaeta"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_eta
  config.cross_section   = :default
end

exMC_etap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaetap"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_etap
  config.cross_section   = :default
end

# ============================================================
# Event selection (BOSS)
# ============================================================
alg_name = "GammaPi0Pi0"
jpsi_alg = Algorithm.new(alg_name)
jpsi_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})   # J/psi center-of-mass energy (GeV)

event_selection = Selection.new
  .select_track {              # charged-track veto: reject any event with a good charged track
    cos_theta 0.93             # |cos(theta)| < 0.93
    Vz        10.0             # |Vz| < 10 cm
    Vr        1.0              # Vr < 1 cm
    nChrp     "==0"            # zero positive charged tracks
    nChrn     "==0"            # zero negative charged tracks
    nNet      "==0"            # zero net charge
  }
  # no PID is applied on top of the veto
  .select_photon {             # photon selection
    tdc_emc_start     0        # EMC timing window start
    tdc_emc_end       14       # EMC timing window end
    angle_to_track    10.0     # at least 10 degrees from any charged track
    energyThreshold_b 0.025    # E > 25 MeV in the barrel
    energyThreshold_e 0.050    # E > 50 MeV in the endcap
    nGam              ">=5"    # at least five good photons
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pair photons -> pi0 (1-C Kalman mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25                               # chi2 < 25 for the pi0 reconstruction
    npi0     ">=2"                            # require at least two pi0 candidates
  }
  .kinematic_fit([:gamma, :pi0, :pi0]) {      # nominal fit of gamma pi0 pi0  -> total 6C (4C + 2 pi0 masses)
    nominal                                   # nominal fit: its four-momenta are used downstream
    constrain_four_momentum                   # constrain total four-momentum to the J/psi four-momentum
    chi2_cut 200                              # loose BOSS-pass cut; the published mass-dependent chi2 cut is applied later
  }                                           # min-chi2 photon/pi0 combination is chosen automatically

# Attach the signal decay card and render the selection into the BOSS algorithm
jpsi_alg.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC, continuum background, and the three exclusive-MC samples
root_files = jpsi_alg.execute_on([jpsi_data, jpsi_incMC, cont_3080, exMC_signal, exMC_eta, exMC_etap])