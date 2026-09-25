# =============================================================================
# J/psi -> gamma eta eta eta' at sqrt(s) = 3.097 GeV
#   Mode I : eta' -> gamma pi+ pi-
#   Mode II: eta' -> pi+ pi- eta
# (eta -> gamma gamma in both modes)
# BOSS part: dataset preparation + event selection up to the kinematic fit.
# The pi0 veto, extra-eta veto and, for Mode I, the M(pi+pi-) > 0.5 GeV/c^2
# requirement are post-fit ROOT-level selections and are therefore not encoded here.
# =============================================================================

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi real data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding inclusive MC

# ---- Decay cards (EvtGen format) ----
# Mode I: eta' -> gamma pi+ pi-
decay_card_modeI = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta eta eta' PHSP;
  Enddecay

  Decay eta'
  1.0000 gamma pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: eta' -> pi+ pi- eta
decay_card_modeII = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta eta eta' PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- Exclusive MC samples (500,000 events each) ----
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_g_eta_eta_etap_modeI"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_g_eta_eta_etap_modeII"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Two independent decay modes with different final states, photon multiplicities and
# kinematic-fit hypotheses -> separate Algorithm objects (Rule T1).
alg_modeI = Algorithm.new("JpsiEtaEtaEtapI")
alg_modeI.set_header(["JpsiEtaEtaEtapIAlg/JpsiEtaEtaEtapI.h"])
        .set_constant({"ECMS" => [:double, 3.097]})   # sqrt(s) = 3.097 GeV

alg_modeII = Algorithm.new("JpsiEtaEtaEtapII")
alg_modeII.set_header(["JpsiEtaEtaEtapIIAlg/JpsiEtaEtaEtapII.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

# Common charged-track selection shared by both modes (no PID is applied in this analysis)
common_track = Selection.new.select_track {
  cos_theta 0.93    # |cos(theta)| < 0.93
  Vz        10.0    # |Vz| < 10 cm
  Vr        1.0     # Vr < 1 cm
  nChrp     ">=1"   # at least one positively charged track
  nChrn     ">=1"   # at least one negatively charged track
}

# ---- Mode I: eta' -> gamma pi+ pi-  (>= 6 photons) ----
sel_modeI = common_track.dup
  .select_photon {
    tdc_emc_start      0        # TDC window 0-14
    tdc_emc_end        14
    angle_to_track     10.0     # at least 10 degrees from any charged track
    energyThreshold_b  0.025    # barrel energy > 25 MeV
    energyThreshold_e  0.050    # endcap energy > 50 MeV
    nGam ">=6"                  # gamma_rad + gamma(eta') + 2 eta -> 4 gamma = 6
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct eta from gamma gamma pairs
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25                               # chi2 < 25
    neta ">=2"                                # Mode I contains two eta mesons
  }
  .assign({:chrgp => :pip, :chrgn => :pim})   # no PID: all charged tracks treated as pions
  # Main fit: 4C (four-momentum) + two eta mass constraints -> 6C; flagged nominal
  .kinematic_fit([:gamma, :gamma, :eta, :eta, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 30                               # chi2 < 30
  }

# ---- Mode II: eta' -> pi+ pi- eta  (>= 7 photons) ----
sel_modeII = common_track.dup
  .select_photon {
    tdc_emc_start      0        # TDC window 0-14
    tdc_emc_end        14
    angle_to_track     10.0     # at least 10 degrees from any charged track
    energyThreshold_b  0.025    # barrel energy > 25 MeV
    energyThreshold_e  0.050    # endcap energy > 50 MeV
    nGam ">=7"                  # gamma_rad + 3 eta -> 6 gamma = 7
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct eta from gamma gamma pairs
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25                               # chi2 < 25
    neta ">=3"                                # Mode II contains three eta mesons
  }
  .assign({:chrgp => :pip, :chrgn => :pim})   # no PID: all charged tracks treated as pions
  # Main fit: 4C (four-momentum) + three eta mass constraints -> 7C; flagged nominal
  .kinematic_fit([:gamma, :eta, :eta, :eta, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 50                               # chi2 < 50
  }

# Attach the decay card and render each algorithm's event selection
alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

# Execute on real data, inclusive MC and the corresponding exclusive MC
root_files_modeI  = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])
root_files_modeII = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])