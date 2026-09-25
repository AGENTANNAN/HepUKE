# =============================================================================
# BESIII observation of J/psi -> phi pi0, first evidence for a doubly OZI
# suppressed electromagnetic J/psi decay  (arXiv:1504.03194v2)
# -- BOSS part (dataset preparation + event selection)
#
# Signal: J/psi -> phi pi0 -> K+ K- gamma gamma, with 1.311 x 10^9 J/psi
# events collected with the BESIII detector.
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets (single energy: the J/psi resonance)
### ---------------------------------------------------------------------------
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi, 3097 MeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # inclusive J/psi MC

### ---------------------------------------------------------------------------
### Decay cards
### ---------------------------------------------------------------------------
# Signal process: J/psi -> phi pi0 -> K+ K- gamma gamma
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 phi pi0 PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Dominant coherent background: J/psi -> K+ K- pi0 (intermediate states decaying
# into K+- pi0 and K+ K-).  This is the process that interferes with phi pi0.
decay_card_bkg_kkm_pi0 = <<~DECAYCARD
  Decay J/psi
  1.0000 K+ K- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# phi-peaking background: J/psi -> phi pi0 pi0
decay_card_bkg_phi_pi0pi0 = <<~DECAYCARD
  Decay J/psi
  1.0000 phi pi0 pi0 PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# phi-peaking background: J/psi -> phi gamma gamma
decay_card_bkg_phi_gg = <<~DECAYCARD
  Decay J/psi
  1.0000 phi gamma gamma PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  End
DECAYCARD

# pi0-peaking background: J/psi -> gamma eta_c(1S) -> gamma K+ K- pi0
decay_card_bkg_etac = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta_c(1S) PHSP;
  Enddecay

  Decay eta_c(1S)
  1.000 K+ K- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Non-phi / non-pi0 background: J/psi -> gamma K+ K-
decay_card_bkg_gkk = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma K+ K- PHSP;
  Enddecay

  End
DECAYCARD

# Non-phi / non-pi0 background: J/psi -> gamma pi0 K+ K-
decay_card_bkg_gpi0kk = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma pi0 K+ K- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Non-phi / non-pi0 background: J/psi -> pi0 pi0 K+ K-
decay_card_bkg_pi0pi0kk = <<~DECAYCARD
  Decay J/psi
  1.0000 pi0 pi0 K+ K- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Additional peaking background used for the K+K- mass resolution study:
# J/psi -> phi eta with eta -> gamma gamma
decay_card_bkg_phi_eta = <<~DECAYCARD
  Decay J/psi
  1.0000 phi eta PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### ---------------------------------------------------------------------------
### Exclusive MC
### ---------------------------------------------------------------------------
# Signal MC for J/psi -> phi pi0 -> K+ K- gamma gamma
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_phi_pi0"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Exclusive background MC samples (equivalent to / larger than the data
# statistics).  The inclusive J/psi sample (jpsi_incMC) covers the remaining
# known and unknown decay modes, generated with BesEvtGen + Lund-Charm.
exMC_bkg_kkm_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "bkg_jpsi_K+K-pi0"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_bkg_kkm_pi0
  config.cross_section   = :default
end

exMC_bkg_phi_pi0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "bkg_jpsi_phi_pi0pi0"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_phi_pi0pi0
  config.cross_section   = :default
end

exMC_bkg_phi_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "bkg_jpsi_phi_gg"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_phi_gg
  config.cross_section   = :default
end

exMC_bkg_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "bkg_jpsi_gamma_etac"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_etac
  config.cross_section   = :default
end

exMC_bkg_gkk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "bkg_jpsi_gamma_KK"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_gkk
  config.cross_section   = :default
end

exMC_bkg_gpi0kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "bkg_jpsi_gamma_pi0_KK"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_gpi0kk
  config.cross_section   = :default
end

exMC_bkg_pi0pi0kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "bkg_jpsi_pi0pi0_KK"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_pi0pi0kk
  config.cross_section   = :default
end

exMC_bkg_phi_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "bkg_jpsi_phi_eta"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_phi_eta
  config.cross_section   = :default
end

### ---------------------------------------------------------------------------
### Algorithm + event selection
### ---------------------------------------------------------------------------
alg_phipi0 = Algorithm.new("JpsiToPhiPi0")
alg_phipi0.set_header(["JpsiToPhiPi0Alg/JpsiToPhiPi0.h"])
          .set_constant({ "ECMS" => [:double, 3.097] })   # J/psi c.m. energy

event_selection = Selection.new

# --- Charged track selection -------------------------------------------------
# Two charged tracks with opposite charge.  Each track must lie inside the MDC
# acceptance |cos(theta)| < 0.93, and its point of closest approach to the
# e+e- interaction point must be within +/-10 cm along the beam direction and
# within 1 cm in the plane perpendicular to the beam.
event_selection.select_track do
  cos_theta 0.93
  Vz        10.0
  Vr        1.0
  nChrp     "==1"   # one positive track  (K+ candidate)
  nChrn     "==1"   # one negative track  (K- candidate)
  nNet      "==0"   # opposite charges
end

# --- Photon selection --------------------------------------------------------
# At least two photons.  E_EMC >= 25 MeV for |cos(theta)| < 0.8 (barrel) and
# >= 50 MeV for 0.86 < |cos(theta)| < 0.92 (endcap).  The angle between the
# shower and the nearest charged track must exceed 20 degrees.  The EMC timing
# window suppresses electronic noise and unrelated energy deposits.
event_selection.select_photon do
  energyThreshold_b 0.025   # 25 MeV in the barrel
  energyThreshold_e 0.050   # 50 MeV in the endcaps
  angle_to_track    20.0    # > 20 degrees from the nearest charged track
  tdc_emc_start     0
  tdc_emc_end       14
  nGam              ">=2"
end

# --- Particle identification -------------------------------------------------
# TOF and dE/dx information are combined to give PID probabilities for the
# pion, kaon and proton hypotheses; a track is identified as a kaon when its
# kaon probability is larger than the pion (and proton) probability.
event_selection.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp "==1"
  nkm "==1"
end

# --- Kinematic fit (endpoint) ------------------------------------------------
# 4C kinematic fit under the K+ K- gamma gamma hypothesis, constraining the
# four-momentum of the final state to that of the colliding beams.  With more
# than two photon candidates the combination with the smallest chi2_4C is
# automatically retained.
event_selection.kinematic_fit([:kp, :km, :gamma, :gamma]) do
  nominal
  constrain_four_momentum
  chi2_cut 200   # loose BOSS-side cut; the published chi2_4C < 30 is applied in ROOT
end

### ---------------------------------------------------------------------------
### Inexpressible BOSS-side procedures (captured for the systematics step)
### ---------------------------------------------------------------------------
alg_phipi0
  .note(:kinematic_fit_cut,
        "the published analysis selects events with chi2_4C(K+ K- gamma gamma) < 30; " \
        "only a loose chi2_cut 200 is applied in BOSS and the tight cut is applied " \
        "in the ROOT stage")
  .note(:pi0_mass_window,
        "after the 4C fit the gamma-gamma invariant mass is required to be in the " \
        "pi0 mass region 0.115 < M(gamma gamma) < 0.155 GeV/c^2 (post-kinematic-fit " \
        "selection, therefore applied in ROOT); the same quantity defines the pi0 " \
        "mass sidebands 0.055-0.095 and 0.175-0.215 GeV/c^2 used for sideband " \
        "subtraction")
  .note(:sideband_subtraction,
        "the incoherent background is estimated from pi0 mass sidebands " \
        "(0.055 < M(gamma gamma) < 0.095 GeV/c^2 and 0.175 < M(gamma gamma) < " \
        "0.215 GeV/c^2) and subtracted from the M(K+ K-) distribution")
  .note(:isr_continuum_background,
        "the continuum phi-peaking background e+e- -> gamma_ISR phi is studied using " \
        "data taken at energies far from any charmonium resonance and cannot be " \
        "written as a plain EvtGen exclusive decay card (initial-state radiation)")
  .note(:inclusive_mc_generator,
        "the inclusive J/psi background sample (1.2 x 10^9 events) is generated with " \
        "BesEvtGen for the known decay modes (measured branching fractions) and with " \
        "Lund-Charm for the unknown decays")
  .note(:eta_c_peaking_background,
        "the pi0-peaking background J/psi -> gamma eta_c(1S) -> gamma K+ K- pi0 is " \
        "neglected in the nominal fit (0.5% of the coherent background) and is " \
        "assigned a systematic uncertainty estimated with a dedicated MC sample")
  .note(:mass_resolution,
        "the K+ K- mass resolution sigma_m = (1.00 +/- 0.02) MeV/c^2 is determined " \
        "from a J/psi -> phi eta sample with 0.50 < M(gamma gamma) < 0.60 GeV/c^2 " \
        "and is used in the maximum-likelihood fit to M(K+ K-)")
  .note(:fit_model,
        "maximum-likelihood fit to the sideband-subtracted M(K+ K-) distribution: " \
        "F_H0 = second-order polynomial (non-phi contribution) and F_H1 = coherent " \
        "sum of a relativistic Breit-Wigner phi resonance and the polynomial, " \
        "convoluted with a Gaussian mass resolution; two solutions are found " \
        "(relative phase delta = 95.9 or -152.1 degrees), both at 6.4 sigma")
  .note(:helix_correction,
        "the detection efficiency (45.1 +/- 0.2)% is evaluated from MC with the " \
        "angular distribution of the J/psi -> phi pi0 decay taken into account")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

### ---------------------------------------------------------------------------
### Execution
### ---------------------------------------------------------------------------
alg_phipi0.execute_on([jpsi_data, jpsi_incMC,
                       exMC_signal,
                       exMC_bkg_kkm_pi0, exMC_bkg_phi_pi0pi0, exMC_bkg_phi_gg,
                       exMC_bkg_etac, exMC_bkg_gkk, exMC_bkg_gpi0kk,
                       exMC_bkg_pi0pi0kk, exMC_bkg_phi_eta])
