# -*- coding: utf-8 -*-
# =============================================================================
# 1707.07042v2
# Measurements of the branching fractions of chi_c0,2 -> eta' eta' and eta eta'
# (BESIII, 448.1 million psi(3686) events)
#
# BOSS-side spec: dataset preparation + event selection up to the final
# 4C (+ eta mass constraint) kinematic fit. The simultaneous fits to the
# M(eta' eta') and M(eta eta') spectra, the background/chi_c1 component
# modelling, the sideband normalisation and the branching-fraction extraction
# are ROOT-level and are not part of this spec.
#
# Analysis strategy: psi(3686) -> gamma chi_cJ with the E1 transition photon
# reconstructed, followed by chi_cJ -> eta' eta' / eta eta'. Five decay modes
# are studied; each has a different final-state multiplicity and a different
# kinematic-fit participant list, so each is an independent Algorithm +
# Selection (Rule T1). The transition photon is part of the fitted final state,
# exactly as in the vetted chi_cJ analyses 1011.6556v2 / 1610.02479v3.
# All charged tracks are assumed to be pions (no PID block).
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets — 448.1 x 10^6 psi(3686) events, plus 48 pb^-1 of 3.65 GeV
### continuum data used as a background control sample.
### ---------------------------------------------------------------------------
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")
cont_data  = DatasetManager.real_data.find("709_3650")   # 48 pb^-1 at sqrt(s) = 3.65 GeV
cont_incMC = DatasetManager.inclusive_mc.find("709_3650")

### ---------------------------------------------------------------------------
### Decay cards (EvtGen). The psi(3686) -> gamma chi_cJ transition is E1 and is
### modelled with HELAMP. chi_c0 and chi_c2 are generated with their PDG
### production branching fractions in a single card, since both are present in
### the data; the subsequent eta'/eta decays use the PHSP uniform phase-space
### assumption of the paper.
### ---------------------------------------------------------------------------

# Mode A: chi_c0,2 -> eta' eta', both eta' -> gamma pi+ pi-
# Final state fitted: gamma(transition) + gamma pi+ pi- + gamma pi+ pi-
decay_card_modeA = <<~DECAYCARD
    Decay psi(2S)
    0.0970 gamma chi_c0                       HELAMP 1.0 0.0 1.0 0.0;
    0.0910 gamma chi_c2                       HELAMP 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay chi_c0
    1.0000 eta' eta'                          PHSP;
    Enddecay

    Decay chi_c2
    1.0000 eta' eta'                          PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi-                      PHSP;
    Enddecay

    End
DECAYCARD

# Mode B: chi_c0,2 -> eta' eta', both eta' -> eta pi+ pi- (eta -> gamma gamma)
# Final state fitted: gamma(transition) + (gamma gamma) pi+ pi- + (gamma gamma) pi+ pi-
decay_card_modeB = <<~DECAYCARD
    Decay psi(2S)
    0.0970 gamma chi_c0                       HELAMP 1.0 0.0 1.0 0.0;
    0.0910 gamma chi_c2                       HELAMP 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay chi_c0
    1.0000 eta' eta'                          PHSP;
    Enddecay

    Decay chi_c2
    1.0000 eta' eta'                          PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi-                        PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma                        PHSP;
    Enddecay

    End
DECAYCARD

# Mode C: chi_c0,2 -> eta' eta', one eta' -> gamma pi+ pi-, the other
# eta' -> eta pi+ pi- (eta -> gamma gamma). The two eta' are distinct here, so
# Alias is used to keep their decay nodes separate (a single Decay eta' block
# cannot carry two different final states).
decay_card_modeC = <<~DECAYCARD
    Alias etap_gamma eta'
    Alias etap_eta eta'

    Decay psi(2S)
    0.0970 gamma chi_c0                       HELAMP 1.0 0.0 1.0 0.0;
    0.0910 gamma chi_c2                       HELAMP 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay chi_c0
    1.0000 etap_gamma etap_eta                PHSP;
    Enddecay

    Decay chi_c2
    1.0000 etap_gamma etap_eta                PHSP;
    Enddecay

    Decay etap_gamma
    1.0000 gamma pi+ pi-                      PHSP;
    Enddecay

    Decay etap_eta
    1.0000 eta pi+ pi-                        PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma                        PHSP;
    Enddecay

    End
DECAYCARD

# Mode I: chi_c0,2 -> eta eta', eta -> gamma gamma, eta' -> gamma pi+ pi-
# Final state fitted: gamma(transition) + (gamma gamma) + gamma pi+ pi-
decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    0.0970 gamma chi_c0                       HELAMP 1.0 0.0 1.0 0.0;
    0.0910 gamma chi_c2                       HELAMP 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay chi_c0
    1.0000 eta eta'                           PHSP;
    Enddecay

    Decay chi_c2
    1.0000 eta eta'                           PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma                        PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi-                      PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: chi_c0,2 -> eta eta', eta -> gamma gamma in both the direct eta and
# in eta' -> eta pi+ pi-
# Final state fitted: gamma(transition) + (gamma gamma) + (gamma gamma) pi+ pi-
decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    0.0970 gamma chi_c0                       HELAMP 1.0 0.0 1.0 0.0;
    0.0910 gamma chi_c2                       HELAMP 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay chi_c0
    1.0000 eta eta'                           PHSP;
    Enddecay

    Decay chi_c2
    1.0000 eta eta'                           PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma                        PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi-                        PHSP;
    Enddecay

    End
DECAYCARD

### ---------------------------------------------------------------------------
### Exclusive signal MC — one sample per decay mode
### ---------------------------------------------------------------------------
exMC_modeA = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chic_etapetap_2gammapipi"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_modeA
  config.cross_section   = :default
end

exMC_modeB = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chic_etapetap_2etapipi"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_modeB
  config.cross_section   = :default
end

exMC_modeC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chic_etapetap_gammapipi_etapipi"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_modeC
  config.cross_section   = :default
end

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chic_etaetap_gammapipi"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chic_etaetap_etapipi"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### ---------------------------------------------------------------------------
### Event selection (BOSS) — chi_c0,2 -> eta' eta', Mode A
###   (both eta' -> gamma pi+ pi-): 4 charged tracks + 3 photons
### ---------------------------------------------------------------------------
alg_name_A = "PsipChicEtapEtapModeA"
alg_A = Algorithm.new(alg_name_A)
alg_A.set_header(["#{alg_name_A}Alg/#{alg_name_A}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

selection_A = Selection.new
  .select_track do
    cos_theta 0.93     # |cos(theta)| < 0.93
    Vz        10.0     # |Vz| < 10 cm
    Vr        1.0      # R_xy < 1 cm
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14      # EMC hit time 0 <= t <= 700 ns
    energyThreshold_b 0.025   # E > 25 MeV in the barrel (|cos(theta)| < 0.8)
    energyThreshold_e 0.050   # E > 50 MeV in the endcaps (0.86 < |cos(theta)| < 0.92)
    nGam              ">=3"   # 1 transition photon + 1 photon per eta' decay
  end
  # All charged tracks are assumed to be pions (no PID requirement).
  # Nominal 4C fit of the full final state gamma + (gamma pi+ pi-) + (gamma pi+ pi-)
  # to the initial beam four-momentum; the combination of photons with the
  # smallest chi2 is retained by the fit itself. The published decay-mode
  # dependent chi2 cut (25-90) is applied at the ROOT stage.
  .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pip, :pim, :pim]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_A
  .note(:track_assumption, "all charged tracks are assumed to be pions.")
  .note(:chi2_cut_published, "the published chi2 requirement is a decay-mode dependent value
    in the range 25-90, obtained by optimising the figure-of-merit
    N_S^MC / sqrt(N_S^data + N_B^data); BOSS applies the loose chi2_cut 200 and the
    optimal cut is applied in ROOT.")
  .note(:best_photon_combination, "if additional photons are found in an event, the
    combination of photons with the least chi2 is retained for further analysis; the
    kinematic fit iterates over the photon combinations automatically.")
  .note(:etap_candidate_selection, "the two eta' candidates are selected by minimising
    Delta = sqrt((M_1 - M_eta')^2 + (M_2 - M_eta')^2) (equivalently (M_i - M_eta')^2 +
    (M_j - M_eta')^2), with M_eta' the nominal eta' mass. This uses post-fit quantities and
    is applied at the ROOT level.")
  .note(:etap_signal_region, "the double-eta' signal region is M(gamma pi+ pi-) in
    (0.943, 0.973) GeV/c^2 for mode A.")
  .note(:background_veto_pi0, "psi(3686) -> pi0 + X background suppressed by requiring the
    invariant mass of any two photons to be out of the pi0 mass region,
    |M(gamma gamma) - M_pi0| > 15 MeV/c^2.")
  .note(:background_veto_jpsi, "psi(3686) -> pi+ pi- J/psi background suppressed by requiring
    the recoil mass of any pi+ pi- pair, M_rec(pi+ pi-), to satisfy
    |M_rec(pi+ pi-) - M_J/psi| > 5 MeV/c^2.")
  .note(:sideband_background, "the eight equal-area boxes surrounding the signal region in the
    M_1 vs M_2 plane are the sidebands: the four corner boxes (type A) estimate the
    background without eta' in the subsequent decay, the other four (type B) the background
    with one eta'. The expected background is M_side^B/2 - M_side^A/4, assuming the
    background is uniform around the signal region. The chi_c0,2 signals are extracted from
    a simultaneous fit of the three M(eta' eta') spectra (modes A, B and C) with common
    B(chi_c0,2 -> eta' eta') parameters, an MC-derived signal shape convolved with a
    Gaussian, an MC-derived chi_c1 peaking background for modes A and C
    (chi_c1 -> gamma J/psi, J/psi -> gamma 2(pi+ pi-) and chi_c1 -> f_0(980) eta') and a
    first-order Chebyshev non-peaking background.")
  .note(:angular_distribution, "the subsequent chi_c0 -> eta' eta' and chi_c2 -> eta' eta'
    decays are generated flat in phase space, while the eta'/eta angular distribution in
    chi_c2 decays follows the pi+- measurement of Ref. [16]; the PHSP signal MC is used for
    the efficiency and the model dependence is a systematic uncertainty.")
  .with_decay_card(decay_card_modeA)
  .apply(selection_A)

### ---------------------------------------------------------------------------
### Event selection (BOSS) — chi_c0,2 -> eta' eta', Mode B
###   (both eta' -> eta pi+ pi-, eta -> gamma gamma): 4 tracks + 5 photons
### ---------------------------------------------------------------------------
alg_name_B = "PsipChicEtapEtapModeB"
alg_B = Algorithm.new(alg_name_B)
alg_B.set_header(["#{alg_name_B}Alg/#{alg_name_B}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

selection_B = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=5"   # 1 transition photon + 2 photons per eta (two eta' decays)
  end
  # Both eta candidates are reconstructed from independent photon pairs, each
  # gamma gamma pair mass-constrained to the nominal eta mass.
  .kalman_kinematic_fit([:gamma, :gamma, :gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=2"
  end
  # Nominal 4C fit of the full final state
  # gamma + (gamma gamma) pi+ pi- + (gamma gamma) pi+ pi- to the initial beam
  # four-momentum, with the two eta mass constraints from the Kalman fit.
  .kinematic_fit([:gamma, :eta, :eta, :pip, :pip, :pim, :pim]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_B
  .note(:track_assumption, "all charged tracks are assumed to be pions.")
  .note(:chi2_cut_published, "the published chi2 requirement is a decay-mode dependent value
    in the range 25-90, obtained by optimising the figure-of-merit
    N_S^MC / sqrt(N_S^data + N_B^data); the loose BOSS chi2_cut 200 is used here and the
    optimal cut is applied in ROOT.")
  .note(:eta_mass_window, "an eta candidate is a photon pair with
    |M(gamma gamma) - M_eta| < 20 MeV/c^2; this is enforced by the Kalman mass constraint and
    the corresponding window is refined at the ROOT level.")
  .note(:best_photon_combination, "if additional photons are found in an event, the
    combination of photons with the least chi2 is retained for further analysis.")
  .note(:etap_candidate_selection, "the two eta' candidates are selected by minimising
    (M_i - M_eta')^2 + (M_j - M_eta')^2 with M_eta' the nominal eta' mass; the eta' -> eta pi+ pi-
    mass M_2 is used for this mode. Applied on post-fit quantities at the ROOT level.")
  .note(:etap_signal_region, "the double-eta' signal region is M(eta pi+ pi-) in
    (0.928, 0.988) GeV/c^2 for mode B.")
  .note(:background_veto_pi0, "psi(3686) -> pi0 + X background suppressed by requiring the
    invariant mass of any two photons to be out of the pi0 mass region,
    |M(gamma gamma) - M_pi0| > 15 MeV/c^2.")
  .note(:background_veto_jpsi, "psi(3686) -> pi+ pi- J/psi background suppressed by requiring
    |M_rec(pi+ pi-) - M_J/psi| > 5 MeV/c^2.")
  .note(:sideband_background, "corner (type A) and edge (type B) sideband boxes are used to
    estimate the non-peaking background as M_side^B/2 - M_side^A/4; the yield is obtained from a
    simultaneous fit of the three M(eta' eta') spectra with common branching fractions.")
  .note(:angular_distribution, "the chi_c0,2 -> eta' eta' decays are generated flat in phase
    space; the angular distribution of the eta'/eta in chi_c2 decays follows Ref. [16] and the
    model dependence is a systematic uncertainty.")
  .with_decay_card(decay_card_modeB)
  .apply(selection_B)

### ---------------------------------------------------------------------------
### Event selection (BOSS) — chi_c0,2 -> eta' eta', Mode C
###   (one eta' -> gamma pi+ pi-, the other -> eta pi+ pi-): 4 tracks + 4 photons
### ---------------------------------------------------------------------------
alg_name_C = "PsipChicEtapEtapModeC"
alg_C = Algorithm.new(alg_name_C)
alg_C.set_header(["#{alg_name_C}Alg/#{alg_name_C}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

selection_C = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"   # 1 transition + 1 (eta' -> gamma pi+ pi-) + 2 (eta -> gamma gamma)
  end
  # The eta from eta' -> eta pi+ pi- is reconstructed from a photon pair
  # mass-constrained to the nominal eta mass.
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  end
  # Nominal 4C fit of the full final state
  # gamma + gamma pi+ pi- + (gamma gamma) pi+ pi- to the initial beam four-momentum,
  # with the eta mass constraint from the Kalman fit.
  .kinematic_fit([:gamma, :gamma, :eta, :pip, :pip, :pim, :pim]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_C
  .note(:track_assumption, "all charged tracks are assumed to be pions.")
  .note(:chi2_cut_published, "the published chi2 requirement is a decay-mode dependent value
    in the range 25-90, obtained by optimising the figure-of-merit
    N_S^MC / sqrt(N_S^data + N_B^data); the loose BOSS chi2_cut 200 is used here.")
  .note(:eta_mass_window, "an eta candidate is a photon pair with
    |M(gamma gamma) - M_eta| < 20 MeV/c^2, enforced by the Kalman mass constraint.")
  .note(:best_photon_combination, "the photon combination with the least chi2 is retained when
    extra photons are present.")
  .note(:etap_candidate_selection, "the two eta' candidates are selected by minimising
    (M_i - M_eta')^2 + (M_j - M_eta')^2, with the index 1 denoting the gamma pi+ pi- (mode A type)
    and index 2 the eta pi+ pi- (mode B type) eta' candidate.")
  .note(:etap_signal_region, "the double-eta' signal region is M_1 in (0.933, 0.983) GeV/c^2 and
    M_2 in (0.943, 0.973) GeV/c^2 for mode C.")
  .note(:identical_particle_factor, "for mode C a factor of two accounts for the identical
    particles (the two eta' decay differently); this enters the branching-fraction extraction in
    ROOT.")
  .note(:background_veto_pi0, "psi(3686) -> pi0 + X background suppressed by requiring the
    invariant mass of any two photons to be out of the pi0 mass region,
    |M(gamma gamma) - M_pi0| > 15 MeV/c^2.")
  .note(:background_veto_jpsi, "psi(3686) -> pi+ pi- J/psi background suppressed by requiring
    |M_rec(pi+ pi-) - M_J/psi| > 5 MeV/c^2.")
  .note(:chi_c1_background, "a chi_c1 peaking background is present in mode C from
    chi_c1 -> f_0(980) eta' and is described by the corresponding MC shape with a floated
    normalisation in the ROOT fit; no chi_c0,2 peaks are observed in the sidebands.")
  .note(:sideband_background, "the non-peaking background is estimated from the eight sideband
    boxes around the signal region and modelled by a first-order Chebyshev polynomial.")
  .note(:angular_distribution, "the chi_c0,2 -> eta' eta' decays are generated flat in phase
    space; the chi_c2 angular distribution follows Ref. [16].")
  .with_decay_card(decay_card_modeC)
  .apply(selection_C)

### ---------------------------------------------------------------------------
### Event selection (BOSS) — chi_c0,2 -> eta eta', Mode I
###   (eta' -> gamma pi+ pi-, eta -> gamma gamma): 2 tracks + 4 photons
### ---------------------------------------------------------------------------
alg_name_I = "PsipChicEtaEtapModeI"
alg_I = Algorithm.new(alg_name_I)
alg_I.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})

selection_I = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"   # 1 transition + 2 (eta -> gamma gamma) + 1 (eta' -> gamma pi+ pi-)
  end
  # The eta is reconstructed from a photon pair mass-constrained to the nominal eta mass.
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  end
  # Nominal 4C fit of the full final state
  # gamma + (gamma gamma) + gamma pi+ pi- to the initial beam four-momentum.
  .kinematic_fit([:gamma, :gamma, :eta, :pip, :pim]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_I
  .note(:track_assumption, "all charged tracks are assumed to be pions.")
  .note(:chi2_cut_published, "the published chi2 requirement is a decay-mode dependent value
    in the range 25-90, obtained by optimising the figure-of-merit
    N_S^MC / sqrt(N_S^data + N_B^data); the loose BOSS chi2_cut 200 is used here and the
    optimal value is applied in ROOT.")
  .note(:eta_mass_window, "an eta candidate is a photon pair with
    |M(gamma gamma) - M_eta| < 20 MeV/c^2, enforced by the Kalman mass constraint.")
  .note(:best_photon_combination, "the photon combination with the least chi2 is retained when
    extra photons are present.")
  .note(:etap_candidate_selection, "the eta' candidate is selected by minimising
    |M(gamma pi+ pi-) - M_eta'| with M_eta' the nominal eta' mass; applied on post-fit
    quantities in ROOT.")
  .note(:etap_signal_region, "the eta' signal region is M_1 in (0.948, 0.968) GeV/c^2 for
    mode I; two sideband regions of the same width are chosen around it.")
  .note(:background_veto_pi0, "psi(3686) -> pi0 + X background suppressed by requiring the
    invariant mass of any two photons to be out of the pi0 mass region,
    |M(gamma gamma) - M_pi0| > 15 MeV/c^2.")
  .note(:background_veto_jpsi, "psi(3686) -> pi+ pi- J/psi background suppressed by requiring
    |M_rec(pi+ pi-) - M_J/psi| > 5 MeV/c^2.")
  .note(:background_veto_gammajpsi_etap, "in the chi_c0,2 -> eta eta' channel the background
    chi_c0,2 -> gamma J/psi, J/psi -> gamma eta' is suppressed by requiring the invariant mass of
    any gamma eta' combination (the gamma being one of the eta daughters) to be out of
    (3.05, 3.16) GeV/c^2 for mode I.")
  .note(:chi_c1_background, "a chi_c1 peaking background from psi(3686) -> gamma chi_c1,
    chi_c1 -> gamma J/psi with J/psi -> gamma gamma pi+ pi- (or gamma eta', eta' -> gamma pi+ pi-)
    is present in mode I and is modelled by its MC shape with a floated normalisation.")
  .note(:simultaneous_fit, "the two M(eta eta') spectra (modes I and II) are fitted
    simultaneously with common B(chi_c0,2 -> eta eta') parameters; the signal shape is the
    MC histogram convolved with a Gaussian, the non-peaking background a first-order Chebyshev
    polynomial.")
  .note(:angular_distribution, "the chi_c0,2 -> eta eta' decays are generated flat in phase
    space; the chi_c2 angular distribution follows the pi+- measurement of Ref. [16].")
  .with_decay_card(decay_card_modeI)
  .apply(selection_I)

### ---------------------------------------------------------------------------
### Event selection (BOSS) — chi_c0,2 -> eta eta', Mode II
###   (eta' -> eta pi+ pi- with eta -> gamma gamma): 2 tracks + 5 photons
### ---------------------------------------------------------------------------
alg_name_II = "PsipChicEtaEtapModeII"
alg_II = Algorithm.new(alg_name_II)
alg_II.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})

selection_II = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=5"   # 1 transition + 2 (direct eta) + 2 (eta from eta')
  end
  # Two eta candidates (the direct eta and the eta from eta' -> eta pi+ pi-), each
  # a photon pair mass-constrained to the nominal eta mass.
  .kalman_kinematic_fit([:gamma, :gamma, :gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=2"
  end
  # Nominal 4C fit of the full final state
  # gamma + (gamma gamma) + (gamma gamma) pi+ pi- to the initial beam four-momentum.
  .kinematic_fit([:gamma, :eta, :eta, :pip, :pim]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_II
  .note(:track_assumption, "all charged tracks are assumed to be pions.")
  .note(:chi2_cut_published, "the published chi2 requirement is a decay-mode dependent value
    in the range 25-90; the loose BOSS chi2_cut 200 is used here and the optimal value is
    applied in ROOT.")
  .note(:eta_mass_window, "an eta candidate is a photon pair with
    |M(gamma gamma) - M_eta| < 20 MeV/c^2, enforced by the Kalman mass constraints.")
  .note(:best_photon_combination, "the photon combination with the least chi2 is retained when
    extra photons are present.")
  .note(:etap_candidate_selection, "the eta' candidate is selected by minimising
    |M(eta pi+ pi-) - M_eta'| with M_eta' the nominal eta' mass.")
  .note(:etap_signal_region, "the eta' signal region is M_2 in (0.943, 0.973) GeV/c^2 for
    mode II; two sideband regions of the same width are chosen around it.")
  .note(:background_veto_pi0, "psi(3686) -> pi0 + X background suppressed by requiring the
    invariant mass of any two photons to be out of the pi0 mass region,
    |M(gamma gamma) - M_pi0| > 15 MeV/c^2.")
  .note(:background_veto_jpsi, "psi(3686) -> pi+ pi- J/psi background suppressed by requiring
    |M_rec(pi+ pi-) - M_J/psi| > 5 MeV/c^2.")
  .note(:background_veto_gammajpsi_etap, "the chi_c0,2 -> gamma J/psi, J/psi -> gamma eta'
    background is suppressed by requiring the invariant mass of any gamma eta' combination to be
    out of (3.049, 3.199) GeV/c^2 for mode II.")
  .note(:simultaneous_fit, "the two M(eta eta') spectra are fitted simultaneously with common
    B(chi_c0,2 -> eta eta') parameters and a first-order Chebyshev non-peaking background.")
  .note(:angular_distribution, "the chi_c0,2 -> eta eta' decays are generated flat in phase
    space; the chi_c2 angular distribution follows Ref. [16].")
  .with_decay_card(decay_card_modeII)
  .apply(selection_II)

### ---------------------------------------------------------------------------
### BOSS-side procedures with no dedicated DSL construct (ROOT-level)
### ---------------------------------------------------------------------------
[alg_A, alg_B, alg_C, alg_I, alg_II].each do |alg|
  alg.note(:detector_and_generator,
          "the psi(3686) production is simulated with KKMC, including beam energy spread and
    initial-state radiation; known decays use EVTGEN with PDG branching fractions and the
    remaining unknown decays use LUNDCHARM. The psi(3686) -> gamma chi_cJ transition is a pure
    E1 process and the chi_c0,2 -> eta' eta' / eta eta' decays are generated flat in phase space.")
  alg.note(:continuum_background,
          "48 pb^-1 of data collected at sqrt(s) = 3.65 GeV, about one fifteenth of the
    psi(3686) integrated luminosity, is used to check the continuum background; almost no events
    survive the selection, so the continuum contribution is neglected.")
  alg.note(:inclusive_mc_background_study,
          "an inclusive MC sample of 3.64 x 10^8 psi(3686) events is used to study the potential
    backgrounds; the cross contamination between the different decay modes is found to be
    negligible.")
  alg.note(:signal_extraction,
          "the number of signal events is obtained from an unbinned maximum-likelihood
    simultaneous fit to the M(eta' eta') and M(eta eta') spectra, with the chi_c0,2 signal shape
    taken from signal MC and convolved with a Gaussian whose parameters are fixed from control
    samples (psi(3686) -> gamma chi_c0,2 with chi_c0,2 -> 2(pi+ pi-) and 2(pi+ pi- pi0)),
    the chi_c1 peaking background from dedicated MC and the non-peaking background from a
    first-order Chebyshev polynomial; the statistical significances are obtained from the
    change in likelihood with and without the corresponding signal component.")
  alg.note(:systematic_uncertainties,
          "the systematic uncertainties comprise the total number of psi(3686) events (0.7%),
    MDC tracking (<1% per charged track), photon detection efficiency (0.6% per photon,
    angularly weighted), eta reconstruction (1.0% per eta), the eta' mass window, the kinematic
    fit (track-parameter correction determined with psi(3686) -> pi+ pi- K+ K-), the fit
    procedure (resolution compensation, fit range, background shape), the branching fractions of
    the intermediate states, and the background-subtraction vetoes.")
end

### ---------------------------------------------------------------------------
### Execution
### ---------------------------------------------------------------------------
datasets_all = [psip_data, psip_incMC, cont_data, cont_incMC]

root_files_A  = alg_A.execute_on(datasets_all + [exMC_modeA])
root_files_B  = alg_B.execute_on(datasets_all + [exMC_modeB])
root_files_C  = alg_C.execute_on(datasets_all + [exMC_modeC])
root_files_I  = alg_I.execute_on(datasets_all + [exMC_modeI])
root_files_II = alg_II.execute_on(datasets_all + [exMC_modeII])
