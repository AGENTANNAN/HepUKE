# =============================================================================
# BESIII measurement of B(psi(3770) -> gamma chi_c1) and search for
# psi(3770) -> gamma chi_c2  (arXiv:1504.07450v2)
# -- BOSS part (dataset preparation + event selection)
#
# The transitions are reconstructed with the chain
#   psi(3770) -> gamma chi_cJ, chi_cJ -> gamma J/psi, J/psi -> l+ l-  (l = e, mu)
# using 2916.94 pb^-1 of data taken at sqrt(s) = 3.773 GeV.
# chi_c1 and chi_c2 share the same final state and the same selection criteria,
# so a single Algorithm per lepton flavour covers both J values (shared final
# state rule).  The e+e- and mu+mu- channels have different lepton identification
# and therefore get separate Algorithm objects.
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets (single energy: the psi(3770) resonance)
### ---------------------------------------------------------------------------
psipp_data  = DatasetManager.real_data.find("712_3773")      # psi(3770), 3773 MeV
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")   # inclusive psi(3770) MC

### ---------------------------------------------------------------------------
### Decay cards
### ---------------------------------------------------------------------------
# Signal: psi(3770) -> gamma chi_c1, chi_c1 -> gamma J/psi, J/psi -> e+ e-
decay_card_chic1_ee = <<~DECAYCARD
  Decay psi(3770)
  1.0000 gamma chi_c1 PHSP;
  Enddecay

  Decay chi_c1
  1.000 gamma J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Signal: psi(3770) -> gamma chi_c2, chi_c2 -> gamma J/psi, J/psi -> e+ e-
decay_card_chic2_ee = <<~DECAYCARD
  Decay psi(3770)
  1.0000 gamma chi_c2 PHSP;
  Enddecay

  Decay chi_c2
  1.000 gamma J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Signal: psi(3770) -> gamma chi_c1, chi_c1 -> gamma J/psi, J/psi -> mu+ mu-
decay_card_chic1_mumu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 gamma chi_c1 PHSP;
  Enddecay

  Decay chi_c1
  1.000 gamma J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 mu+ mu- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Signal: psi(3770) -> gamma chi_c2, chi_c2 -> gamma J/psi, J/psi -> mu+ mu-
decay_card_chic2_mumu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 gamma chi_c2 PHSP;
  Enddecay

  Decay chi_c2
  1.000 gamma J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 mu+ mu- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# The expected angular distributions of psi(3770) -> gamma chi_cJ are taken
# into account in the signal MC generation (PHSP is used here as the fallback
# generator since the dedicated model is not expressible in the DSL).

### ---------------------------------------------------------------------------
### Exclusive MC
### ---------------------------------------------------------------------------
exMC_chic1_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psipp_gamma_chic1_ee"
  config.related_dataset = psipp_data
  config.events          = 500_000
  config.decay_card      = decay_card_chic1_ee
  config.cross_section   = :default
end

exMC_chic2_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psipp_gamma_chic2_ee"
  config.related_dataset = psipp_data
  config.events          = 500_000
  config.decay_card      = decay_card_chic2_ee
  config.cross_section   = :default
end

exMC_chic1_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psipp_gamma_chic1_mumu"
  config.related_dataset = psipp_data
  config.events          = 500_000
  config.decay_card      = decay_card_chic1_mumu
  config.cross_section   = :default
end

exMC_chic2_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psipp_gamma_chic2_mumu"
  config.related_dataset = psipp_data
  config.events          = 500_000
  config.decay_card      = decay_card_chic2_mumu
  config.cross_section   = :default
end

### ===========================================================================
### Algorithm I -- gamma gamma e+ e- channel
### ===========================================================================
alg_ee = Algorithm.new("PsippToGammaChiC_EE")
alg_ee.set_header(["PsippToGammaChiC_EEAlg/PsippToGammaChiC_EE.h"])
      .set_constant({ "ECMS" => [:double, 3.773] })   # psi(3770) c.m. energy

sel_ee = Selection.new

# --- Charged track selection -------------------------------------------------
# Charged tracks are reconstructed from hits in the MDC; |cos(theta)| < 0.93,
# distance of closest approach to the average e+e- interaction point less than
# 1.0 cm in the plane perpendicular to the beam and less than 15.0 cm along the
# beam direction.  Exactly two oppositely charged tracks are required.
sel_ee.select_track do
  cos_theta 0.93
  Vz        15.0
  Vr        1.0
  nChrp     "==1"
  nChrn     "==1"
  nNet      "==0"
end

# --- Photon selection --------------------------------------------------------
# Two good photon candidates; the deposited energy of a neutral EMC cluster must
# exceed 50 MeV, the EMC timing suppresses electronic noise and unrelated energy
# deposits, and the angle to the nearest charged track must exceed 10 degrees.
sel_ee.select_photon do
  energyThreshold_b 0.050   # > 50 MeV
  energyThreshold_e 0.050   # > 50 MeV
  angle_to_track    10.0    # > 10 degrees from the nearest charged track
  tdc_emc_start     0
  tdc_emc_end       14
  nGam              ">=2"
end

# --- Lepton identification ---------------------------------------------------
# Electrons and muons are separated by the ratio E/p (E = energy deposited in
# the EMC, p = momentum measured in the MDC): E/p > 0.7 identifies an electron
# or positron; otherwise an EMC energy between 0.05 and 0.35 GeV identifies a
# muon.  The J/psi candidates are built from lepton pairs with momenta between
# 1.2 and 1.9 GeV/c.
sel_ee.pid(method: :probability) do
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.2,
                                 treat_as_electron_if_energy_above: 0.6
  nlp "==1"
  nlm "==1"
end

# --- Kinematic fit (endpoint) ------------------------------------------------
# 5C kinematic fit under the gamma gamma e+ e- hypothesis: the total energy and
# the three-momentum components are constrained to the expected centre-of-mass
# values (the small beam crossing angle is taken into account), plus the
# invariant mass of the e+ e- pair is constrained to the J/psi mass.
# The pi0 and eta vetoes on M(gamma gamma) are applied on the photon pair.
sel_ee.kinematic_fit([:gamma, :gamma, :lp, :lm]) do
  nominal
  invariant_mass_of(:gamma, :gamma).out_of(0.124, 0.146).out_of(0.537, 0.558)
  invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
  constrain_four_momentum
  chi2_cut 200   # loose BOSS-side cut; the published chi2_5C < 25 is applied in ROOT
end

alg_ee
  .note(:lepton_id_method,
        "electron/muon separation uses the ratio E/p: E/p > 0.7 identifies an " \
        "electron or positron, otherwise an EMC energy deposit between 0.05 and " \
        "0.35 GeV identifies a muon. identify_high_momentum_leptons with " \
        "treat_as_electron_if_energy_above is used as the closest DSL equivalent")
  .note(:lepton_momentum_window,
        "the J/psi candidates are reconstructed from lepton pairs with momenta in " \
        "the range 1.2-1.9 GeV/c")
  .note(:electron_angular_cut,
        "for the gamma gamma e+ e- mode the positron and electron polar angles are " \
        "required to satisfy cos(theta_e+) < 0.5 and cos(theta_e-) > -0.5 to reduce " \
        "radiative Bhabha background; charge-dependent cos(theta) cuts are not " \
        "expressible with the per-candidate DSL properties")
  .note(:pi0_eta_veto,
        "events with the two-photon invariant mass inside the pi0 mass window " \
        "(0.124, 0.146) GeV/c^2 or the eta mass window (0.537, 0.558) GeV/c^2 are " \
        "excluded to remove J/psi pi0 and J/psi eta background with pi0/eta -> gamma gamma")
  .note(:higher_energy_photon,
        "the chi_c1 and chi_c2 are reconstructed from the invariant mass of " \
        "gamma^H J/psi, where gamma^H denotes the higher-energy photon of the " \
        "gamma gamma l+ l- final state; this candidate choice is applied in ROOT")
  .note(:background_generators,
        "background MC samples generated for this analysis: inclusive psi(3770) " \
        "decays (EvtGen + LundCharm for unknown modes), e+e- -> (gamma) J/psi, " \
        "e+e- -> (gamma) psi(3686), e+e- -> q qbar (q = u, d, s), e+e- -> tau+ tau- " \
        "(kkmc), and e+e- -> (gamma) e+e- / (gamma) mu+ mu- (Babayaga)")
  .note(:isr_psi3686_peaking_background,
        "e+e- -> (gamma_ISR) psi(3686) with psi(3686) -> gamma chi_cJ gives a peaking " \
        "background with the same topology as the signal; its yield is fixed in the " \
        "fit to the predicted number of events (5.3 for chi_c0, 225.4 for chi_c1, " \
        "158.4 for chi_c2) computed from the observed ISR cross section, the " \
        "luminosity, the branching fractions and the MC misidentification rates")
  .note(:helix_correction,
        "the uncertainty associated with the kinematic fit is determined by " \
        "comparing the chi2 distributions and the chi2 < 25 efficiencies for data " \
        "and MC using psi(3686) -> gamma gamma l+ l- events at sqrt(s) = 3.686 GeV")
  .with_decay_card(decay_card_chic1_ee)
  .apply(sel_ee)

### ===========================================================================
### Algorithm II -- gamma gamma mu+ mu- channel
### ===========================================================================
alg_mumu = Algorithm.new("PsippToGammaChiC_MuMu")
alg_mumu.set_header(["PsippToGammaChiC_MuMuAlg/PsippToGammaChiC_MuMu.h"])
        .set_constant({ "ECMS" => [:double, 3.773] })   # psi(3770) c.m. energy

# The selection is identical to the e+ e- channel except for the lepton flavour.
sel_mumu = Selection.new

sel_mumu.select_track do
  cos_theta 0.93
  Vz        15.0
  Vr        1.0
  nChrp     "==1"
  nChrn     "==1"
  nNet      "==0"
end

sel_mumu.select_photon do
  energyThreshold_b 0.050   # > 50 MeV
  energyThreshold_e 0.050   # > 50 MeV
  angle_to_track    10.0    # > 10 degrees from the nearest charged track
  tdc_emc_start     0
  tdc_emc_end       14
  nGam              ">=2"
end

sel_mumu.pid(method: :probability) do
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.2,
                                 treat_as_electron_if_energy_above: 0.6
  nlp "==1"
  nlm "==1"
end

sel_mumu.kinematic_fit([:gamma, :gamma, :lp, :lm]) do
  nominal
  invariant_mass_of(:gamma, :gamma).out_of(0.124, 0.146).out_of(0.537, 0.558)
  invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
  constrain_four_momentum
  chi2_cut 200   # loose BOSS-side cut; the published chi2_5C < 25 is applied in ROOT
end

alg_mumu
  .note(:lepton_id_method,
        "muon identification: E/p < 0.7 together with an EMC energy deposit in the " \
        "range 0.05-0.35 GeV; the MUC provides additional muon information for " \
        "momenta above 0.5 GeV/c")
  .note(:lepton_momentum_window,
        "the J/psi candidates are reconstructed from lepton pairs with momenta in " \
        "the range 1.2-1.9 GeV/c")
  .note(:pi0_eta_veto,
        "events with the two-photon invariant mass inside the pi0 mass window " \
        "(0.124, 0.146) GeV/c^2 or the eta mass window (0.537, 0.558) GeV/c^2 are " \
        "excluded to remove J/psi pi0 and J/psi eta background with pi0/eta -> gamma gamma")
  .note(:higher_energy_photon,
        "the chi_c1 and chi_c2 are reconstructed from the invariant mass of " \
        "gamma^H J/psi, where gamma^H denotes the higher-energy photon of the " \
        "gamma gamma l+ l- final state; this candidate choice is applied in ROOT")
  .note(:cross_contamination,
        "cross contamination between the e+e- and mu+mu- modes of the signal events " \
        "was studied with MC and found to be negligible")
  .note(:isr_psi3686_peaking_background,
        "e+e- -> (gamma_ISR) psi(3686) with psi(3686) -> gamma chi_cJ is an " \
        "irreducible peaking background with a yield fixed in the fit from the " \
        "predicted number of background events")
  .with_decay_card(decay_card_chic1_mumu)
  .apply(sel_mumu)

### ---------------------------------------------------------------------------
### Execution
### ---------------------------------------------------------------------------
alg_ee.execute_on([psipp_data, psipp_incMC,
                   exMC_chic1_ee, exMC_chic2_ee])

alg_mumu.execute_on([psipp_data, psipp_incMC,
                     exMC_chic1_mumu, exMC_chic2_mumu])
