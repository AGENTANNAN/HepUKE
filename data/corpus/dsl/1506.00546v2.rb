# ============================================================================
# arXiv:1506.00546v2 (BESIII)
# Mass-independent amplitude analysis of the pi0 pi0 system in radiative
# J/psi decays: J/psi -> gamma pi0 pi0, pi0 -> gamma gamma.
#
# Data: (1.311 +/- 0.011) x 10^9 J/psi events collected at sqrt(s) = 3.097 GeV.
# The amplitude analysis (unbinned extended maximum-likelihood fits in bins of
# M(pi0 pi0)) is the ROOT stage; here we cover dataset preparation, the signal
# and background decay cards, the exclusive MC and the selection chain up to
# and including the kinematic fit.
#
# BOSS part only: decay cards, exclusive MC, event selection up to the 6C
# kinematic fit.
# ============================================================================

### Dataset description ###
# J/psi data at 3.097 GeV and the corresponding inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Continuum data at 3.080 GeV, used to estimate the e+e- -> gamma pi0 pi0
# continuum background
cont_data_3080 = DatasetManager.real_data.find("708_3080")

# ----------------------------------------------------------------------------
# Signal decay card: J/psi -> gamma pi0 pi0 with pi0 -> gamma gamma.
# The five-photon final state (one radiated photon plus two pi0 daughters)
# is generated according to phase space.
# ----------------------------------------------------------------------------
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ----------------------------------------------------------------------------
# Background decay cards for the two dominant peaking backgrounds
# J/psi -> gamma eta and J/psi -> gamma eta', which are simulated with
# exclusive MC according to the PDG branching fractions.
# ----------------------------------------------------------------------------
# J/psi -> gamma eta, eta -> pi0 pi0 pi0, pi0 -> gamma gamma
decay_card_bg_eta = <<~DECAYCARD
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

# J/psi -> gamma eta', eta' -> eta pi0 pi0, eta -> gamma gamma, pi0 -> gamma gamma
decay_card_bg_etap = <<~DECAYCARD
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

# ----------------------------------------------------------------------------
# Exclusive MC: the signal sample (phase-space J/psi -> gamma pi0 pi0) and the
# two exclusive background samples J/psi -> gamma eta(').
# ----------------------------------------------------------------------------
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "JpsiToGammapi0pi0"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

exMC_bg_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "JpsiToGammaEta"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_bg_eta
  config.cross_section   = :default
end

exMC_bg_etap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "JpsiToGammaEtaPrime"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_bg_etap
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiToGammapi0pi0"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.097]})
   .set_alias({"std::vector<double>" => "Vdouble"})

# Selection: no charged tracks, at least five good photons, then the kinematic
# fit of the final state gamma pi0 pi0.
event_selection = Selection.new
event_selection.select_track {
                  cos_theta   0.93   # |cos(theta)| < 0.93 in the MDC
                  Vz          10.0   # |Vz| < 10 cm along the beam direction
                  Vr          1.0    # Vr < 1 cm in the plane perpendicular to the beam
                  nChrp       "==0"  # no positively charged tracks (charged-track veto)
                  nChrn       "==0"  # no negatively charged tracks
                  nNet        "==0"  # net charge zero
                }
               .select_photon {
                  tdc_emc_start     0     # EMC timing window
                  tdc_emc_end       14
                  angle_to_track    10.0  # photon at least 10 deg from any charged track
                  energyThreshold_b 0.025 # barrel minimum EMC energy: 25 MeV
                  energyThreshold_e 0.050 # endcap minimum EMC energy: 50 MeV
                  nGam              ">=5" # at least five good photon candidates
                }
               # Reconstruct the two pi0 candidates from four of the five photons,
               # each photon pair mass-constrained to the nominal pi0 mass (1C per
               # pi0); the fifth photon is the radiated photon of the signal.
               .kalman_kinematic_fit([:gamma, :gamma, :gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # first pi0
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # second pi0
                  chi2_cut 25
                  npi0 ">=2"     # two reconstructed pi0 candidates
                }
               # Nominal kinematic fit: constrain the four-momentum of the detected
               # gamma pi0 pi0 system to the initial J/psi four-momentum. Together
               # with the two pi0 mass constraints this corresponds to the
               # published 6C kinematic fit.
               .kinematic_fit([:gamma, :pi0, :pi0]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200   # loose BOSS-pass cut; the published chi2 < 20 / < 60 is applied in ROOT
                }

# Notes for BOSS-side criteria without a dedicated DSL construct (all of them
# act on the 6C-fit-corrected quantities or on the fit permutations and are
# therefore applied at the ROOT stage)
alg
  .note(:chi2_selection, "The published selection uses a six-constraint (6C) " \
        "kinematic fit of each photon permutation to the final state gamma pi0 pi0: " \
        "a constraint on the four-momentum of the final state to that of the " \
        "initial J/psi (4C) plus an additional constraint (1C) on each photon pair " \
        "to have the invariant mass of the pi0. The chi2 requirement is " \
        "mass dependent: chi2 < 20 for M(pi0 pi0) below the KK threshold and " \
        "chi2 < 60 above it. The nominal BOSS pass uses the loose default " \
        "chi2_cut 200 and the mass-dependent values are applied in ROOT, where the " \
        "systematic effect of loosening them (to < 60 below and < 125 above the KK " \
        "threshold) is also evaluated (0.1% on the branching fraction).")
  .note(:combination_selection, "If more than one permutation of the five photons " \
        "satisfies the selection criteria, only the permutation with the minimum " \
        "chi2 from the 6C kinematic fit is retained. Additionally, the invariant " \
        "mass of any photon pair associated with a pi0 must fall within " \
        "13 MeV/c^2 of the pi0 mass before the fit; the smallest-chi2 permutation " \
        "choice is made on the fitted quantities in the ROOT analysis.")
  .note(:background_veto, "The J/psi -> omega pi0 (omega -> gamma pi0) background " \
        "is reduced by requiring the invariant mass of each gamma pi0 pair to be at " \
        "least 50 MeV/c^2 away from the nominal omega mass; an alternative treatment " \
        "that instead includes an explicit omega pi0 amplitude in the fit gives " \
        "consistent results (0.8% on the branching fraction).")
  .note(:background_veto, "The misreconstructed background, in which the radiated " \
        "photon is paired with a pi0 daughter photon to form a pi0, is reduced by " \
        "requiring the invariant mass of the radiated photon paired with any pi0 " \
        "daughter photon to be greater than 0.15 GeV/c^2.")
  .note(:background_estimation, "Backgrounds from J/psi -> gamma eta (eta -> pi0 pi0 pi0) " \
        "and J/psi -> gamma eta' (eta' -> eta pi0 pi0) are generated with exclusive " \
        "MC according to the PDG branching fractions and are included in the " \
        "unbinned extended maximum-likelihood fit with a negative weight, which " \
        "approximately cancels their residual contribution in the data. All other " \
        "backgrounds are estimated with the inclusive MC sample (about 1.5% " \
        "contamination, assigned a 100% systematic uncertainty), plus the " \
        "misreconstructed background from a dedicated exclusive MC sample that " \
        "resembles the data.")
  .note(:efficiency_curve, "The reconstruction efficiency of the selection is " \
        "determined to be 28.7%. The branching fraction " \
        "B(J/psi -> gamma pi0 pi0) = (N_gamma pi0 pi0 - N_bkg) / " \
        "(epsilon_gamma * N_J/psi) uses the acceptance-corrected signal yield, " \
        "N_bkg = 35,951, epsilon_gamma = 0.9993 for the extrapolation of the " \
        "pi0 pi0 spectrum down to a radiative photon energy of zero (+0.07% of " \
        "events), and N_J/psi = (1.311 +/- 0.011) x 10^9. Continuum backgrounds " \
        "are estimated from the 3.080 GeV data sample scaled by luminosity and by " \
        "the cross-section ratio as a function of the centre-of-mass energy.")

# Generate the BOSS algorithm package (the angular-distribution amplitude fits
# in bins of M(pi0 pi0) are performed in ROOT)
alg.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC, continuum data and the exclusive MC samples
root_files = alg.execute_on([jpsi_data, jpsi_incMC, cont_data_3080,
                             exMC_signal, exMC_bg_eta, exMC_bg_etap])
