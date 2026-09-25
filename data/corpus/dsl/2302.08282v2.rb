# Measurement of the branching fraction of J/psi -> gamma eta
#   [arXiv:2302.08282]
#
# J/psi -> gamma eta with eta -> pi+ pi- pi0 and eta -> 3pi0
# Two independent decay modes -> two Algorithm objects (Rule T1)
# Single energy: sqrt(s) = 3.097 GeV, (10087 +/- 44) x 10^6 J/psi events

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ---- Decay cards ----
# Mode I: J/psi -> gamma eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_mode1 = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta   VSP_PWAVE;
  Enddecay

  Decay eta
  1.0000 pi+ pi- pi0   PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: J/psi -> gamma eta, eta -> pi0 pi0 pi0, pi0 -> gamma gamma
decay_card_mode2 = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta   VSP_PWAVE;
  Enddecay

  Decay eta
  1.0000 pi0 pi0 pi0   PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay

  End
DECAYCARD

# ---- Exclusive MC samples ----
exMC_mode1 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_gamma_eta_to_pippimpi0"
  c.related_dataset = jpsi_data
  c.events          = 1_000_000
  c.decay_card      = decay_card_mode1
  c.cross_section   = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_gamma_eta_to_3pi0"
  c.related_dataset = jpsi_data
  c.events          = 1_000_000
  c.decay_card      = decay_card_mode2
  c.cross_section   = :default
end

# =============================================================================
# ALGORITHM 1: Mode I — J/psi -> gamma eta, eta -> pi+ pi- pi0
# =============================================================================
alg1 = Algorithm.new("JpsiToGammaEta_Charged")
alg1.set_header(["JpsiToGammaEtaChargedAlg/JpsiToGammaEtaCharged.h"])
    .set_constant({"ECMS" => [:double, 3.097]})

sel1 = Selection.new

# Exactly 2 good charged tracks with net charge zero
sel1.select_track {
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     ">=1"
       nChrn     ">=1"
       nTot      2
       nNet      0
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :pip, :pim, against: [:kaon]
     }
     # At least 3 photons: radiative photon + 2 photons from pi0
     .select_photon {
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       angle_to_track    10.0
       tdc_emc_start     0
       tdc_emc_end       14
       nGam              ">=3"
     }

# Reconstruct pi0 from gamma gamma (1C Kalman fit, chi2 < 25)
sel1.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
end

# The radiative photon is the one with maximum energy among the remaining photons.
# 6C kinematic fit: 4C + 1C (pi0 mass) + 1C (eta mass), chi2_6C < 80
# Participants: gamma_radiative + pi+ + pi- + pi0    (eta = pi+ pi- pi0)
sel1.kinematic_fit([:gamma, :pip, :pim, :pi0]) do
  constrain_four_momentum
  invariant_mass_of(:pip, :pim, :pi0).constrain_to_nominal_mass_of(:eta)
  chi2_cut 80
  nominal
end

alg1.with_decay_card(decay_card_mode1).apply(sel1)

alg1.note(:radiative_photon,
        "The radiative photon is identified as the photon with the maximum energy " \
        "among those not used in the pi0 reconstruction.")
  .note(:background_veto_mode1,
        "Veto for J/psi -> 3gamma (pi0 pi0 gamma): if there are exactly 3 photons " \
        "and the event passes a 4C kinematic fit to the 3gamma hypothesis with " \
        "better probability than the 2gamma + pi+pi- hypothesis, the event is rejected. " \
        "Similarly, if there are 4 photons and the event passes a 4C fit to the 4gamma " \
        "hypothesis, it is rejected. These vetoes are applied at the ROOT level after " \
        "storing the 4C chi2 values.")
  .note(:photon_selection,
        "Photon selection: E > 25 MeV in barrel (|cos(theta)| < 0.80), " \
        "E > 50 MeV in endcap (0.86 < |cos(theta)| < 0.92). " \
        "EMC timing: 0 <= T <= 14 (in units of 50 ns). " \
        "Opening angle between photon and nearest charged track > 10 degrees.")
  .note(:background_mode1,
        "Dominant backgrounds: J/psi -> pi+ pi- pi0 (radiative Bhabhas misidentified), " \
        "J/psi -> gamma pi0 pi0 (one pi0 misreconstructed), " \
        "J/psi -> gamma pi+ pi- (non-resonant). " \
        "Continuum e+e- -> gamma pi+ pi- pi0 estimated from off-peak data. " \
        "Background shape modeled with MC and sideband data (ROOT level).")
  .note(:signal_extraction,
        "Signal yield extracted from an unbinned maximum likelihood fit to the " \
        "pi+ pi- pi0 invariant mass spectrum in the eta mass region [0.50, 0.60] GeV/c^2. " \
        "The signal is modeled with a Crystal Ball function; background with a " \
        "second-order Chebyshev polynomial. The fit yields 631,686 +/- 1,123 signal events.")
  .note(:efficiency,
        "Efficiency for Mode I: 25.57% (determined from signal MC). " \
        "Systematic uncertainties include tracking (1.0% per track), " \
        "PID (1.0% per track), photon detection (1.0% per photon), " \
        "pi0 reconstruction (1.0%), kinematic fit (1.0%), and fit model (0.5%).")

alg1.execute_on([jpsi_data, jpsi_incMC, exMC_mode1])

# =============================================================================
# ALGORITHM 2: Mode II — J/psi -> gamma eta, eta -> 3pi0
# =============================================================================
alg2 = Algorithm.new("JpsiToGammaEta_Neutral")
alg2.set_header(["JpsiToGammaEtaNeutralAlg/JpsiToGammaEtaNeutral.h"])
    .set_constant({"ECMS" => [:double, 3.097]})

sel2 = Selection.new

# All-neutral: no charged tracks
sel2.select_track {
       nTot 0
     }
     # At least 7 photons: radiative photon + 6 photons from 3 pi0
     .select_photon {
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       angle_to_track    10.0
       tdc_emc_start     0
       tdc_emc_end       14
       nGam              ">=7"
     }

# Reconstruct three pi0 from gamma gamma pairs (1C Kalman fit, chi2 < 25 each)
# |cos(theta_decay)| < 0.95 for each pi0
sel2.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0 ">=3"
end

# The radiative photon is the one with maximum energy among the remaining photons.
# 8C kinematic fit: 4C + 3 x 1C (pi0 masses) + 1C (eta mass), chi2_8C < 70
# Participants: gamma_radiative + pi0 + pi0 + pi0   (eta = 3pi0)
sel2.kinematic_fit([:gamma, :pi0, :pi0, :pi0]) do
  constrain_four_momentum
  invariant_mass_of(:pi0, :pi0, :pi0).constrain_to_nominal_mass_of(:eta)
  chi2_cut 70
  nominal
end

alg2.with_decay_card(decay_card_mode2).apply(sel2)

alg2.note(:all_neutral_mode,
        "All-neutral final state: no charged tracks required.")
  .note(:radiative_photon_mode2,
        "The radiative photon is identified as the photon with the maximum energy " \
        "among those not used in pi0 reconstruction.")
  .note(:pi0_decay_angle,
        "For each pi0 candidate, the decay angle |cos(theta_decay)| < 0.95 is required, " \
        "where theta_decay is the angle between the photon direction in the pi0 rest " \
        "frame and the pi0 boost direction. Applied at ROOT level after Kalman fit.")
  .note(:photon_selection_mode2,
        "Same photon selection as Mode I: E > 25 MeV (barrel), E > 50 MeV (endcap), " \
        "EMC timing 0-14 units, angle to nearest charged track > 10 degrees " \
        "(though there are no charged tracks in this mode).")
  .note(:background_mode2,
        "Dominant backgrounds: J/psi -> gamma pi0 pi0 (only 5 gammas, negligible), " \
        "J/psi -> 3pi0 (missing a radiative photon), " \
        "J/psi -> omega eta (omega -> gamma pi0), " \
        "e+e- -> gamma gamma (ISR + continuum). " \
        "Continuum estimated from off-peak data (ROOT level).")
  .note(:signal_extraction_mode2,
        "Signal yield extracted from an unbinned ML fit to the 3pi0 invariant mass " \
        "distribution in the eta mass region [0.50, 0.60] GeV/c^2. " \
        "The fit yields 272,322 +/- 1,033 signal events.")
  .note(:efficiency_mode2,
        "Efficiency for Mode II: 8.15% (determined from signal MC). " \
        "The lower efficiency relative to Mode I is mainly due to the requirement of " \
        "six photons from three pi0 decays. Systematic uncertainties include photon " \
        "detection (1.0% per photon), pi0 reconstruction (1.0% per pi0), " \
        "kinematic fit (1.5%), and fit model (0.5%).")

alg2.execute_on([jpsi_data, jpsi_incMC, exMC_mode2])