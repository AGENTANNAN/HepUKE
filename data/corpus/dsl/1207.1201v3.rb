# BESIII paper arXiv:1207.1201v3
# First observation of the isospin violating decay J/psi -> Lambda Sigma0bar + c.c.,
# search for Lambda(1520) -> gamma Lambda and measurement of eta_c -> Lambda Lambdabar
# in J/psi -> gamma eta_c, based on (225.2 +/- 2.8) x 10^6 J/psi events.
# Final state: gamma Lambda Lambdabar -> gamma p pi- anti-p pi+ (two photons are not
# required; only one photon from Sigma0 -> gamma Lambda or from the radiative J/psi decay).
# Nominal fit: 4C kinematic fit to the gamma Lambda Lambdabar hypothesis.

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi data, 225.2 M events
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # J/psi inclusive MC (225 x 10^6 events, kkmc + evtgen + lundcharm)

### Decay cards (EvtGen format) ###
# Signal: J/psi -> Lambda Sigma0bar (Sigma0bar -> gamma Lambdabar), Lambdabar -> anti-p- pi+
# The J/psi -> Lambda Sigma0bar + c.c. events are generated with an angular distribution
# 1 + alpha cos^2(theta), where theta is the polar angle of the baryon in the J/psi rest
# frame and alpha = 0.38 as extracted from fits to the data.
decay_card_lambda_sigma0bar = <<~DECAYCARD
    Decay J/psi
    1.0000  Lambda0  anti-Sigma0        PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000  gamma  anti-Lambda0         PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                     HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-  pi+                HypWK;
    Enddecay

    End
DECAYCARD

# Signal (charge conjugate): J/psi -> Lambdabar Sigma0 (Sigma0 -> gamma Lambda), Lambda -> p+ pi-
decay_card_lambdabar_sigma0 = <<~DECAYCARD
    Decay J/psi
    1.0000  anti-Lambda0  Sigma0        PHSP;
    Enddecay

    Decay Sigma0
    1.0000  gamma  Lambda0              PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                     HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-  pi+                HypWK;
    Enddecay

    End
DECAYCARD

# Signal: J/psi -> gamma eta_c, eta_c -> Lambda Lambdabar
# The J/psi -> gamma eta_c decays are generated with a 1 + cos^2(theta_gamma) angular
# distribution and a phase-space distribution for eta_c -> Lambda Lambdabar (no spin
# correlation effects are considered).
decay_card_gamma_etac = <<~DECAYCARD
    Decay J/psi
    1.0000  gamma  eta_c                PHSP;
    Enddecay

    Decay eta_c
    1.0000  Lambda0  anti-Lambda0       PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                     HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-  pi+                HypWK;
    Enddecay

    End
DECAYCARD

# Signal (search): J/psi -> Lambda anti-Lambda(1520) + c.c., with Lambda(1520) -> gamma Lambda.
# The detection efficiency for this channel is obtained with a phase-space MC simulation.
decay_card_lambda_1520_a = <<~DECAYCARD
    Decay J/psi
    1.0000  Lambda0  anti-Lambda(1520)0  PHSP;
    Enddecay

    Decay anti-Lambda(1520)0
    1.0000  gamma  anti-Lambda0          PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                      HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-  pi+                 HypWK;
    Enddecay

    End
DECAYCARD

decay_card_lambda_1520_b = <<~DECAYCARD
    Decay J/psi
    1.0000  Lambda(1520)0  anti-Lambda0  PHSP;
    Enddecay

    Decay Lambda(1520)0
    1.0000  gamma  Lambda0               PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                      HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-  pi+                 HypWK;
    Enddecay

    End
DECAYCARD

# Background (studied with the inclusive MC sample): the primary background sources are
# J/psi -> Lambda Lambdabar, J/psi -> Sigma0 Sigma0bar and J/psi -> Lambda Lambdabar pi0,
# where either an unrelated EMC cluster is misidentified as a photon or one of the photons
# from the Sigma0 / pi0 decay is undetected.
decay_card_bkg_lambda_lambdabar = <<~DECAYCARD
    Decay J/psi
    1.0000  Lambda0  anti-Lambda0       PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                     HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-  pi+                HypWK;
    Enddecay

    End
DECAYCARD

decay_card_bkg_sigma0_sigma0bar = <<~DECAYCARD
    Decay J/psi
    1.0000  Sigma0  anti-Sigma0         PHSP;
    Enddecay

    Decay Sigma0
    1.0000  gamma  Lambda0              PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000  gamma  anti-Lambda0         PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                     HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-  pi+                HypWK;
    Enddecay

    End
DECAYCARD

decay_card_bkg_lambda_lambdabar_pi0 = <<~DECAYCARD
    Decay J/psi
    1.0000  Lambda0  anti-Lambda0  pi0   PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                     HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-  pi+                HypWK;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples ###
exMC_lambda_sigma0bar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_lambda_sigma0bar"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_lambda_sigma0bar
  config.cross_section   = :default
end
exMC_lambda_sigma0bar.save_to_config(format: :yaml, file_path: 'exMC_jpsi_lambda_sigma0bar')

exMC_lambdabar_sigma0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_lambdabar_sigma0"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_lambdabar_sigma0
  config.cross_section   = :default
end
exMC_lambdabar_sigma0.save_to_config(format: :yaml, file_path: 'exMC_jpsi_lambdabar_sigma0')

exMC_gamma_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_gamma_etac_etac_to_llbar"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_gamma_etac
  config.cross_section   = :default
end
exMC_gamma_etac.save_to_config(format: :yaml, file_path: 'exMC_jpsi_gamma_etac')

exMC_lambda_1520_a = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_lambda_antilambda1520"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_lambda_1520_a
  config.cross_section   = :default
end

exMC_lambda_1520_b = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_lambda1520_antilambda"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_lambda_1520_b
  config.cross_section   = :default
end

exMC_bkg_lambda_lambdabar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_lambda_lambdabar_bkg"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_lambda_lambdabar
  config.cross_section   = :default
end

exMC_bkg_sigma0_sigma0bar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_sigma0_sigma0bar_bkg"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_sigma0_sigma0bar
  config.cross_section   = :default
end

exMC_bkg_lambda_lambdabar_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_lambda_lambdabar_pi0_bkg"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_lambda_lambdabar_pi0
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# All four signal channels (Lambda Sigma0bar, Lambdabar Sigma0, gamma eta_c with
# eta_c -> Lambda Lambdabar, and Lambda Lambdabar(1520) + c.c.) share the identical final
# state gamma Lambda Lambdabar and the identical selection chain below, so a single
# Algorithm instance is used for all of them.
alg_name = "JpsiToGammaLambdaLambdaBar"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.097]})
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93     # track directions within the MDC fiducial volume
                  Vz        20.0     # within +/- 20 cm of the IP along the beam direction
                  Vr        10.0     # within 10 cm of the IP in the plane perpendicular to the beam
                  nChrp     "==2"    # four charged tracks with net charge zero
                  nChrn     "==2"
                  nNet      "==0"
                }
               .select_photon {
                  tdc_emc_start     0      # EMC cluster timing suppresses electronic noise and
                  tdc_emc_end      14      # energy deposits unrelated to the event
                  angle_to_track    5.0    # photon direction at least 5 degrees from the nearest
                                           # proton / charged pion track
                  energyThreshold_b 0.025  # E_min = 25 MeV for barrel showers (|cos theta| < 0.80)
                  energyThreshold_e 0.050  # E_min = 50 MeV for end-cap showers (0.86 < |cos theta| < 0.92)
                  nGam             ">=1"   # at least one photon candidate (radiative photon)
                }
               .select_isolated_photon {
                  angle_to_prp_track  5.0   # photon direction at least 5 degrees from the nearest proton
                  angle_to_prm_track 30.0   # and at least 30 degrees from the nearest anti-proton,
                                            # since more EMC showers are found near the anti-proton
                  nGam               ">=1"
                }
               # Lambda -> p pi- from all combinations of positive and negative charged
               # track pairs; the vertex-finding algorithm must succeed.
               .secondary_vertex_fit([:prp, :pim]) {
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               # Lambdabar -> anti-p- pi+ reconstructed from the remaining tracks.
               .secondary_vertex_fit([:prm, :pip]) {
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               # 4C energy-momentum conservation kinematic fit to the gamma Lambda Lambdabar
               # hypothesis; for events with more than one photon candidate the combination
               # with the minimum chi2_4C is selected automatically.
               .kinematic_fit([:gamma, :Lambda, :Lambda_bar]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200          # loose BOSS cut; the paper requires chi2_4C < 45 (applied in ROOT)
               }

alg.note(:lambda_mass_window,
         "The Lambda and Lambdabar candidates are selected with " \
         "|M(p pi-) - M(Lambda)| < 5 MeV/c^2 and |M(anti-p- pi+) - M(Lambda_bar)| < 5 MeV/c^2. " \
         "When more than one accepted (p pi-)(anti-p- pi+) combination exists in an event, the " \
         "candidate minimizing (M(p pi-) - M(Lambda))^2 + (M(anti-p- pi+) - M(Lambda_bar))^2 is " \
         "chosen. These windows are applied on kinematic-fit-corrected four-momenta and therefore " \
         "belong to the ROOT-level analysis.")
   .note(:low_momentum_proton_veto,
         "Events containing any proton or anti-proton track candidate with momentum below " \
         "0.3 GeV/c are rejected, because the detection efficiencies for low-momentum (anti-)protons " \
         "differ between data and MC simulation. This per-track momentum cut is not expressible " \
         "through the select_track keywords.")
   .note(:sigma0_signal_region,
         "The Sigma0 (Sigma0bar) signal region is defined as within +/- 3 sigma of the nominal " \
         "Sigma0 mass, using the M(gamma Lambda) spectrum. An unbinned maximum-likelihood fit with a " \
         "double-Gaussian signal shape (Gaussian widths floating, other parameters from MC), " \
         "MC-fixed J/psi -> Lambda Lambdabar and J/psi -> Sigma0 Sigma0bar backgrounds and a " \
         "floatable second-order polynomial for the remaining background is performed over " \
         "1.165-1.30 GeV/c^2. Fitted yields: 308 +/- 24 (J/psi -> Lambdabar Sigma0) and " \
         "234 +/- 21 (J/psi -> Lambda Sigma0bar); expected backgrounds 105 +/- 10 and 95 +/- 9 events; " \
         "efficiencies 21.7% and 17.6%. This fit is a ROOT-level procedure.")
   .note(:lambda_1520_search,
         "For the search for Lambda(1520) -> gamma Lambda, the requirement M(Lambda Lambdabar) < " \
         "2.9 GeV/c^2 suppresses combinatorial backgrounds from J/psi -> Lambda Lambdabar, " \
         "Sigma0 Sigma0bar, Lambda Lambdabar pi0 and J/psi -> gamma eta_c (eta_c -> Lambda Lambdabar). " \
         "The combined M(gamma Lambda) and M(gamma Lambdabar) spectrum is fitted over " \
         "1.35-1.70 GeV/c^2 with a Breit-Wigner signal shape convolved with a double-Gaussian " \
         "resolution function and a second-order polynomial background, giving a yield of " \
         "31 +/- 24 events and an upper limit of 62.5 signal events at the 90% confidence level " \
         "(efficiency 18.8%). Both the mass requirement and the fit use post-fit quantities.")
   .note(:etac_fit,
         "For the eta_c analysis, the dominant remaining backgrounds after event selection are " \
         "incoherent contributions from J/psi -> Sigma0 Sigma0bar and J/psi -> Lambda Sigmabar0 + c.c. " \
         "(637 +/- 52 expected events in the signal region, fixed to the MC shape and level) plus an " \
         "irreducible nonresonant J/psi -> gamma Lambda Lambdabar background described by a " \
         "second-order polynomial with floating yield and shape. The eta_c line shape is " \
         "(E_gamma^3 x BW(m) x damping(E_gamma)) convolved with a Gaussian resolution, with the KEDR " \
         "damping function as default and the CLEO form as an alternative; the eta_c mass and width " \
         "are fixed to the BESIII values M = 2984.3 +/- 0.8 MeV/c^2 and Gamma = 32.0 +/- 1.6 MeV. " \
         "The fit over 2.76-3.06 GeV/c^2 gives 360 +/- 38 signal events with an efficiency of 19.8%.")
   .note(:angular_distribution_alpha,
         "The angular distribution of the baryon in J/psi -> B8 B8bar is expected to follow " \
         "1 + alpha cos^2(theta). A simultaneous fit to the efficiency-corrected cos(theta) " \
         "distributions of Lambda and Lambdabar gives alpha = 0.38 +/- 0.39; the signal MC is " \
         "generated with alpha = 0.38.")
   .note(:systematic_uncertainties,
         "Systematic uncertainties on the branching fractions (%): photon detection 1, tracking 4, " \
         "Lambda and Lambdabar vertex fits 2, 4C kinematic fit 2.3, signal shape 1.3 / 2.6 / 4.8 / 7.6 " \
         "(Lambda Sigmabar0 / Lambdabar Sigma0 / Lambda(1520) / eta_c), fitting range 1.6 / 0.9 / 1.4 / 1.4, " \
         "alpha 5.5 / 5.1 / 10.2 / -, fixed backgrounds 0.6 / 0.4 / - / 12.8, nonresonant background " \
         "shape 0.3 / 0.1 / 1.9 / 1.7, QED correction factor 0.1 / 0.1 / - / -, cited branching " \
         "fractions 0.8, number of J/psi events 1.3; total 8.0 / 7.9 / 12.6 / 16.0.")
   .note(:qed_continuum_correction,
         "The relative ratio of the QED background from e+e- -> gamma* -> Lambda Sigmabar0 + c.c. is " \
         "estimated to be (5.4 +/- 0.1)% of the measured J/psi -> Lambda Sigmabar0 + c.c. yield, so " \
         "the result is scaled by a factor 0.946; 0.1% is assigned as the systematic uncertainty of " \
         "this correction.")

alg.with_decay_card(decay_card_lambda_sigma0bar).apply(event_selection)

root_files = alg.execute_on([jpsi_data, jpsi_incMC, exMC_lambda_sigma0bar, exMC_lambdabar_sigma0,
                             exMC_gamma_etac, exMC_lambda_1520_a, exMC_lambda_1520_b,
                             exMC_bkg_lambda_lambdabar, exMC_bkg_sigma0_sigma0bar,
                             exMC_bkg_lambda_lambdabar_pi0])
