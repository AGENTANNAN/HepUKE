# =============================================================================
# BESIII arXiv:1510.02610v1
# First absolute measurement of B(Lambda_c+ -> Lambda e+ nu_e).
#
# Data: 567 pb^-1 at sqrt(s) = 4.599 GeV (sample 703_4600), just above the
# Lambda_c+ Lambda_c- threshold.
#
# Method: single-tag (ST) / double-tag (DT) technique with the tagged-D
# approach of the Mark III collaboration.  The anti-Lambda_c- is reconstructed
# in eleven hadronic decay modes (ST sample) and the semileptonic decay
# Lambda_c+ -> Lambda e+ nu_e is searched for in the system recoiling against
# it.  The neutrino is undetected; its kinematics is inferred from
# four-momentum conservation through U_miss = E_miss - c|p_miss|, which peaks
# at zero for signal.  The absolute branching fraction follows from
#   B(Lambda_c+ -> Lambda e+ nu_e) =
#       N_semi / (N_tot x epsilon_semi x B(Lambda -> p pi-)),
# so the ST-side uncertainties largely cancel.
#
# Charge-conjugate processes are included throughout.
# This is the electron channel of the same measurement whose muon channel is
# arXiv:1611.04382v1; the tag side and the analysis structure are identical.
# =============================================================================

### Dataset description ###
data_4600  = DatasetManager.real_data.find("703_4600")    # sqrt(s) = 4.599 GeV, 567 pb^-1
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600") # Lambda_c+ Lambda_c-, D(s), ISR psi, continuum, QED

### Decay card for the exclusive signal MC ###
# The Lambda_c+ decays only to Lambda e+ nu_e and the anti-Lambda_c- decays to
# the tag modes.  The signal is generated with the form-factor predictions of
# Heavy Quark Effective Theory and QCD sum rules (Ref. [13] of the paper); the
# form-factor shape is not expressible in EvtGen syntax here, so phase space is
# used and the model dependence is recorded as a systematic uncertainty through
# a note.  ISR and FSR effects are included in the production and decay.
decay_card_leptonic = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-                        PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  Lambda0  e+  nu_e                                PHOTOS PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                                         PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                                          PHSP;
    Enddecay

    End
DECAYCARD

exMC_leptonic = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "LambdacToLambdaENu"
  config.related_dataset = data_4600
  config.events          = 500_000
  config.decay_card      = decay_card_leptonic
  config.cross_section   = :default
end
exMC_leptonic.save_to_config(format: :yaml, file_path: 'temp_for_test')

datasets = [data_4600, incMC_4600]

### Event selection (BOSS) — ST anti-Lambda_c- + Lambda_c+ -> Lambda e+ nu_e ###
alg_name = "LcTagLambdaENu"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.599] })
   .set_alias({ "std::vector<double>" => "Vdouble" })
   .note(:tag_side_selection,
         "The anti-Lambda_c- is reconstructed in the eleven hadronic modes listed " \
         "below, from the pre-stored tag collection.  Charged tracks require " \
         "|cos(theta)| < 0.93 with respect to the beam direction and a distance of " \
         "closest approach to the interaction point below 10 cm along the beam and " \
         "below 1 cm in the perpendicular plane; tracks originating from K_S0 and " \
         "Lambda decays are exempt from the distance requirements.  Pion/kaon " \
         "separation combines the MDC dE/dx and TOF information into likelihoods " \
         "L_pi and L_K, with L_pi > L_K for pions and L_K > L_pi for kaons; proton " \
         "identification combines dE/dx, TOF and EMC into L'_p and requires " \
         "L'_p > L'_pi and L'_p > L'_K.")
   .note(:tag_photon_pi0_selection,
         "Photon candidates are isolated EMC clusters in the barrel " \
         "(|cos(theta)| <= 0.80) or end cap (0.86 <= |cos(theta)| <= 0.92) regions " \
         "with deposited energy above 25 MeV (barrel) or 50 MeV (end cap), more than " \
         "10 deg from the nearest charged track and with the EMC time minus the event " \
         "start time within (0, 700) ns.  pi0 candidates are photon pairs with " \
         "0.110 < M(gamma gamma) < 0.155 GeV/c^2, mass-constrained to the nominal " \
         "pi0 mass with chi2 < 20; the fitted pi0 momenta are used downstream.")
   .note(:tag_v0_selection,
         "K_S0 and Lambda candidates are rebuilt with a secondary vertex fit and are " \
         "required to have a positive decay length.  The invariant-mass windows are " \
         "0.485 < M(pi+ pi-) < 0.510 GeV/c^2 for K_S0, " \
         "1.110 < M(pbar pi+) < 1.121 GeV/c^2 for Lambda, " \
         "1.179 < M(gamma Lambda) < 1.205 GeV/c^2 for Sigmabar0 and " \
         "1.173 < M(pbar pi0) < 1.200 GeV/c^2 for Sigma-.  Sigmabar0 -> gammabar " \
         "Lambda and Sigma- -> pbar pi0 / pbar pi- pi0 are reconstructed through " \
         "these intermediate states.")
   .note(:tag_background_veto,
         "For the pbar K_S0 pi0 tag mode, Lambda and Sigmabar- backgrounds are " \
         "rejected by vetoing events with M(pbar pi+) inside (1.105, 1.125) GeV/c^2 " \
         "or M(pbar pi0) inside (1.173, 1.200) GeV/c^2.  For the Lambda pi+ pi- pi- " \
         "and Sigmabar- pi+ pi- modes, K_S0 backgrounds are suppressed by requiring " \
         "M(pi+ pi-) outside (0.480, 0.520) GeV/c^2, and Lambda backgrounds are " \
         "removed by requiring M(pbar pi+) outside (1.105, 1.125) GeV/c^2.")
   .note(:tag_mbc_deltae,
         "Tag candidates are identified with the beam-constrained mass " \
         "M_BC = sqrt(E_beam^2 - |p_anti-Lambda_c-|^2) and the energy difference " \
         "Delta E = E_beam - E_anti-Lambda_c-.  Each mode keeps the approximately " \
         "3 sigma_DeltaE window around its Delta E peak; the per-mode bounds from " \
         "Table I of the paper are, in the order of the mode list below, " \
         "[-0.025, 0.028], [-0.019, 0.023], [-0.035, 0.049], [-0.044, 0.052], " \
         "[-0.029, 0.032], [-0.033, 0.035], [-0.037, 0.052], [-0.028, 0.030], " \
         "[-0.029, 0.032], [-0.038, 0.062] and [-0.049, 0.054] GeV.  The per-mode " \
         "yields are extracted from unbinned maximum-likelihood fits to the M_BC " \
         "spectra (MC signal shape convoluted with a double-Gaussian resolution " \
         "function, ARGUS background) in the signal region " \
         "(2.280, 2.296) GeV/c^2; summed over the eleven modes the total ST yield is " \
         "N_tot = 14415 +- 159.  M_BC and Delta E are stored and windowed at the ROOT " \
         "stage.")
   .note(:signal_lambda_selection,
         "The signal Lambda is formed from a p pi- pair recoiling against the ST " \
         "anti-Lambda_c-, using the same secondary-vertex and mass-window criteria " \
         "as the ST selection.  The charged track identified as the positron must be " \
         "the only other unused track in the event.")
   .note(:signal_electron_pid,
         "Electron identification uses the probabilities L'_e computed from the MDC " \
         "dE/dx, the TOF and the EMC, requiring L'_e > 0.001 and " \
         "L'_e/(L'_e + L'_pi + L'_K) > 0.8.  The energy loss due to bremsstrahlung " \
         "photons is partially recovered by adding showers within a 5 deg cone about " \
         "the positron momentum.  The tag vocabulary classifies an `ep` track through " \
         "SimplePIDSvc with fixed v1 lepton thresholds, which are not DSL-tunable; " \
         "the working-point difference is absorbed in the electron PID systematic " \
         "uncertainty (1.0%, studied with e+e- -> (gamma) e+ e-).")
   .note(:missing_neutrino_observable,
         "The neutrino is undetected.  U_miss = E_miss - c|p_miss| is formed from " \
         "E_miss = E_beam - E_Lambda - E_e+ and " \
         "p_miss = p_Lambda_c+ - p_Lambda - p_e+, where the Lambda_c+ momentum is " \
         "p_Lambda_c+ = -p_hat_tag sqrt(E_beam^2 - m^2_anti-Lambda_c-) using the ST " \
         "anti-Lambda_c- direction and its nominal mass.  Signal events peak at " \
         "U_miss = 0.  The yield is extracted in ROOT from a fit to the U_miss " \
         "distribution with a Gaussian core plus two power-law tails (tail " \
         "parameters fixed from signal MC) plus a polynomial background, giving " \
         "N_semi = 109.4 +- 10.9 before background subtraction and " \
         "N_semi = 103.5 +- 10.9 after subtracting the two peaking backgrounds " \
         "(1.4 +- 0.8 from non-Lambda semileptonic decays, estimated from the Lambda " \
         "sideband, and 4.5 +- 0.5 from Lambda mu+ nu_mu and hadronic decays such as " \
         "Lambda pi+ pi0, Lambda pi+ and Sigma0 pi+ from MC).")
   .note(:signal_mc_model,
         "The Lambda_c+ -> Lambda e+ nu_e signal is generated with the form-factor " \
         "predictions obtained from Heavy Quark Effective Theory and QCD sum rules; " \
         "varying the form-factor parameterisation and accounting for the q^2 " \
         "dependence observed in data gives a 4.5% model uncertainty.")
   .note(:absolute_branching_fraction,
         "B(Lambda_c+ -> Lambda e+ nu_e) = " \
         "N_semi / (N_tot x epsilon_semi x B(Lambda -> p pi-)) with the overall " \
         "efficiency epsilon_semi = (30.92 +- 0.26)% (weighted by the ST yields of " \
         "data for each tag, excluding B(Lambda -> p pi-)).  The result is " \
         "(3.63 +- 0.38 (stat.) +- 0.20 (syst.))%, with a total systematic " \
         "uncertainty of 5.6% dominated by the Lambda reconstruction efficiency " \
         "(2.5%) and the signal model (4.5%).")
   .note(:results,
         "B(Lambda_c+ -> Lambda e+ nu_e) = (3.63 +- 0.38 +- 0.20)% -- the first " \
         "absolute measurement, more than twofold more precise than the previous " \
         "world average, and the benchmark for all other Lambda_c+ semileptonic " \
         "channels.")

# Tag side: the eleven hadronic anti-Lambda_c- modes.
#   pbar K_S0, pbar K+ pi-, pbar K+ pi- pi0, pbar K_S0 pi0, pbar K_S0 pi+ pi-,
#   Lambdabar pi-, Lambdabar pi- pi0, Lambda pi- pi+ pi-, Sigmabar0 pi-,
#   Sigma- pi0, Sigma- pi+ pi-.
# charm -1 pins the tagged anti-Lambda_c-.
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP,          # anti-Lambda_c- -> pbar K_S0
          :LambdacPtoKPiP,         # anti-Lambda_c- -> pbar K+ pi-
          :LambdacPtoKsPi0P,       # anti-Lambda_c- -> pbar K_S0 pi0
          :LambdacPtoKPiPi0P,      # anti-Lambda_c- -> pbar K+ pi- pi0
          :LambdacPtoKsPiPiP,      # anti-Lambda_c- -> pbar K_S0 pi+ pi-
          :LambdacPtoLambdaPi,     # anti-Lambda_c- -> Lambdabar pi-
          :LambdacPtoLambdaPiPi0,  # anti-Lambda_c- -> Lambdabar pi- pi0
          :LambdacPtoLambdaPiPiPi, # anti-Lambda_c- -> Lambda pi- pi+ pi-
          :LambdacPtoSigma0Pi,     # anti-Lambda_c- -> Sigmabar0 pi-
          :LambdacPtoSigmaPPi0,    # anti-Lambda_c- -> Sigma- pi0
          :LambdacPtoSigmaPPiPi    # anti-Lambda_c- -> Sigma- pi+ pi-
  t.charm -1
end

# Signal side: the Lambda -> p pi- pair recoiling against the tag, plus the
# positron, plus the undetected neutrino.
alg.signal_side do |s|
  s.photons 0
  s.charged(prp: 1, pim: 1, ep: 1)
  s.require_charge 1                 # +1 (p) - 1 (pi-) + 1 (e+) = +1
  s.min_photon_angle 10.0
  s.missing :nu_e                    # massless (semileptonic)
end

# Kinematic fit: 4-momentum conservation including the missing neutrino, with
# the Lambda mass constraint over its two daughters.
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
  f.chi2_cut 200
end

alg.with_decay_card(decay_card_leptonic).apply
alg.execute_on(datasets + [exMC_leptonic])
