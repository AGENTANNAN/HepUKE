# BESIII: observation of Lambda_c+ -> n K_S0 pi+ (567 pb^-1 at sqrt(s) = 4.599 GeV)
#
# The analysis is a Lambda_c single-tag (ST) / double-tag (DT) measurement:
#   * ST sample: anti-Lambda_c- reconstructed in eleven hadronic decay modes,
#     recoiling against the Lambda_c+ that decays to the signal final state.
#   * DT sample: Lambda_c+ -> n K_S0 pi+ searched in the system recoiling
#     against the ST anti-Lambda_c-.  The neutron is undetected and its
#     kinematics is inferred from four-momentum conservation (missing mass).
#
# DTagTool reconstructs the tag side from pre-stored tag candidates; the
# signal side is built from the tracks and showers the tag did not use.  The
# signal-side content is identical for every tag mode, so a single TagAnalysis
# declares all eleven tag modes.

### Dataset description ###
data_4599  = DatasetManager.real_data.find("703_4600")    # 567 pb^-1 at sqrt(s) = 4.599 GeV
incMC_4599 = DatasetManager.inclusive_mc.find("703_4600") # generic Lc+ Lc-, D(s)*, ISR psi, QED

### Decay card for the exclusive signal MC ###
# Lambda_c+ decays only to the signal mode while the anti-Lambda_c- (the tag
# side) decays generically, since the tag modes are the eleven hadronic
# channels already handled by DTagAlg.
decay_card_nKsPi = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-                        PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  n0  K_S0  pi+                                    PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                                         PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                                         PHSP;
    Enddecay

    End
DECAYCARD

exMC_nKsPi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "LambdacToNKsPi"
  config.related_dataset = data_4599
  config.events          = 500000
  config.decay_card      = decay_card_nKsPi
  config.cross_section   = :default
end
exMC_nKsPi.save_to_config(format: :yaml, file_path: 'temp_for_test')

datasets = [data_4599, incMC_4599]

### Event selection (BOSS) — ST anti-Lambda_c- + Lambda_c+ -> n K_S0 pi+ ###
alg_name = "LcTagNKsPiMissN"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.599] })
   .set_alias({ "std::vector<double>" => "Vdouble" })
   .note(:tag_side_selection,
         "Tag anti-Lambda_c- reconstructed in eleven hadronic modes.  Charged tracks: " \
         "|cos(theta)| < 0.93, |Vz| < 10 cm, |Vr| < 1 cm; the distance requirements are " \
         "not applied to tracks originating from K_S0 / Lambda decays.  PID: pion " \
         "(L_pi > L_K), kaon (L_K > L_pi), proton (L'_p > L'_pi and L'_p > L'_K) using " \
         "dE/dx, TOF and EMC.  Photons: E > 25 MeV barrel (|cos(theta)| <= 0.80) or " \
         "E > 50 MeV endcap (0.86 <= |cos(theta)| <= 0.92), angle to the nearest charged " \
         "track > 10 deg, EMC time within (0, 700) ns.  pi0: M(gamma gamma) in " \
         "(0.110, 0.155) GeV/c^2 with a 1C mass constraint and chi2 < 20, fitted momenta " \
         "used downstream.  K_S0 / Lambda: vertex-constrained fit to pi+ pi- / pbar pi+ " \
         "with the fitted track parameters used downstream and signed decay length L > 0; " \
         "M(pi+ pi-) in (0.485, 0.510) GeV/c^2 for K_S0, M(pbar pi+) in " \
         "(1.110, 1.121) GeV/c^2 for Lambda, M(gamma Lambdabar) in (1.179, 1.205) GeV/c^2 " \
         "for Sigma0 and M(pbar pi0) in (1.173, 1.200) GeV/c^2 for Sigma-.")
   .note(:tag_background_veto,
         "For the tag mode pbar K_S0 pi0 the Lambda / Sigma- backgrounds are rejected by " \
         "vetoing M(pbar pi+) in (1.105, 1.125) GeV/c^2 and M(pbar pi0) in " \
         "(1.173, 1.200) GeV/c^2.  For the tag modes Lambda pi+ pi- pi- and " \
         "Sigma- pi+ pi- pi+ the K_S0 / Lambda backgrounds are suppressed by requiring " \
         "M(pi+ pi-) not in (0.480, 0.520) GeV/c^2 and M(pbar pi+) not in " \
         "(1.105, 1.125) GeV/c^2.")
   .note(:tag_mbc_deltae,
         "Tag candidates are identified with M_BC = sqrt(E_beam^2/c^4 - |p_Lambdabar_c-|^2/c^2) " \
         "and Delta E = E_beam - E_Lambdabar_c-, each mode keeping the approximately " \
         "3 sigma_DeltaE window around the Delta E peak (mode-dependent bounds spanning " \
         "-0.025..0.062 GeV).  The per-mode tag yields are extracted from fits to M_BC in " \
         "the signal region (2.280, 2.296) GeV/c^2, giving a total ST yield " \
         "N_tot = 14415 +- 159 summed over the eleven modes.  Stored and applied at the " \
         "ROOT stage.")
   .note(:signal_ks0_selection,
         "The signal-side K_S0 uses the same criteria as the tag side (vertex-constrained " \
         "fit with the fitted track parameters used downstream and L > 0) but without the " \
         "M(pi+ pi-) mass requirement; when several candidates are formed the one with the " \
         "largest decay-length significance L/sigma_L is retained.  The signal pion has " \
         "charge opposite to the ST Lambdabar_c-.")
   .note(:missing_neutron_observable,
         "The neutron is undetected: M_miss^2 = E_miss^2/c^4 - |p_miss|^2/c^2 is formed " \
         "from E_miss = E_beam - E_K_S0 - E_pi+ and p_miss = p_Lambda_c+ - p_K_S0 - p_pi+, " \
         "with p_Lambda_c+ = -p_hat_tag sqrt(E_beam^2/c^2 - m^2_Lambdabar_c- c^2) using the " \
         "ST Lambda_c- direction and the nominal Lambda_c- mass.  Signal peaks at the " \
         "nominal neutron mass squared; signal extraction is a simultaneous 2D unbinned " \
         "maximum-likelihood fit to (M_miss^2, M(pi+ pi-)) in the M_BC signal and sideband " \
         "regions, giving N_obs = 83.2 +- 10.6.  Handled in the ROOT analysis.")
   .note(:background_veto,
         "Dominant backgrounds Lambda_c+ -> Sigma- pi+ pi+ and Lambda_c+ -> Sigma+ pi+ pi- " \
         "with Sigma+- -> n pi+- peak in M_miss^2 but are flat in M(pi+ pi-); non-Lambda_c+ " \
         "backgrounds are estimated from the ST candidates in the M_BC sideband " \
         "(2.252, 2.272) GeV/c^2.")
   .note(:signal_mc_model,
         "The signal is generated with a phase-space model since the M(n pi+), M(n K_S0) " \
         "and M(K_S0 pi+) spectra show no obvious structure in data; the 1.3% model " \
         "uncertainty is estimated from the statistical variations of these spectra.")
   .note(:absolute_branching_fraction,
         "B(Lambda_c+ -> n K_S0 pi+) = N_obs / (N_tot x epsilon x B(K_S0 -> pi+ pi-)) with " \
         "epsilon = (45.9 +- 0.3)% (DT/ST efficiency ratio weighted over the tag modes, " \
         "excluding K_S0 -> pi+ pi-); result (1.82 +- 0.23 +- 0.11)%.")

# Tag side: the eleven hadronic anti-Lambda_c- modes.  charm -1 pins the tagged
# anti-Lambda_c- baryon (the Lambda_c+ carries the signal n K_S0 pi+).
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP,           # p K_S0
          :LambdacPtoKPiP,          # p K- pi+
          :LambdacPtoKsPi0P,        # p K_S0 pi0
          :LambdacPtoKsPiPiP,       # p K_S0 pi+ pi-
          :LambdacPtoKPiPi0P,       # p K- pi+ pi0
          :LambdacPtoPiPiP,         # p pi+ pi-
          :LambdacPtoLambdaPi,      # Lambda pi+
          :LambdacPtoLambdaPiPi0,   # Lambda pi+ pi0
          :LambdacPtoLambdaPiPiPi,  # Lambda pi+ pi+ pi-
          :LambdacPtoLambdaPiOmega, # Lambda pi+ omega
          :LambdacPtoPiSIGMA0LambdaGam # Sigma0 pi+ (Sigma0 -> gamma Lambda)
  t.charm -1
end

# Signal side: the K_S0 -> pi+ pi- pair (reconstructed as a V0) plus the signal
# pion, and the undetected neutron.
alg.signal_side do |s|
  s.photons 0
  s.charged(pip: 2, pim: 1)
  s.require_charge 1
  s.missing :n0
end

# Kinematic fit: 4-momentum conservation including the missing neutron, with
# the K_S0 mass constraint over its two pion daughters.
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut 200
end

alg.with_decay_card(decay_card_nKsPi).apply
alg.execute_on(datasets + [exMC_nKsPi])
