# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset preparation — sqrt(s) = 4.599 GeV ###
data_4600  = DatasetManager.real_data.find("703_4600")      # 567 pb^-1 of real data at 4.599 GeV
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")   # matching inclusive MC (sample 703_4600)

# Decay card for the signal process: e+e- -> Lambda_c+ anti-Lambda_c-,
# Lambda_c+ -> Lambda e+ nu_e (PHOTOS + phase space), Lambda -> p pi-.
# The anti-Lambda_c- is left to decay inclusively ("anything") through the
# EvtGen default decay table, so that the hadronic tag modes are populated.
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c-  PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Lambda e+ nu_e  PHOTOS PHSP;
  Enddecay

  Decay Lambda
  1.0000 p+ pi-  PHSP;
  Enddecay

  End
DECAYCARD

# 500k exclusive signal MC events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4600_lambdac_to_lambda_e_nu"
  config.related_dataset = data_4600
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) — tag-based analysis ###
alg_name = "LcToLambdaENu"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.599]})

# -----------------------------------------------------------------------------
# Tag side: anti-Lambda_c- reconstructed in the eleven hadronic modes
# (charge-conjugate modes of p K_S0, p K- pi+, ... -> charm -1 pins anti-Lambda_c-)
# -----------------------------------------------------------------------------
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP,          # pbar K_S0
          :LambdacPtoKPiP,         # pbar K+ pi-
          :LambdacPtoKPiPi0P,      # pbar K+ pi- pi0
          :LambdacPtoKsPi0P,       # pbar K_S0 pi0
          :LambdacPtoKsPiPiP,      # pbar K_S0 pi+ pi-
          :LambdacPtoLambdaPi,     # Lambdabar pi-
          :LambdacPtoLambdaPiPi0,  # Lambdabar pi- pi0
          :LambdacPtoLambdaPiPiPi, # Lambda pi- pi+ pi-
          :LambdacPtoSigma0Pi,     # Sigmabar0 pi-
          :LambdacPtoSigmamPi0,    # Sigma- pi0
          :LambdacPtoSigmamPiPi    # Sigma- pi+ pi-
  t.charm(-1)                    # the tagged side is the anti-Lambda_c-
  # Tag-side selection windows explicitly requested by the analysis
  # (per-mode DeltaE windows lie inside this envelope)
  t.window :deltaE, min: -0.049, max: 0.062
  t.window :mBC,    min: 2.280,  max: 2.296
end

# -----------------------------------------------------------------------------
# Signal side: Lambda_c+ -> Lambda e+ nu_e with Lambda -> p pi-
# exactly one p, one pi- and one e+ (net charge +1), no photons, missing nu_e
# -----------------------------------------------------------------------------
alg.signal_side do |s|
  s.charged(prp: 1, pim: 1, ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

# -----------------------------------------------------------------------------
# Kinematic fit: 4-momentum conservation (with the missing neutrino) plus a
# nominal-Lambda mass constraint on the p pi- pair; loose chi2 < 200 in BOSS.
# -----------------------------------------------------------------------------
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
  f.chi2_cut 200
end

# -----------------------------------------------------------------------------
# BOSS-side procedures that cannot be expressed in the DSL surface
# -----------------------------------------------------------------------------
alg.note(:track_selection,
  "charged tracks require |cos theta| < 0.93, |Vz| < 10 cm and Vr < 1 cm; " \
  "tracks originating from K_S0 / Lambda secondary vertices are exempt from the Vr/Vz cuts")

alg.note(:photon_pi0_selection,
  "photons are barrel (|cos theta| <= 0.80) EMC clusters with E > 25 MeV or endcap " \
  "(0.86 <= |cos theta| <= 0.92) clusters with E > 50 MeV, isolated by more than 10 degrees " \
  "from the nearest charged track, with EMC time in (0, 700) ns; pi0 candidates are gamma gamma " \
  "pairs with 0.110 < M(gamma gamma) < 0.155 GeV/c^2 mass-constrained with chi2 < 20")

alg.note(:pid_criteria,
  "pi/K separation from MDC dE/dx and TOF likelihoods (L_pi > L_K for pions, L_K > L_pi for kaons); " \
  "proton identification from combined dE/dx+TOF+EMC likelihoods (L'_p > L'_pi and L'_p > L'_K); " \
  "electron candidates require L'_e > 0.001 and L'_e/(L'_e+L'_pi+L'_K) > 0.8")

alg.note(:bremsstrahlung_recovery,
  "bremsstrahlung photon recovery applied to electron candidates within a 5-degree cone")

alg.note(:secondary_vertex_selection,
  "K_S0, Lambda, Sigmabar0 and Sigma- are formed with secondary-vertex fits requiring a positive " \
  "decay length; mass windows: K_S0 -> pi+pi- 0.485-0.510, Lambda -> pbar pi+ 1.110-1.121, " \
  "Sigmabar0 -> gamma Lambdabar 1.179-1.205, Sigma- -> pbar pi0 1.173-1.200 GeV/c^2")

alg.note(:background_veto,
  "pbar K_S0 pi0 tag mode vetoes M(pbar pi+) in (1.105,1.125) and M(pbar pi0) in (1.173,1.200) GeV/c^2; " \
  "Lambda pi+ pi- pi- and Sigma- pi+ pi- tag modes veto M(pi+pi-) in (0.480,0.520) and " \
  "M(pbar pi+) in (1.105,1.125) GeV/c^2")

alg.note(:tag_windows,
  "the beam-constrained-mass and DeltaE windows are applied per tag mode (the declared windows are " \
  "the envelope spanning -0.049 to +0.062 GeV in DeltaE and 2.280-2.296 GeV/c^2 in mBC); " \
  "the resulting single-tag yield is N_tot = 14415 +/- 159")

alg.note(:neutrino_kinematics,
  "the missing-neutrino kinematics is inferred from U_miss = E_miss - c|p_miss|, which peaks at zero for signal")

alg.apply

root_files = alg.execute_on([data_4600, incMC_4600, exMC_signal])