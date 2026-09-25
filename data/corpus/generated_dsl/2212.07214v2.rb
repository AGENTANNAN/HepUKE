# =============================================================================
# BOSS-side spec: search for Lambda_c+ -> Sigma+ gamma
#   (Sigma+ -> p pi0, pi0 -> gamma gamma)
# with the oppositely charged anti-Lambda_c- reconstructed in single-tag
# hadronic modes, at sqrt(s) = 4.60 - 4.70 GeV (seven energy points).
# Tag-based analysis  ->  TagAnalysis (no Algorithm + Selection).
# =============================================================================

### Dataset description ###
# Seven energy points; sample name convention = "<BOSS version>_<CMS energy in MeV>".
data_4600 = DatasetManager.real_data.find("703_4600")   # 4.59953 GeV
data_4610 = DatasetManager.real_data.find("706_4610")   # 4.61184 GeV
data_4620 = DatasetManager.real_data.find("706_4620")   # 4.62800 GeV
data_4640 = DatasetManager.real_data.find("706_4640")   # 4.64067 GeV
data_4660 = DatasetManager.real_data.find("706_4660")   # 4.66122 GeV
data_4680 = DatasetManager.real_data.find("706_4680")   # 4.68184 GeV
data_4700 = DatasetManager.real_data.find("706_4700")   # 4.69857 GeV
data_points = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]

# Matching inclusive MC samples (same seven energy points).
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_points = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

# Decay card for the signal process e+e- -> Lambda_c+ anti-Lambda_c-
# (EvtGen names; psi(4260) is the BESIII convention for the KKMC top mother).
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000  Lambda_c+  anti-Lambda_c-   PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000  anti-p-  K+  pi-            PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000  Sigma+  gamma               PHSP;
  Enddecay

  Decay Sigma+
  1.0000  p+  pi0                     PHSP;
  Enddecay

  Decay pi0
  1.0000  gamma  gamma                PHSP;
  Enddecay

  End
DECAYCARD

# 500k-event exclusive MC at every energy point (same card, different dataset).
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lambdacbar_tag_sigma_gamma"   # -> ..._703_4600, ..._706_4610, ...
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) — TagAnalysis ###
alg_name = "LambdacTagSigmaGamma"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   # ECMS is only a fallback: the tag fit builds the CMS four-vector per run from
   # MeasuredEcmsSvc (measured beam energy + measured boost). Nominal scan centre.
   .set_constant({ "ECMS" => [:double, 4.649] })
   .with_decay_card(decay_card_signal)

# --- Tag side: single tag on the anti-Lambda_c- (charm = -1) ----------------
alg.tag_side(:Lambdac) do |t|
  # Only three tag modes are available (the other seven paper modes are not in DTagAlg).
  t.modes :LambdacPtoKPiP, :LambdacPtoKsP, :LambdacPtoLambdaPiP   # anti-p K+ pi-, anti-p K_S0, anti-Lambda pi-
  t.charm -1                                                      # pin the anti-Lambda_c- side
  # Explicit tag-side mBC window requested by the analysis (only sanctioned tag-side cut form).
  t.window :mBC, min: 2.275, max: 2.310
end

# --- Signal side: what the tag did not use ---------------------------------
alg.signal_side do |s|
  s.photons 3                       # radiative photon + the two pi0 photons
  s.charged(prp: 1)                 # exactly one additional proton
  s.require_charge(1)               # the extra charged track is a proton (p+), not an anti-proton
  s.min_photon_angle 10.0           # pi0 photons must be at least 10 deg from charged tracks
  s.min_photon_energy 0.025         # pi0 photon energy floor (25 MeV)
end

# --- Kinematic fit: 4C, four-momentum conservation, chi2 < 200 -------------
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)        # M(gamma gamma) pi0 window
  f.invariant_mass_of(:prp, :gamma, :gamma).between(1.176, 1.200)  # M(p pi0) Sigma+ window
  f.chi2_cut 200
end

# --- Inexpressible BOSS-side procedures ------------------------------------
alg.note(:tag_side_selection,
  "tag-side charged tracks are required to have |cos(theta)| < 0.93, |Vxy| < 1 cm and " \
  "|Vz| < 10 cm (except the K_S0 / Lambda daughters), with PID from dE/dx + TOF probability; " \
  "these are applied inside DTagAlg's own single-tag reconstruction and are not DSL-tunable. " \
  "K_S0 -> pi+ pi- uses |Vz| < 20 cm, secondary-vertex chi2 < 100, decay length > 2 sigma and " \
  "M(pi+ pi-) in (0.487, 0.511) GeV/c^2; Lambda -> anti-p pi+ uses M(anti-p pi+) in " \
  "(1.111, 1.121) GeV/c^2.")
alg.note(:tag_candidate_ranking,
  "when more than one tag candidate passes the mBC window, the candidate with the minimum " \
  "|DeltaE| is retained (candidate ordering on the stored single-tag collection)")
alg.note(:signal_side_selection,
  "signal-side extra proton required with |Vz| < 20 cm; the radiative photon is required to " \
  "have E > 0.65 GeV and the two pi0 photons 25 MeV < E < 0.45 GeV. The per-photon energy " \
  "windows and the proton |Vz| cut have no representation in signal_side.")
alg.note(:signal_deltaE,
  "signal-side DeltaE_sig is required to lie in (-0.038, 0.026) GeV and the combination with " \
  "the minimum |DeltaE_sig| is kept; DeltaE_sig is stored and windowed on the ROOT side")

# --- Render and run ---------------------------------------------------------
alg.apply
root_files = alg.execute_on(data_points + incMC_points + exMCs)