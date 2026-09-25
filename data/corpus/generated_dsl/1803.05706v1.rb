# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# 4.6 GeV real data (~567 pb^-1) and the corresponding inclusive MC sample (BOSS 7.0.3, 4599.53 MeV)
data_4600  = DatasetManager.real_data.find("703_4600")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

# Decay card for the exclusive signal MC (EvtGen format):
#   psi(4260) -> Lambda_c+ anti-Lambda_c-,
#   Lambda_c+ -> Lambda pi+ pi-,  Lambda -> p+ pi-,
#   anti-Lambda_c- -> anti-p- K_S0,  K_S0 -> pi+ pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000   Lambda_c+   anti-Lambda_c-       PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000   Lambda0   pi+   pi-              PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000   anti-p-   K_S0                   PHSP;
    Enddecay

    Decay Lambda0
    1.0000   p+   pi-                         PHSP;
    Enddecay

    Decay K_S0
    1.0000   pi+   pi-                        PHSP;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive signal MC associated with the 4.6 GeV data set
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4600_lambdac_to_lambdax"
  config.related_dataset = data_4600
  config.events          = 200_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: "temp_for_test")

### Event selection / tag analysis (BOSS) ###
alg_name = "LambdacTagIncLambda"
tag_alg = TagAnalysis.new(alg_name)
tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 4.600]})
       .with_decay_card(decay_card_signal)

# Tag side: the anti-Lambda_c- is taken from the pre-stored tag candidates in two hadronic
# modes, anti-p K_S0 and anti-p K+ pi-.  One tag_side call = single tag: the other
# Lambda_c+ is measured inclusively on the signal side (double-tag method).
tag_alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP, :LambdacPtoKPiP   # anti-p K_S0 and anti-p K+ pi-
  t.charm -1                                # pin the tagged side to the anti-Lambda_c-
  t.window :deltaE, abs: 0.02               # deltaE window (mode-dependent +-2.5 sigma, see note)
end

# Signal side: inclusive Lambda -> p+ pi- built from the tracks the tag did not use
tag_alg.signal_side do |s|
  s.charged(prp: 1, pim: 1, at_least: true)  # Lambda -> p pi-; extra (inclusive) tracks allowed
  s.require_charge(0)                        # net charge of the signal side is zero
end

# Kinematic fit: four-momentum conservation with a loose chi2 cut (tight cut applied in ROOT)
tag_alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:prp, :pim).between(1.111, 1.121)   # Lambda mass window [1111, 1121] MeV/c^2
  f.chi2_cut 200
end

# ---- BOSS-side procedures that the tag layer cannot express ----
tag_alg
  .note(:tag_track_quality, "charged tracks used by the tag (DTagAlg) are required to satisfy
    |cos(theta)| < 0.93, |Vz| < 10 cm and Vr < 1 cm; these track-quality cuts live inside the
    pre-stored tag reconstruction and are not tunable from the tag DSL")
  .note(:tag_pid_criteria, "likelihood-based PID: protons required to satisfy L(p) > L(K) and
    L(p) > L(pi), kaons required to satisfy L(K) > L(pi); no PID is applied to the pions from
    Lambda or K_S0 decays. Enforced inside the tag-side reconstruction (DTagAlg), hence not
    expressible in the tag DSL")
  .note(:deltaE_window, "deltaE windows are mode dependent and defined at +-2.5 sigma
    (sigma ~ 8 MeV, i.e. different widths for the anti-p K_S0 and anti-p K+ pi- tag modes);
    the single symmetric window declared on the tag side is only an approximation, the final
    per-mode windows are applied in ROOT")
  .note(:ks_vertex_fit, "K_S0 -> pi+ pi- reconstructed with a vertex-fit chi2 < 100 and a
    secondary-vertex constraint pointing back to the IP, requiring a decay length larger than
    twice its resolution; only one track pair per event is kept; the K_S0 mass is restricted to
    [487, 511] MeV/c^2")
  .note(:lambda_vertex_fit, "the Lambda -> p pi- candidate is required to have a vertex-fit
    chi2 < 100 and a secondary-vertex constraint pointing back to the IP with a decay length
    larger than twice its resolution; only one track pair per event is kept (the Lambda mass
    window is expressed as the invariant-mass window inside the kinematic fit)")

tag_alg.apply
root_files = tag_alg.execute_on([data_4600, incMC_4600, exMC_signal])