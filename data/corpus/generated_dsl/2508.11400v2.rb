### Dataset preparation ###
# Single-tag analysis at 13 centre-of-mass energies from 4.600 to 4.951 GeV (6.4 fb^-1 total).
# Sample-name convention: [BOSS version]_[CMS energy in MeV].
energy_sample_names = %w[
  703_4600 706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840 707_4914 707_4946
]
lc_data  = energy_sample_names.map { |name| DatasetManager.real_data.find(name) }     # real data at each point
lc_incMC = energy_sample_names.map { |name| DatasetManager.inclusive_mc.find(name) } # inclusive MC at each point

# Decay card: e+e- -> Lambda_c+ anti-Lambda_c-.
# Lambda_c+ decays with equal 20% probability to each of the five signal modes;
# anti-Lambda_c- is left to decay inclusively (not fixed here).
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  0.200 p+ K_S0      PHSP;
  0.200 Lambda0 pi+  PHSP;
  0.200 Sigma0 pi+   PHSP;
  0.200 Sigma+ pi0   PHSP;
  0.200 p+ K- pi+    PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- PHSP;
  Enddecay

  Decay Sigma0
  1.000 Lambda0 gamma PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 200k-event exclusive signal MC at every energy point
exMC_signal = DatasetManager.create_exclusive_mc_for(lc_data) do |config|
  config.sample_name   = "LcLcbar_signal_mc"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (tag-based analysis) ###
alg_name = "LcTagPolarization"
tag_analysis = TagAnalysis.new(alg_name)
tag_analysis.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.775]}) # representative point; per-run beam energy read from DB
            .with_decay_card(decay_card_signal)

# Tag side: single tag of Lambda_c+ with only the three modes present in the tagger vocabulary.
# Sigma0 pi+ and Sigma+ pi0 have no DTagAlg counterpart and are therefore dropped.
tag_analysis.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP, :LambdacPtoKsP, :LambdacPtoLambdaPiP
  t.charm 1 # pin the tagged Lambda_c+
end

# Signal side: the anti-Lambda_c- is not reconstructed (decays inclusively); require zero photons.
tag_analysis.signal_side do |s|
  s.photons 0
  s.missing :"Lambda_c-" # anti-Lambda_c- treated as one missing particle
end

# Kinematic fit: four-momentum conservation plus the tag Lambda_c+ mass constraint.
tag_analysis.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
  f.store_fitted_momenta # fitted four-momenta fed to the ROOT angular analysis
end

tag_analysis.note(:inclusive_antilambdac_decay,
  "anti-Lambda_c- decays inclusively and is not fixed in the decay card; it is treated as a single " \
  "missing particle (mass = m(Lambda_c)) in the 4C kinematic fit")

tag_analysis.apply
tag_analysis.execute_on(lc_data + lc_incMC + exMC_signal)