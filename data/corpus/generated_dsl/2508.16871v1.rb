# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# e+e- -> Sigma_c Sigma_cbar (4.918, 4.951 GeV)
# e+e- -> Lambda_c+ Sigma_cbar- (4.750, 4.781, 4.843, 4.918, 4.951 GeV)
# BOSS 707 scan samples (sample name = [BOSS]_[CMS energy in MeV])
data_4750 = DatasetManager.real_data.find("707_4750")   # 4.750 GeV
data_4780 = DatasetManager.real_data.find("707_4780")   # 4.781 GeV
data_4840 = DatasetManager.real_data.find("707_4840")   # 4.843 GeV
data_4914 = DatasetManager.real_data.find("707_4914")   # 4.918 GeV
data_4946 = DatasetManager.real_data.find("707_4946")   # 4.951 GeV

incMC_4750 = DatasetManager.inclusive_mc.find("707_4750")
incMC_4780 = DatasetManager.inclusive_mc.find("707_4780")
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

data_points  = [data_4750, data_4780, data_4840, data_4914, data_4946]
incMC_points = [incMC_4750, incMC_4780, incMC_4840, incMC_4914, incMC_4946]

# Decay card (EvtGen): psi(4520) -> Lambda_c+ anti-Sigma_c- phase space,
# with the four Lambda_c+ modes at 25% each plus the K_S0 / Lambda / Sigma0 cascades
decay_card_signal = <<~DECAYCARD
    Decay psi(4520)
    1.0000 Lambda_c+ anti-Sigma_c- PHSP;
    Enddecay

    Decay Lambda_c+
    0.25 p+ K_S0        PHSP;
    0.25 p+ K- pi+      PHSP;
    0.25 Lambda0 pi+    PHSP;
    0.25 Sigma0 pi+     PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi-       PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi-        HypWK;
    Enddecay

    Decay Sigma0
    1.000 Lambda0 gamma PHSP;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive signal MC, one sample per scan energy point
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_Lambdac_Sigmacbar"
  config.events        = 200000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) — single-tag (ST) Lambda_c+ ###
alg_name = "LambdacTag"
my_Algorithm = TagAnalysis.new(alg_name)          # TagAnalysis replaces Algorithm + Selection
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.918]})
            .with_decay_card(decay_card_signal)
            .note(:tag_mode_unavailable, "Lambda_c+ -> Sigma0 pi+ tag mode dropped: no matching
              DTagAlg decay-mode symbol in the vocabulary")
            .note(:tag_sample_scope, "ST, ST+pi and ST+pi pi tag samples are all foreseen;
              only the single-tag sample is encoded here")

# Tag side: the tagging Lambda_c+ is taken from the pre-stored DTag candidates (single tag)
my_Algorithm.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP, :LambdacPtoKsP, :LambdacPtoLamPi   # p K- pi+ , p K_S0 , Lambda pi+
  t.charm 1                                                    # pin the tagged side to Lambda_c+
end

# Signal side: no photons required at BOSS level
my_Algorithm.signal_side do |s|
  s.photons 0
end

# 4-momentum kinematic fit with the tag Lambda_c+ invariant mass constrained to its nominal value
my_Algorithm.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

my_Algorithm.apply   # takes no Selection argument

# Execute on real data, inclusive MC and the exclusive signal MC
root_files = my_Algorithm.execute_on(data_points + incMC_points + exMC_signal)