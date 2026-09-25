### Dataset description ###
# 12 real-data energy points spanning 4.6119 - 4.9509 GeV  (BOSS 706 & 707)
data_scan = [
  DatasetManager.real_data.find("706_4610"),  # 4.6119 GeV
  DatasetManager.real_data.find("706_4620"),  # 4.6280 GeV
  DatasetManager.real_data.find("706_4640"),  # 4.6409 GeV
  DatasetManager.real_data.find("706_4660"),  # 4.6612 GeV
  DatasetManager.real_data.find("706_4680"),  # 4.6819 GeV
  DatasetManager.real_data.find("706_4700"),  # 4.6988 GeV
  DatasetManager.real_data.find("707_4740"),  # 4.7397 GeV
  DatasetManager.real_data.find("707_4750"),  # 4.7501 GeV
  DatasetManager.real_data.find("707_4780"),  # 4.7805 GeV
  DatasetManager.real_data.find("707_4840"),  # 4.8431 GeV
  DatasetManager.real_data.find("707_4914"),  # 4.9180 GeV
  DatasetManager.real_data.find("707_4946")   # 4.9509 GeV
]

# 4.5995 GeV sample kept for the double-tag study (ROOT-level 2D simultaneous fit)
data_dt = DatasetManager.real_data.find("703_4600")

# Inclusive MC at 4.682 GeV
incMC = DatasetManager.inclusive_mc.find("706_4680")

### Decay card: e+e- -> Lambda_c+ Lambda_c- , tag mode Lambda_c+ -> p K- pi+ ###
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 p+ K- pi+ PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC: 200k events per energy point for Lambda_c+ -> p K- pi+ (PHSP) ###
exMCs = DatasetManager.create_exclusive_mc_for(data_scan + [data_dt]) do |config|
  config.sample_name   = "exmc_lambdac_pKpi"   # auto-suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) — tag-based single-tag analysis ###
alg_name = "LambdacTagPKPi"
tag_alg = TagAnalysis.new(alg_name)
tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({ "ECMS" => [:double, 4.682] })   # fit uses MeasuredEcmsSvc per run
       .with_decay_card(decay_card_signal)

# Tag side: a single Lambda_c+ -> p K- pi+ candidate (single-tag method; charm = +1 pins Lambda_c+)
tag_alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP
  t.charm 1
end

# Signal side: the other Lambda_c- is left unreconstructed (one massive missing
# particle); no photon is allowed on the signal side.  M_BC is stored
# unconditionally (store-not-cut) and windowed later in ROOT.
tag_alg.signal_side do |s|
  s.photons 0
  s.missing :Lambdac
end

# 4C kinematic fit with chi2 < 200
tag_alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

tag_alg.apply
# Execute on all data points, the 4.5995 GeV DT sample, inclusive MC and the signal MC set
tag_alg.execute_on(data_scan + [data_dt, incMC] + exMCs)