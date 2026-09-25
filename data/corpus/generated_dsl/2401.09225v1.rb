# =============================================================================
# e+e- -> anti-Lambda_c- Lambda_c(2625)+ , Lambda_c(2625)+ -> Lambda_c+ pi+ pi-
# Single-tag (Lambda_c+) partial-reconstruction analysis at
# sqrt(s) = 4.918 GeV and 4.950 GeV (analysed separately and combined)
# =============================================================================

### Dataset preparation ###
data_4918  = DatasetManager.real_data.find("707_4914")     # 4.918 GeV real data (219.78 pb^-1)
incMC_4918 = DatasetManager.inclusive_mc.find("707_4914")  # inclusive MC at 4.918 GeV
data_4950  = DatasetManager.real_data.find("707_4946")     # 4.950 GeV real data (148.70 pb^-1)
incMC_4950 = DatasetManager.inclusive_mc.find("707_4946")  # inclusive MC at 4.950 GeV

# Signal decay card (EvtGen): e+e- -> anti-Lambda_c- Lambda_c(2625)+
# Lambda_c(2625)+ -> Lambda_c+ pi+ pi- in phase space,
# Lambda_c+ -> p+ K- pi+ (representative tag mode)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 anti-Lambda_c- Lambda_c(2625)+ PHSP;
    Enddecay

    Decay Lambda_c(2625)+
    1.0000 Lambda_c+ pi+ pi- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 p+ K- pi+ PHSP;
    Enddecay

    End
DECAYCARD

# Charge-conjugate decay card: e+e- -> Lambda_c+ anti-Lambda_c(2625)-
decay_card_conjugate = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c(2625)- PHSP;
    Enddecay

    Decay anti-Lambda_c(2625)-
    1.0000 anti-Lambda_c- pi- pi+ PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for the signal, and another 500k for its charge conjugate
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4918_antiLc_Lc2625_signal"
  config.related_dataset = data_4918
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

exMC_conjugate = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4918_antiLc_Lc2625_conjugate"
  config.related_dataset = data_4918
  config.events          = 500_000
  config.decay_card      = decay_card_conjugate
  config.cross_section   = :default
end

### Tag analysis (BOSS) ###
alg_name = "Lc2625Tag"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.918]})   # c.m. energy constant = 4.918 GeV
   .with_decay_card(decay_card_signal)
   .note(:delta_m_peak, "Delta M = M(Lambda_c+ pi+ pi-) - M(Lambda_c+) peak expected near 0.341 GeV; extracted offline in the ROOT analysis")
   .note(:mrecoil_window, "M_recoil window around the anti-Lambda_c- nominal mass (2.28646 GeV/c2) applied offline for the yield extraction")
   .note(:mbc_window, "beam-constrained-mass window of 5 MeV/c2 for the p K- pi+ tag mode applied offline for the yield extraction")

# Tag side: single tag on Lambda_c+ through three hadronic tag modes
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP, :LambdacPtoKsP, :LambdacPtoLambdaPi   # p K- pi+, p K_S0, Lambda pi+
  t.charm 1                                                     # tag the Lambda_c+ side
end

# Signal side: pi+ pi- recoiling against the tag; anti-Lambda_c- is missing
alg.signal_side do |s|
  s.photons 0                     # zero photons on the signal side
  s.charged(pip: 1, pim: 1)       # exactly one pi+ and one pi-
  s.require_charge 0              # net signal-side charge zero
  s.missing :Lambdac              # missing anti-Lambda_c-, nominal mass 2.28646 GeV/c2
end

# 4C kinematic fit: constrain the total four-momentum to the c.m. energy
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# Render the tag spec (no Selection argument) and run on the datasets
alg.apply
alg.execute_on([data_4918, incMC_4918, data_4950, incMC_4950, exMC_signal, exMC_conjugate])