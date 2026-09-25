# =====================================================================
# Double-tagged hadronic D0 / D+ analysis at psi(3770), 3.773 GeV
# Tag-based (DTagTool) analysis: tag side from standard hadronic ST modes,
# signal side rebuilt from the remaining tracks.
# =====================================================================

### Datasets ###
data_3773  = DatasetManager.real_data.find("712_3773")     # psi(3770) real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # matching inclusive MC

### Decay cards (one exclusive-MC signal process per mode) ###
# Mode 1: D0 -> K- 3pi+ 2pi-,  tag anti-D0 -> K+ pi-
decay_card_mode1 = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.000 K- pi+ pi+ pi+ pi- pi- PHSP;
  Enddecay

  Decay anti-D0
  1.000 K+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode 2: D0 -> K- 2pi+ pi- 2pi0 (pi0 -> gamma gamma),  tag anti-D0 -> K+ pi-
decay_card_mode2 = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.000 K- pi+ pi+ pi- pi0 pi0 PHSP;
  Enddecay

  Decay anti-D0
  1.000 K+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode 3: D+ -> K- 3pi+ pi- pi0 (pi0 -> gamma gamma),  tag D- -> K+ pi- pi-
decay_card_mode3 = <<~DECAYCARD
  Decay psi(3770)
  1.000 D+ D- PHSP;
  Enddecay

  Decay D+
  1.000 K- pi+ pi+ pi+ pi- pi0 PHSP;
  Enddecay

  Decay D-
  1.000 K+ pi- pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC (100k events per signal mode) ###
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_D0_K3pi2pi_tag_antiD0Kpi"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_mode1
  config.cross_section   = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_D0_K2pi1pi2pi0_tag_antiD0Kpi"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_mode2
  config.cross_section   = :default
end

exMC_mode3 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_Dplus_K3pi1pi1pi0_tag_DminusKpipi"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_mode3
  config.cross_section   = :default
end

### Mode 1: D0 -> K- 3pi+ 2pi-  (tagged by anti-D0) ###
alg_mode1 = TagAnalysis.new("DTagD0K3Pi2Pi")
alg_mode1.set_header(["DTagD0K3Pi2PiAlg/DTagD0K3Pi2Pi.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_mode1)

# Tag side: anti-D0 from the standard hadronic ST modes (charge conjugates handled internally)
alg_mode1.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1                       # pin the tagged D to the anti-D0
end

# Signal side: remaining tracks with the required charged multiplicity (0 photons)
alg_mode1.signal_side do |s|
  s.charged(km: 1, pip: 3, pim: 2)
end

# Four-momentum-constraint fit (4C), loose chi2
alg_mode1.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mode1.apply
root_files_mode1 = alg_mode1.execute_on([data_3773, incMC_3773, exMC_mode1])

### Mode 2: D0 -> K- 2pi+ pi- 2pi0 (pi0 -> gamma gamma)  (tagged by anti-D0) ###
alg_mode2 = TagAnalysis.new("DTagD0K2Pi1Pi2Pi0")
alg_mode2.set_header(["DTagD0K2Pi1Pi2Pi0Alg/DTagD0K2Pi1Pi2Pi0.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_mode2)

alg_mode2.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: 1 K-, 2 pi+, 1 pi- and the 4 photons from the two pi0
alg_mode2.signal_side do |s|
  s.photons 4
  s.charged(km: 1, pip: 2, pim: 1)
end

alg_mode2.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mode2.apply
root_files_mode2 = alg_mode2.execute_on([data_3773, incMC_3773, exMC_mode2])

### Mode 3: D+ -> K- 3pi+ pi- pi0 (pi0 -> gamma gamma)  (tagged by D-) ###
alg_mode3 = TagAnalysis.new("DTagDpK3Pi1Pi1Pi0")
alg_mode3.set_header(["DTagDpK3Pi1Pi1Pi0Alg/DTagDpK3Pi1Pi1Pi0.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_mode3)

# Tag side: D- from the standard hadronic ST modes of the charged D
alg_mode3.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1                       # pin the tagged D to the D-
end

# Signal side: 1 K-, 3 pi+, 1 pi- and the 2 photons from the pi0
alg_mode3.signal_side do |s|
  s.photons 2
  s.charged(km: 1, pip: 3, pim: 1)
end

alg_mode3.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mode3.apply
root_files_mode3 = alg_mode3.execute_on([data_3773, incMC_3773, exMC_mode3])