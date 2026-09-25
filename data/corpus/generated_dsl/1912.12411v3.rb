### ======================================================================
###  BOSS DSL — psi(3770) (3.773 GeV) double-tag analyses
###    (1) D0 -> eta pi+ pi-   (eta -> gamma gamma)          , tag anti-D0
###    (2) D+ -> eta eta pi+   (eta -> gamma gamma x2)       , tag D-
###    (3) D+ -> eta pi+ pi0   (eta -> gamma gamma, pi0 -> gamma gamma), tag D-
###  The tag side comes from the pre-stored DTag collection (DTagTool);
###  the signal side is built from the tracks/showers the tag did not use.
### ======================================================================

### ---------------------- Dataset preparation ----------------------
data_3773  = DatasetManager.real_data.find("712_3773")      # psi(3770) real data @ 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC

# ---- Decay cards (EvtGen), top mother = psi(3770) ----
# Mode 1 : psi(3770) -> D0 anti-D0 ; D0 -> eta pi+ pi- ; eta -> gamma gamma
#          (anti-D0 tag decays through the Kpi / Kpipi0 / Kpipipi channels)
decay_card_D0_etapipi = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0 eta pi+ pi- PHSP;
  Enddecay
  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay
  Decay anti-D0
  0.334 K+ pi- PHSP;
  0.333 K+ pi- pi0 PHSP;
  0.333 K+ pi- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Mode 2 : psi(3770) -> D+ D- ; D+ -> eta eta pi+ ; eta -> gamma gamma
#          (D- tag decays through the six hadronic channels)
decay_card_Dp_etaetapi = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D- PHSP;
  Enddecay
  Decay D+
  1.0 eta eta pi+ PHSP;
  Enddecay
  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay
  Decay D-
  0.167 K+ pi- pi- PHSP;
  0.167 K_S0 pi- PHSP;
  0.167 K+ pi- pi- pi0 PHSP;
  0.167 K_S0 pi- pi0 PHSP;
  0.167 K_S0 pi- pi+ pi- PHSP;
  0.167 K+ K- pi- PHSP;
  Enddecay
  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Mode 3 : psi(3770) -> D+ D- ; D+ -> eta pi+ pi0 ; eta -> gamma gamma ; pi0 -> gamma gamma
#          (same six D- tag channels)
decay_card_Dp_etapipi0 = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D- PHSP;
  Enddecay
  Decay D+
  1.0 eta pi+ pi0 PHSP;
  Enddecay
  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay
  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay
  Decay D-
  0.167 K+ pi- pi- PHSP;
  0.167 K_S0 pi- PHSP;
  0.167 K+ pi- pi- pi0 PHSP;
  0.167 K_S0 pi- pi0 PHSP;
  0.167 K_S0 pi- pi+ pi- PHSP;
  0.167 K+ K- pi- PHSP;
  Enddecay
  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# 100k-event exclusive MC for each of the three signal modes
exMC_D0_etapipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_D0_etapipi_dtag"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_D0_etapipi
  config.cross_section   = :default
end

exMC_Dp_etaetapi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_Dp_etaetapi_dtag"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_Dp_etaetapi
  config.cross_section   = :default
end

exMC_Dp_etapipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_Dp_etapipi0_dtag"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_Dp_etapipi0
  config.cross_section   = :default
end

### --------------------------------------------------------------------
### Analysis 1 : D0 -> eta pi+ pi-   (eta -> gamma gamma), tag = anti-D0
### --------------------------------------------------------------------
alg_mode1 = TagAnalysis.new("D0ToEtaPiPiDT")
alg_mode1.set_header(["D0ToEtaPiPiDTAlg/D0ToEtaPiPiDT.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_D0_etapipi)

# Tag side : anti-D0 reconstructed in Kpi, Kpipi0, Kpipipi
alg_mode1.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1                              # pin the anti-D0 tag
end

# Signal side : the tag's unused tracks/showers -> 2 photons + pi+ + pi-
alg_mode1.signal_side do |s|
  s.photons 2
  s.charged(pip: 1, pim: 1)
end

# 4C kinematic fit (chi2 < 200) with the eta mass constraint
alg_mode1.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_mode1.apply
alg_mode1.execute_on([data_3773, incMC_3773, exMC_D0_etapipi])

### --------------------------------------------------------------------
### Analysis 2 : D+ -> eta eta pi+   (eta -> gamma gamma x2), tag = D-
### --------------------------------------------------------------------
alg_mode2 = TagAnalysis.new("DpToEtaEtaPiDT")
alg_mode2.set_header(["DpToEtaEtaPiDTAlg/DpToEtaEtaPiDT.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_Dp_etaetapi)

# Tag side : D- reconstructed in the six hadronic modes
alg_mode2.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1                              # pin the D- tag
end

# Signal side : the tag's unused tracks/showers -> 4 photons + pi+
alg_mode2.signal_side do |s|
  s.photons 4
  s.charged(pip: 1)
end

# 4C kinematic fit (chi2 < 200) with two eta mass constraints
alg_mode2.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_mode2.apply
alg_mode2.execute_on([data_3773, incMC_3773, exMC_Dp_etaetapi])

### --------------------------------------------------------------------
### Analysis 3 : D+ -> eta pi+ pi0   (eta,pi0 -> gamma gamma), tag = D-
### --------------------------------------------------------------------
alg_mode3 = TagAnalysis.new("DpToEtaPiPi0DT")
alg_mode3.set_header(["DpToEtaPiPi0DTAlg/DpToEtaPiPi0DT.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_Dp_etapipi0)

# Tag side : same six D- hadronic modes
alg_mode3.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1                              # pin the D- tag
end

# Signal side : the tag's unused tracks/showers -> 4 photons + pi+
alg_mode3.signal_side do |s|
  s.photons 4
  s.charged(pip: 1)
end

# 4C kinematic fit (chi2 < 200) with eta and pi0 mass constraints
alg_mode3.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_mode3.apply
alg_mode3.execute_on([data_3773, incMC_3773, exMC_Dp_etapipi0])