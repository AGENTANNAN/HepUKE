### Paper 1912.12411v3: D -> eta eta pi+, eta pi+ pi0, eta pi+ pi- at psi(3770)
### TAG-BASED: ST/DT method using TagAnalysis
### Sample: psi(3770) at 3.773 GeV (712_3773)
### Three signal decays measured via double-tag (DT) method

data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

### ====================================================================
### Signal 1: D0 -> eta pi+ pi-  (tagged by anti-D0)
### Tag side: D0bar with 3 hadronic modes
### ====================================================================

decay_card_sig1 = <<~DECAYCARD
    Decay anti-D0
    1  eta  pi+  pi-    PHSP;
    Enddecay

    Decay eta
    1  gamma  gamma    PHSP;
    Enddecay

    End
DECAYCARD

exMC_sig1 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_D0toEtapipim"
  c.related_dataset = data_3773
  c.events = 100000
  c.decay_card = decay_card_sig1
  c.cross_section = :default
end

alg_sig1 = TagAnalysis.new("D0toEtapipim")
alg_sig1.set_header(["D0toEtapipimAlg/D0toEtapipim.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_sig1)

### Tag side: anti-D0 -> K+pi-, K+pi-pi0, K+pi-pi-pi+
alg_sig1.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1    # anti-D0 tag
end

### Signal side: eta -> gamma gamma, pi+, pi-
alg_sig1.signal_side do |s|
  s.photons 2       # eta -> gamma gamma
  s.charged(pip: 1, pim: 1)
end

### Kinematic fit: 4C + eta mass constraint
alg_sig1.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_sig1.apply
alg_sig1.execute_on([data_3773, incMC_3773, exMC_sig1])

### ====================================================================
### Signal 2: D+ -> eta eta pi+  (tagged by D-)
### Tag side: D- with 6 hadronic modes
### ====================================================================

decay_card_sig2 = <<~DECAYCARD
    Decay D+
    1  eta  eta  pi+    PHSP;
    Enddecay

    Decay eta
    1  gamma  gamma    PHSP;
    Enddecay

    End
DECAYCARD

exMC_sig2 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_DptoEtaEtapip"
  c.related_dataset = data_3773
  c.events = 100000
  c.decay_card = decay_card_sig2
  c.cross_section = :default
end

alg_sig2 = TagAnalysis.new("DptoEtaEtapip")
alg_sig2.set_header(["DptoEtaEtapipAlg/DptoEtaEtapip.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_sig2)

### Tag side: D- with 6 hadronic modes
alg_sig2.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1    # D- tag
end

### Signal side: 2 eta -> 4 gamma, pi+
alg_sig2.signal_side do |s|
  s.photons 4       # 2 eta -> 4 gamma
  s.charged(pip: 1)
end

### Kinematic fit: 4C + two eta mass constraints
alg_sig2.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)   # first eta
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)   # second eta
  f.chi2_cut 200
end

alg_sig2.apply
alg_sig2.execute_on([data_3773, incMC_3773, exMC_sig2])

### ====================================================================
### Signal 3: D+ -> eta pi+ pi0  (tagged by D-)
### Tag side: D- with 6 hadronic modes (same as Signal 2)
### ====================================================================

decay_card_sig3 = <<~DECAYCARD
    Decay D+
    1  eta  pi+  pi0    PHSP;
    Enddecay

    Decay eta
    1  gamma  gamma    PHSP;
    Enddecay

    Decay pi0
    1  gamma  gamma    PHSP;
    Enddecay

    End
DECAYCARD

exMC_sig3 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_DptoEtapipPi0"
  c.related_dataset = data_3773
  c.events = 100000
  c.decay_card = decay_card_sig3
  c.cross_section = :default
end

alg_sig3 = TagAnalysis.new("DptoEtapipPi0")
alg_sig3.set_header(["DptoEtapipPi0Alg/DptoEtapipPi0.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_sig3)

### Tag side: D- with 6 hadronic modes
alg_sig3.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1    # D- tag
end

### Signal side: eta -> 2 gamma, pi0 -> 2 gamma, pi+
alg_sig3.signal_side do |s|
  s.photons 4       # 2 from eta + 2 from pi0
  s.charged(pip: 1)
end

### Kinematic fit: 4C + eta mass constraint + pi0 mass constraint
alg_sig3.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_sig3.apply
alg_sig3.execute_on([data_3773, incMC_3773, exMC_sig3])