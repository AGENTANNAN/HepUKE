# Core DSL classes and dependencies will be loaded automatically at execution

### Datasets ###
data_3773  = DatasetManager.real_data.find("712_3773")     # psi(3770) real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

### Signal decay cards (EvtGen format) ###
# Channel 1: D+ -> K+ pi0 pi0 (doubly Cabibbo-suppressed), tagged on the D- side
decay_card_KPi0Pi0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K+ pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Channel 2: D+ -> K+ pi0 eta (doubly Cabibbo-suppressed), tagged on the D- side
decay_card_KPi0Eta = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K+ pi0 eta PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Channel 1: D- tag + D+ -> K+ pi0 pi0 ###
alg_KPi0Pi0 = TagAnalysis.new("DpToKPi0Pi0")
alg_KPi0Pi0.set_header(["DpToKPi0Pi0Alg/DpToKPi0Pi0.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })
           .with_decay_card(decay_card_KPi0Pi0)

# D- tagged through three hadronic modes; charm = -1 pins the tagged side to D-
alg_KPi0Pi0.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0   # D- -> K+ pi- pi-, K_S0 pi-, K+ pi- pi- pi0
  t.charm -1
end

# Signal side: one K+ and four photons (two pi0)
alg_KPi0Pi0.signal_side do |s|
  s.charged(kp: 1)
  s.photons 4
end

# 4C kinematic fit: four-momentum conservation, chi2 < 200
alg_KPi0Pi0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# D-tag reconstruction: beam energy from the conditions DB, local mode
alg_KPi0Pi0.dtag_reconstruction do |d|
  d.beam_energy :db
  d.local true
end

alg_KPi0Pi0.apply
alg_KPi0Pi0.execute_on([data_3773, incMC_3773])

### Channel 2: D- tag + D+ -> K+ pi0 eta ###
alg_KPi0Eta = TagAnalysis.new("DpToKPi0Eta")
alg_KPi0Eta.set_header(["DpToKPi0EtaAlg/DpToKPi0Eta.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
          .with_decay_card(decay_card_KPi0Eta)

# Same D- tag modes with charm = -1
alg_KPi0Eta.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0   # D- -> K+ pi- pi-, K_S0 pi-, K+ pi- pi- pi0
  t.charm -1
end

# Signal side: one K+ and two photons
alg_KPi0Eta.signal_side do |s|
  s.charged(kp: 1)
  s.photons 2
end

# 4C kinematic fit: four-momentum conservation, chi2 < 200
alg_KPi0Eta.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# D-tag reconstruction: beam energy from the conditions DB, local mode
alg_KPi0Eta.dtag_reconstruction do |d|
  d.beam_energy :db
  d.local true
end

alg_KPi0Eta.apply
alg_KPi0Eta.execute_on([data_3773, incMC_3773])