# Companion paper to 2607.23945: same D -> pi lepton nu channels, decay dynamics measurement.
# BOSS-level selection is identical.
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

decay_card_D0Pienu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0                    VSS_BMIX dm;
  Enddecay

  Decay D0
  1.0000 pi- e+ nu_e                   PHOTOS  ISGW2;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi-                        PHSP;
  Enddecay

  End
DECAYCARD

decay_card_D0Pimunu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0                    VSS_BMIX dm;
  Enddecay

  Decay D0
  1.0000 pi- mu+ nu_mu                 PHOTOS  ISGW2;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi-                        PHSP;
  Enddecay

  End
DECAYCARD

decay_card_DpPi0enu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-                         PHSP;
  Enddecay

  Decay D+
  1.0000 pi0 e+ nu_e                   PHOTOS  ISGW2;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi-                    PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                   PHSP;
  Enddecay

  End
DECAYCARD

decay_card_DpPi0munu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-                         PHSP;
  Enddecay

  Decay D+
  1.0000 pi0 mu+ nu_mu                 PHOTOS  ISGW2;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi-                    PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                   PHSP;
  Enddecay

  End
DECAYCARD

exMC_D0Pienu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "D0_pi_e_nu_dyn"
  c.related_dataset = data_3773
  c.events          = 500000
  c.decay_card      = decay_card_D0Pienu
  c.cross_section   = :default
end

exMC_D0Pimunu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "D0_pi_mu_nu_dyn"
  c.related_dataset = data_3773
  c.events          = 500000
  c.decay_card      = decay_card_D0Pimunu
  c.cross_section   = :default
end

exMC_DpPi0enu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Dp_pi0_e_nu_dyn"
  c.related_dataset = data_3773
  c.events          = 500000
  c.decay_card      = decay_card_DpPi0enu
  c.cross_section   = :default
end

exMC_DpPi0munu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Dp_pi0_mu_nu_dyn"
  c.related_dataset = data_3773
  c.events          = 500000
  c.decay_card      = decay_card_DpPi0munu
  c.cross_section   = :default
end

# Channel 1: D0 -> pi- e+ nu_e
alg1 = TagAnalysis.new("D0DynPiENu")
alg1.set_header(["D0DynPiENuAlg/D0DynPiENu.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })
    .with_decay_card(decay_card_D0Pienu)

alg1.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPi0Pi0, :D0toKPiPiPiPi0
  t.charm -1
end

alg1.signal_side do |s|
  s.charged(pim: 1, ep: 1)
  s.require_charge 0
  s.missing :nu_e
end

alg1.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg1.apply
alg1.execute_on([data_3773, incMC_3773, exMC_D0Pienu])

# Channel 2: D0 -> pi- mu+ nu_mu
alg2 = TagAnalysis.new("D0DynPiMuNu")
alg2.set_header(["D0DynPiMuNuAlg/D0DynPiMuNu.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })
    .with_decay_card(decay_card_D0Pimunu)

alg2.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPi0Pi0, :D0toKPiPiPiPi0
  t.charm -1
end

alg2.signal_side do |s|
  s.charged(pim: 1, mup: 1)
  s.require_charge 0
  s.missing :nu_mu
end

alg2.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg2.apply
alg2.execute_on([data_3773, incMC_3773, exMC_D0Pimunu])

# Channel 3: D+ -> pi0 e+ nu_e
alg3 = TagAnalysis.new("DpDynPi0ENu")
alg3.set_header(["DpDynPi0ENuAlg/DpDynPi0ENu.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })
    .with_decay_card(decay_card_DpPi0enu)

alg3.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

alg3.signal_side do |s|
  s.photons 2
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

alg3.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg3.apply
alg3.execute_on([data_3773, incMC_3773, exMC_DpPi0enu])

# Channel 4: D+ -> pi0 mu+ nu_mu
alg4 = TagAnalysis.new("DpDynPi0MuNu")
alg4.set_header(["DpDynPi0MuNuAlg/DpDynPi0MuNu.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })
    .with_decay_card(decay_card_DpPi0munu)

alg4.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

alg4.signal_side do |s|
  s.photons 2
  s.charged(mup: 1)
  s.require_charge 1
  s.missing :nu_mu
end

alg4.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg4.apply
alg4.execute_on([data_3773, incMC_3773, exMC_DpPi0munu])
