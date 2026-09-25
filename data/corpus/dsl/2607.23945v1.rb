# Dataset preparation: BESIII 20.3 fb^-1 psi(3770) data collected at sqrt(s)=3.773 GeV
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Signal decay cards — one per channel:
# D0 -> pi- e+ nu_e ; D0 -> pi- mu+ nu_mu ; D+ -> pi0 e+ nu_e ; D+ -> pi0 mu+ nu_mu
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
  c.sample_name     = "D0_pi_e_nu"
  c.related_dataset = data_3773
  c.events          = 500000
  c.decay_card      = decay_card_D0Pienu
  c.cross_section   = :default
end

exMC_D0Pimunu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "D0_pi_mu_nu"
  c.related_dataset = data_3773
  c.events          = 500000
  c.decay_card      = decay_card_D0Pimunu
  c.cross_section   = :default
end

exMC_DpPi0enu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Dp_pi0_e_nu"
  c.related_dataset = data_3773
  c.events          = 500000
  c.decay_card      = decay_card_DpPi0enu
  c.cross_section   = :default
end

exMC_DpPi0munu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Dp_pi0_mu_nu"
  c.related_dataset = data_3773
  c.events          = 500000
  c.decay_card      = decay_card_DpPi0munu
  c.cross_section   = :default
end

# ---------------------------------------------------------------------------
# Channel 1: D0 -> pi- e+ nu_e   (tag: anti-D0 -> multi hadronic modes; missing nu_e)
# ---------------------------------------------------------------------------
alg_D0e = TagAnalysis.new("D0TagPiENu")
alg_D0e.set_header(["D0TagPiENuAlg/D0TagPiENu.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .with_decay_card(decay_card_D0Pienu)

alg_D0e.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPi0Pi0, :D0toKPiPiPiPi0
  t.charm -1
end

alg_D0e.signal_side do |s|
  s.charged(pim: 1, ep: 1)
  s.require_charge 0
  s.missing :nu_e                 # massless
end

alg_D0e.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_D0e.apply
alg_D0e.execute_on([data_3773, incMC_3773, exMC_D0Pienu])

# ---------------------------------------------------------------------------
# Channel 2: D0 -> pi- mu+ nu_mu (tag: anti-D0 -> multi hadronic modes; missing nu_mu)
# ---------------------------------------------------------------------------
alg_D0mu = TagAnalysis.new("D0TagPiMuNu")
alg_D0mu.set_header(["D0TagPiMuNuAlg/D0TagPiMuNu.h"])
        .set_constant({ "ECMS" => [:double, 3.773] })
        .with_decay_card(decay_card_D0Pimunu)

alg_D0mu.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPi0Pi0, :D0toKPiPiPiPi0
  t.charm -1
end

alg_D0mu.signal_side do |s|
  s.charged(pim: 1, mup: 1)
  s.require_charge 0
  s.missing :nu_mu
end

alg_D0mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_D0mu.apply
alg_D0mu.execute_on([data_3773, incMC_3773, exMC_D0Pimunu])

# ---------------------------------------------------------------------------
# Channel 3: D+ -> pi0 e+ nu_e (tag: D- -> multi hadronic modes; missing nu_e; pi0 -> gamma gamma)
# ---------------------------------------------------------------------------
alg_Dpe = TagAnalysis.new("DpTagPi0ENu")
alg_Dpe.set_header(["DpTagPi0ENuAlg/DpTagPi0ENu.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .with_decay_card(decay_card_DpPi0enu)

alg_Dpe.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

alg_Dpe.signal_side do |s|
  s.photons 2
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

alg_Dpe.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_Dpe.apply
alg_Dpe.execute_on([data_3773, incMC_3773, exMC_DpPi0enu])

# ---------------------------------------------------------------------------
# Channel 4: D+ -> pi0 mu+ nu_mu (tag: D- -> hadronic modes; missing nu_mu; pi0 -> gamma gamma)
# ---------------------------------------------------------------------------
alg_Dpmu = TagAnalysis.new("DpTagPi0MuNu")
alg_Dpmu.set_header(["DpTagPi0MuNuAlg/DpTagPi0MuNu.h"])
        .set_constant({ "ECMS" => [:double, 3.773] })
        .with_decay_card(decay_card_DpPi0munu)

alg_Dpmu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

alg_Dpmu.signal_side do |s|
  s.photons 2
  s.charged(mup: 1)
  s.require_charge 1
  s.missing :nu_mu
end

alg_Dpmu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_Dpmu.apply
alg_Dpmu.execute_on([data_3773, incMC_3773, exMC_DpPi0munu])
