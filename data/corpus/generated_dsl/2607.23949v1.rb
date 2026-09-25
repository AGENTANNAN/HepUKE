### Dataset description ###
# ψ(3770) at 3.773 GeV — single-tag semileptonic D decays
data_3773 = DatasetManager.real_data.find("712_3773")       # 3.773 GeV real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # matching inclusive MC

### Decay cards (EvtGen) — one per signal channel ###
# The signal D undergoes the semileptonic decay; the opposite (tag) D decays via
# a representative tag mode (the full tag-mode menu is declared in tag_side).
decay_card_d0_pi_enu = <<~DECAYCARD
    Decay psi(3770)
    1.0 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0 pi- e+ nu_e PHSP;
    Enddecay

    Decay anti-D0
    1.0 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_d0_pi_munu = <<~DECAYCARD
    Decay psi(3770)
    1.0 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0 pi- mu+ nu_mu PHSP;
    Enddecay

    Decay anti-D0
    1.0 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_dp_pi0_enu = <<~DECAYCARD
    Decay psi(3770)
    1.0 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0 pi0 e+ nu_e PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    Decay D-
    1.0 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_dp_pi0_munu = <<~DECAYCARD
    Decay psi(3770)
    1.0 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0 pi0 mu+ nu_mu PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    Decay D-
    1.0 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC — 500k events for each of the four signal modes ###
exMC_d0_pi_enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_D0ToPiEneu"
  config.related_dataset = data_3773
  config.events = 500_000
  config.decay_card = decay_card_d0_pi_enu
  config.cross_section = :default
end

exMC_d0_pi_munu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_D0ToPiMunu"
  config.related_dataset = data_3773
  config.events = 500_000
  config.decay_card = decay_card_d0_pi_munu
  config.cross_section = :default
end

exMC_dp_pi0_enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_DpToPi0Eneu"
  config.related_dataset = data_3773
  config.events = 500_000
  config.decay_card = decay_card_dp_pi0_enu
  config.cross_section = :default
end

exMC_dp_pi0_munu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_DpToPi0Munu"
  config.related_dataset = data_3773
  config.events = 500_000
  config.decay_card = decay_card_dp_pi0_munu
  config.cross_section = :default
end

### Single-tag analyses (TagAnalysis) ###

# ------------------------------------------------------------------
# Channel 1: tag a D0 (charm -1), signal = D0 -> pi- e+ nu_e
# ------------------------------------------------------------------
alg_d0_e = TagAnalysis.new("D0TagPiEneu")
alg_d0_e.set_header(["D0TagPiEneuAlg/D0TagPiEneu.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .note(:companion_selection, "companion-side (tag) selection follows the dynamics paper")
        .with_decay_card(decay_card_d0_pi_enu)

# Tag side: D0 single tag, charm pinned to -1, hadronic tag modes
alg_d0_e.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPi0Pi0, :D0toKPiPiPiPi0
  t.charm(-1)
end

# Signal side: exactly one pi- and one e+, net charge 0, neutrino missing
alg_d0_e.signal_side do |s|
  s.charged(pim: 1, ep: 1)
  s.require_charge 0
  s.missing :nu_e
end

# 4C kinematic fit
alg_d0_e.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_e.apply

# ------------------------------------------------------------------
# Channel 2: tag a D0 (charm -1), signal = D0 -> pi- mu+ nu_mu
# ------------------------------------------------------------------
alg_d0_mu = TagAnalysis.new("D0TagPiMunu")
alg_d0_mu.set_header(["D0TagPiMunuAlg/D0TagPiMunu.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .note(:companion_selection, "companion-side (tag) selection follows the dynamics paper")
         .with_decay_card(decay_card_d0_pi_munu)

alg_d0_mu.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPi0Pi0, :D0toKPiPiPiPi0
  t.charm(-1)
end

# Signal side: exactly one pi- and one mu+, net charge 0, neutrino missing
alg_d0_mu.signal_side do |s|
  s.charged(pim: 1, mup: 1)
  s.require_charge 0
  s.missing :nu_mu
end

alg_d0_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_mu.apply

# ------------------------------------------------------------------
# Channel 3: tag a D+ (charm -1), signal = D+ -> pi0 e+ nu_e (pi0 -> gamma gamma)
# ------------------------------------------------------------------
alg_dp_e = TagAnalysis.new("DpTagPi0Eneu")
alg_dp_e.set_header(["DpTagPi0EneuAlg/DpTagPi0Eneu.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .note(:companion_selection, "companion-side (tag) selection follows the dynamics paper")
        .with_decay_card(decay_card_dp_pi0_enu)

# Tag side: D+ single tag, charm pinned to -1, hadronic tag modes
alg_dp_e.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm(-1)
end

# Signal side: two photons + one e+, net charge +1, neutrino missing
alg_dp_e.signal_side do |s|
  s.photons 2
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

# 4C fit + pi0 mass constraint on the two-photon system
alg_dp_e.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_dp_e.apply

# ------------------------------------------------------------------
# Channel 4: tag a D+ (charm -1), signal = D+ -> pi0 mu+ nu_mu (pi0 -> gamma gamma)
# ------------------------------------------------------------------
alg_dp_mu = TagAnalysis.new("DpTagPi0Munu")
alg_dp_mu.set_header(["DpTagPi0MunuAlg/DpTagPi0Munu.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .note(:companion_selection, "companion-side (tag) selection follows the dynamics paper")
         .with_decay_card(decay_card_dp_pi0_munu)

alg_dp_mu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm(-1)
end

# Signal side: two photons + one mu+, net charge +1, neutrino missing
alg_dp_mu.signal_side do |s|
  s.photons 2
  s.charged(mup: 1)
  s.require_charge 1
  s.missing :nu_mu
end

alg_dp_mu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_dp_mu.apply

### Execute each channel over its data, inclusive MC and matching exclusive MC ###
root_files_d0_e  = alg_d0_e.execute_on([data_3773, incMC_3773, exMC_d0_pi_enu])
root_files_d0_mu = alg_d0_mu.execute_on([data_3773, incMC_3773, exMC_d0_pi_munu])
root_files_dp_e  = alg_dp_e.execute_on([data_3773, incMC_3773, exMC_dp_pi0_enu])
root_files_dp_mu = alg_dp_mu.execute_on([data_3773, incMC_3773, exMC_dp_pi0_munu])