# BOSS Ruby DSL for arXiv:2007.10563v1
# "Measurements of absolute branching fractions of D0(+) -> KKpipi decays"
#
# Analysis: e+e- -> psi(3770) -> D Dbar pairs at sqrt(s) = 3.773 GeV
#   Double-tag method to measure absolute BFs of 9 D -> KKpipi decay modes.
#   D0 modes: K+K-pi0pi0, KS0KS0pi+pi-, KS0K-pi+pi0, KS0K+pi-pi0
#   D+ modes: K+K-pi+pi0, KS0K+pi0pi0, KS0K-pi+pi+, KS0K+pi+pi-, KS0KS0pi+pi0
#   Tag D0-bar: K+pi-, K+pi-pi0, K+pi-pi-pi+
#   Tag D-: K+pi-pi-, KS0pi-, K+pi-pi-pi0, KS0pi-pi0, KS0pi+pi-pi-, K+K-pi-

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ===========================================================================
# Shared decay cards and MC configurations
# ===========================================================================

# D0 -> K+ K- pi0 pi0  (decay card with D0-bar tags)
decay_card_d0_kk_pi0pi0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K+ K- pi0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    0.400 K+ pi- PHSP;
    0.350 K+ pi- pi0 PHSP;
    0.250 K+ pi- pi- pi+ PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
End
DECAYCARD

# D0 -> K_S0 K_S0 pi+ pi-
decay_card_d0_ksks_pipi = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K_S0 K_S0 pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    0.400 K+ pi- PHSP;
    0.350 K+ pi- pi0 PHSP;
    0.250 K+ pi- pi- pi+ PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
End
DECAYCARD

# D+ -> K+ K- pi+ pi0  (decay card with D- tags)
decay_card_dp_kk_pipi0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K+ K- pi+ pi0 PHSP;
    Enddecay

    Decay D-
    0.350 K+ pi- pi- PHSP;
    0.250 K_S0 pi- PHSP;
    0.400 K+ pi- pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay
End
DECAYCARD

# D+ -> K_S0 K_S0 pi+ pi0
decay_card_dp_ksks_pipi0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K_S0 K_S0 pi+ pi0 PHSP;
    Enddecay

    Decay D-
    0.350 K+ pi- pi- PHSP;
    0.250 K_S0 pi- PHSP;
    0.400 K+ pi- pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay
End
DECAYCARD

# ===========================================================================
# Exclusive MC samples (representative subset of the 9 channels)
# ===========================================================================
exMC_d0_kk_pi0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "d0_kk_pi0pi0_mc"
  config.related_dataset = psi3770_data
  config.events = 500_000
  config.decay_card = decay_card_d0_kk_pi0pi0
  config.cross_section = :default
end

exMC_d0_ksks_pipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "d0_ksks_pipi_mc"
  config.related_dataset = psi3770_data
  config.events = 500_000
  config.decay_card = decay_card_d0_ksks_pipi
  config.cross_section = :default
end

exMC_dp_kk_pipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "dp_kk_pipi0_mc"
  config.related_dataset = psi3770_data
  config.events = 500_000
  config.decay_card = decay_card_dp_kk_pipi0
  config.cross_section = :default
end

exMC_dp_ksks_pipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "dp_ksks_pipi0_mc"
  config.related_dataset = psi3770_data
  config.events = 500_000
  config.decay_card = decay_card_dp_ksks_pipi0
  config.cross_section = :default
end

all_exMC = [exMC_d0_kk_pi0pi0, exMC_d0_ksks_pipi, exMC_dp_kk_pipi0, exMC_dp_ksks_pipi0]

# ===========================================================================
# TagAnalysis 1: D0 -> K+ K- pi0 pi0  (tag: D0-bar -> K+pi-, K+pi-pi0, K+pi-pi-pi+)
#   Signal side: K+, K-, pi0, pi0 (4 photons from 2 pi0)
# ===========================================================================
alg_d0_kk_pi0pi0 = TagAnalysis.new("D0toKKPi0Pi0")
alg_d0_kk_pi0pi0.set_header(["D0toKKPi0Pi0Alg/D0toKKPi0Pi0.h"])
                 .set_constant({"ECMS" => [:double, 3.773]})
                 .with_decay_card(decay_card_d0_kk_pi0pi0)
                 .note(:signal_channels_not_all, "9 signal modes in paper. Representative subset expressed in DSL due to BOSS tag vocabulary limitations. Remaining signal modes: D0->KS0K-pi+pi0, D0->KS0K+pi-pi0, D+->KS0K+pi0pi0, D+->KS0K-pi+pi+, D+->KS0K+pi+pi-")

alg_d0_kk_pi0pi0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

alg_d0_kk_pi0pi0.signal_side do |s|
  s.photons 4
  s.charged(kp: 1, km: 1)
  s.require_charge 0
end

alg_d0_kk_pi0pi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_d0_kk_pi0pi0.apply

# ===========================================================================
# TagAnalysis 2: D0 -> K_S0 K_S0 pi+ pi-  (tag: D0-bar CF modes)
#   Signal side: KS0->pi+pi-, KS0->pi+pi-, pi+, pi- (no photons)
# ===========================================================================
alg_d0_ksks_pipi = TagAnalysis.new("D0toKSKSPiPi")
alg_d0_ksks_pipi.set_header(["D0toKSKSPiPiAlg/D0toKSKSPiPi.h"])
                .set_constant({"ECMS" => [:double, 3.773]})
                .with_decay_card(decay_card_d0_ksks_pipi)
                .note(:ks0_reconstruction, "K_S0 reconstructed via secondary vertex fit on remaining pi+pi- pairs. KS0 mass window (0.486,0.510) GeV/c2 and decay-length > 2sigma vertex resolution applied in ROOT")

alg_d0_ksks_pipi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

alg_d0_ksks_pipi.signal_side do |s|
  s.charged(pip: 3, pim: 3)
  s.require_charge 0
end

alg_d0_ksks_pipi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_ksks_pipi.apply

# ===========================================================================
# TagAnalysis 3: D+ -> K+ K- pi+ pi0  (tag: D- -> K+pi-pi-, KS0pi-, K+pi-pi-pi0)
#   Signal side: K+, K-, pi+, pi0 -> gamma gamma
# ===========================================================================
alg_dp_kk_pipi0 = TagAnalysis.new("DptoKKPiPi0")
alg_dp_kk_pipi0.set_header(["DptoKKPiPi0Alg/DptoKKPiPi0.h"])
               .set_constant({"ECMS" => [:double, 3.773]})
               .with_decay_card(decay_card_dp_kk_pipi0)

alg_dp_kk_pipi0.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0
  t.charm -1
end

alg_dp_kk_pipi0.signal_side do |s|
  s.photons 2
  s.charged(kp: 1, km: 1, pip: 1)
  s.require_charge 1
end

alg_dp_kk_pipi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_dp_kk_pipi0.apply

# ===========================================================================
# TagAnalysis 4: D+ -> K_S0 K_S0 pi+ pi0  (tag: D- CF modes)
#   Signal side: KS0->pi+pi-, KS0->pi+pi-, pi+, pi0->gamma gamma
# ===========================================================================
alg_dp_ksks_pipi0 = TagAnalysis.new("DptoKSKSPiPi0")
alg_dp_ksks_pipi0.set_header(["DptoKSKSPiPi0Alg/DptoKSKSPiPi0.h"])
                 .set_constant({"ECMS" => [:double, 3.773]})
                 .with_decay_card(decay_card_dp_ksks_pipi0)

alg_dp_ksks_pipi0.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0
  t.charm -1
end

alg_dp_ksks_pipi0.signal_side do |s|
  s.photons 2
  s.charged(pip: 3, pim: 3)
  s.require_charge 1
end

alg_dp_ksks_pipi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_dp_ksks_pipi0.apply

# Execute all tag analyses
alg_d0_kk_pi0pi0.execute_on([psi3770_data, psi3770_incMC, exMC_d0_kk_pi0pi0])
alg_d0_ksks_pipi.execute_on([psi3770_data, psi3770_incMC, exMC_d0_ksks_pipi])
alg_dp_kk_pipi0.execute_on([psi3770_data, psi3770_incMC, exMC_dp_kk_pipi0])
alg_dp_ksks_pipi0.execute_on([psi3770_data, psi3770_incMC, exMC_dp_ksks_pipi0])