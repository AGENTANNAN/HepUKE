# BOSS Ruby DSL for arXiv:2007.07959v1
# "Improved model-independent determination of the strong-phase difference between
#  D0 and anti-D0 -> K_{S,L}^0 K+ K- decays"
#
# Analysis: e+e- -> psi(3770) -> D0 anti-D0 at sqrt(s) = 3.773 GeV
#   Double-tag method exploiting quantum correlations of the D0 anti-D0 pair.
#   Signal channels: D -> K_S^0 K+ K-  and  D -> K_L^0 K+ K-
#   Multiple tag types: flavor, CP-even, CP-odd, mixed-CP
#   Dalitz plot binned in N=2,3,4 equal-delta_delta_D bins for strong-phase extraction.
#   ROOT-level: unbinned max likelihood fit to extract c_i and s_i parameters.
#   All tag-side mBC/deltaE cuts, DT yields, and Dalitz-plot binning done in ROOT.

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ===========================================================================
# Decay card: D0 -> K_S0 K+ K-  vs various tag modes
# Tag modes list:
#   Flavor:  D0 -> K- pi+, D0 -> K- pi+ pi0, D0 -> K- e+ nu_e
#   CP-even: D0 -> K+ K-, D0 -> pi+ pi-, D0 -> K_S0 pi0, D0 -> pi+ pi- pi0
#   CP-odd:  D0 -> K_S0 pi0, D0 -> K_S0 eta, D0 -> K_S0 omega, D0 -> K_S0 eta'
#   Mixed:   D0 -> K_S0 pi+ pi-, D0 -> K_S0 K+ K-
#   Plus KL0 tag modes where KL0 is reconstructed via EMC shower.
# ===========================================================================
decay_card_kskk = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K_S0 K+ K- PHSP;
    Enddecay

    Decay anti-D0
    0.250 K+ pi- PHSP;
    0.250 K+ pi- pi0 PHSP;
    0.125 K+ K- PHSP;
    0.125 pi+ pi- PHSP;
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

# Decay card for D0 -> K_L0 K+ K-  (KL0 not explicitly generated, handled as missing)
# For simulation, KL0 is produced via D0 -> K0bar K+ K- with K0bar -> K_L0
# In real analysis, KL0 inferred from EMC shower cluster + missing mass
decay_card_klkk = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K_L0 K+ K- PHSP;
    Enddecay

    Decay anti-D0
    0.250 K+ pi- PHSP;
    0.250 K+ pi- pi0 PHSP;
    0.125 K+ K- PHSP;
    0.125 pi+ pi- PHSP;
    0.250 K+ pi- pi- pi+ PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
End
DECAYCARD

# Exclusive MC for D0 -> K_S0 K+ K-
exMC_kskk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "d0_kskk_signal_mc"
  config.related_dataset = psi3770_data
  config.events = 1_000_000
  config.decay_card = decay_card_kskk
  config.cross_section = :default
end

# Exclusive MC for D0 -> K_L0 K+ K-
exMC_klkk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "d0_klkk_signal_mc"
  config.related_dataset = psi3770_data
  config.events = 1_000_000
  config.decay_card = decay_card_klkk
  config.cross_section = :default
end

# ===========================================================================
# TagAnalysis: D0 -> K_S0 K+ K- signal vs D0-bar tags
#   Signal side: K_S0 -> pi+ pi- + K+ K-
# ===========================================================================
alg_kskk = TagAnalysis.new("D0toKSKK")
alg_kskk.set_header(["D0toKSKKAlg/D0toKSKK.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .with_decay_card(decay_card_kskk)
        .note(:tag_modes_full_list, "Full tag mode list from paper: flavor (K-pi+, K-pi+pi0, K-e+nu), CP-even (K+K-, pi+pi-, KSpi0, pi+pi-pi0, KLpi0, KLetap, KLeta), CP-odd (KSpi0, KSeta, KSomega, KSetap), mixed-CP (KSpi+pi-, KLpi+pi-, KSK+K-). Only modes available in BOSS tag-mode vocabulary are declared; others applied in ROOT")
        .note(:kl0_tag_reconstruction, "K_L0 tag modes reconstruct KL0 via EMC shower + missing mass; KL0 is treated as missing particle in kinematic fit. Not expressible in TagAnalysis DSL")
        .note(:dalitz_plot_analysis, "Dalitz plot binned in N=2,3,4 equal-delta_delta_D bins. Bin yields extracted via sideband estimation. c_i and s_i strong-phase parameters from unbinned max likelihood fit. ROOT-level only")
        .note(:dcs_correction, "DCS contamination correction applied to flavor-tag yields using amplitude model from BaBar. Bin migration correction via unfolding matrix. ROOT-level")
        .note(:kinematic_fit_improvement, "Kinematic fit constraining D0 mass improves m_pm^2 resolution by 35-40%. For KL0 modes, missing KL0 treatment + mass constraint applied. ROOT-level")

# Tag side: D0-bar reconstructed via available flavor + CP tag modes
# Note: many tag modes from the paper are not in the frozen BOSS tag vocabulary
alg_kskk.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKK, :D0toPiPi, :D0toKsPi0, :D0toPiPiPi0, :D0toKsPiPi
  t.charm -1
end

# Signal side: D0 -> K_S0 K+ K- with K_S0 -> pi+ pi-
# K_S0 reconstructed via secondary vertex on remaining tracks
alg_kskk.signal_side do |s|
  s.charged(kp: 1, km: 1, pip: 1, pim: 1)
  s.require_charge 0
end

alg_kskk.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_kskk.apply
alg_kskk.execute_on([psi3770_data, psi3770_incMC, exMC_kskk])

# ===========================================================================
# TagAnalysis: D0 -> K_L0 K+ K- signal vs D0-bar tags
#   Signal side: K+ K- + missing K_L0
#   K_L0 is massive (mass from PDG used in missing-particle kinematic fit)
# ===========================================================================
alg_klkk = TagAnalysis.new("D0toKLKK")
alg_klkk.set_header(["D0toKLKKAlg/D0toKLKK.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .with_decay_card(decay_card_klkk)
        .note(:tag_modes_full_list, "Full tag mode list from paper: see note on D0toKSKK algorithm")
        .note(:kl0_emc_shower, "K_L0 reconstructed via EMC shower cluster + missing momentum. 3-momentum from shower position, mass fixed to nominal KL0 for kinematic fit. ROOT-level")
        .note(:dalitz_plot_analysis, "Dalitz plot analysis for KL0 K+ K-: binning, sideband estimation, unfolding. ROOT-level")

# Tag side: D0-bar reconstructed via available modes
alg_klkk.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKK, :D0toPiPi, :D0toKsPi0, :D0toPiPiPi0, :D0toKsPiPi
  t.charm -1
end

# Signal side: D0 -> K_L0 K+ K- with K_L0 as missing particle (massive)
alg_klkk.signal_side do |s|
  s.charged(kp: 1, km: 1)
  s.require_charge 0
  s.missing :K_L0
end

alg_klkk.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_klkk.apply
alg_klkk.execute_on([psi3770_data, psi3770_incMC, exMC_klkk])