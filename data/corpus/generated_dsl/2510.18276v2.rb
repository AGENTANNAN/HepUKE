# =============================================================================
#  psi(3770) -> D Dbar double-tag analysis: absolute branching fractions of
#  D0 / D+ into K K K pi final states   (BOSS / TagAnalysis part)
#  Both D mesons are taken from the pre-stored DTag collection (DTagTool),
#  so the whole selection lives in TagAnalysis (no Selection object).
# =============================================================================

### Dataset description ###
psipp_data  = DatasetManager.real_data.find("712_3773")     # psi(3770) data, sqrt(s) = 3.773 GeV, ~20.3 fb^-1
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC sample

# --- decay card for the D0 D0bar signal process (K pi / K_S0 tag chains) ---
decay_card_D0D0bar = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# --- decay card for the D+ D- signal process (K pi pi / K_S0 tag chains) ---
decay_card_DpDm = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 K- pi+ pi+ PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# --- exclusive MC samples: 200k events for each of the two D Dbar processes ---
exMC_D0D0bar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0D0bar"
  config.related_dataset = psipp_data
  config.events          = 200000
  config.decay_card      = decay_card_D0D0bar
  config.cross_section   = :default
end

exMC_DpDm = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpDm"
  config.related_dataset = psipp_data
  config.events          = 200000
  config.decay_card      = decay_card_DpDm
  config.cross_section   = :default
end

# =============================================================================
#  D0 analysis:  anti-D0 (tag)  x  D0 -> K K K pi (signal side, DT measurement)
# =============================================================================
alg_D0 = TagAnalysis.new("D0ToKKKpiDT")
alg_D0.set_header(["D0ToKKKpiDTAlg/D0ToKKKpiDT.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_D0D0bar)

# Tag side 1: anti-D0 -> K+ pi-, K+ pi- pi0, K+ pi- pi- pi+   (charm = -1 pins the D0bar side)
alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm(-1)
end

# Tag side 2: the D0 side with the same modes, charm = +1, candidates ranked by invariant mass
alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm(1)
  t.rank_by :inv
end

# Kinematic fit: 4-momentum conservation, both D tag masses constrained to the
# nominal D0 mass, chi2 < 200.  Tag mBC / deltaE are stored unconditionally
# (store-not-cut) and windowed later in the ROOT analysis.
alg_D0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# Signal D0 -> K K K pi mode list, reconstructed from the tracks/showers not used
# by the tag (K_S0 -> pi+ pi-, phi -> K+ K-).  These signal mode lists have no
# DTagAlg channel name, so they cannot be declared as tag modes.
alg_D0.note(:signal_mode_reconstruction,
  "signal D0 is reconstructed from the tracks left unused by the tag in the K K K pi modes " \
  "K_S0 K+ K- pi0, K_S0 K_S0 K- pi+, K_S0 K_S0 K+ pi-, K+ K- K- pi+, together with the " \
  "phi sub-modes phi K_S0 pi0 and phi K- pi+ (phi -> K+ K-, K_S0 -> pi+ pi-); " \
  "only the D Dbar pair is present in the final state at the psi(3770), no extra particles")

alg_D0.apply
alg_D0.execute_on([psipp_data, psipp_incMC, exMC_D0D0bar])

# =============================================================================
#  D+ analysis:  D- (tag)  x  D+ -> K K K pi (signal side, DT measurement)
# =============================================================================
alg_Dp = TagAnalysis.new("DpToKKKpiDT")
alg_Dp.set_header(["DpToKKKpiDTAlg/DpToKKKpiDT.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_DpDm)

# Tag side 1: D- -> K+ pi- pi-, K_S0 pi-, K+ pi- pi- pi0, K_S0 pi- pi0, K_S0 pi- pi- pi+, K+ K- pi-
alg_Dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm(-1)
end

# Tag side 2: charge-conjugate D+ tag with the same modes, ranked by invariant mass
alg_Dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm(1)
  t.rank_by :inv
end

# Kinematic fit: 4-momentum conservation, both D tag masses constrained to the
# nominal D+ mass, chi2 < 200.
alg_Dp.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"D+")
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:"D+")
  f.chi2_cut 200
end

# Signal D+ -> K K K pi mode list, reconstructed from the tracks not used by the
# tag (K_S0 -> pi+ pi-, phi -> K+ K-); no DTagAlg channel name exists for them.
alg_Dp.note(:signal_mode_reconstruction,
  "signal D+ is reconstructed from the tracks left unused by the tag in " \
  "K_S0 K+ K- pi+ and phi K_S0 pi+ (phi -> K+ K-, K_S0 -> pi+ pi-); " \
  "only the D Dbar pair is present in the final state at the psi(3770), no extra particles")

alg_Dp.apply
alg_Dp.execute_on([psipp_data, psipp_incMC, exMC_DpDm])