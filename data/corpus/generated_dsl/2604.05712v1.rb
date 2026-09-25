# =====================================================================
# BESIII BOSS DSL — psi(3770) quantum-correlated D0 D0bar, CKM angle gamma
# via an unbinned optimal Fourier method (strong-phase measurement,
# following JHEP 06 (2025) 086).
#
# This is a tag-based analysis: each event is a psi(3770) -> D0 D0bar pair,
# one D reconstructed from the pre-stored DTag candidates (the "tag"),
# the other reconstructed in the signal mode (K_S0 h+ h-).  The two
# signal modes are independent double-tag selections, so each gets its own
# TagAnalysis object.  No Selection / select_track / pid / photons are used —
# the tag already carries selected, PID'd tracks and showers.
# =====================================================================

### Dataset description ###
# psi(3770): real data (2010-2011 = round03/04, 2021-2022 = round15/16)
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card (signal mode I): D0 -> K_S0 pi+ pi- modeled with a Dalitz-plot
# amplitude; the opposite D0bar -> K+ pi- (phase space); K_S0 -> pi+ pi-.
decay_card_kspipi = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0000 K_S0 pi+ pi- D_DALITZ;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi- PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Decay card (signal mode II): D0 -> K_S0 K+ K- modeled with a Dalitz-plot
# amplitude; the opposite D0bar -> K+ pi- (phase space); K_S0 -> pi+ pi-.
decay_card_kskk = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0000 K_S0 K+ K- D_DALITZ;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi- PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC: 500k events for each of the two signal modes
exMC_kspipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "D0toKsPiPi_3773_exMC"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_kspipi
  config.cross_section   = :default
end

exMC_kskk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "D0toKsKK_3773_exMC"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_kskk
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# =====================================================================
# Double-tag selection I : signal D0 -> K_S0 pi+ pi-
# =====================================================================
alg_kspipi = TagAnalysis.new("D0ToKsPiPiTag")
alg_kspipi.set_header(["D0ToKsPiPiTagAlg/D0ToKsPiPiTag.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .with_decay_card(decay_card_kspipi)

# signal side: D0 (charm +1) reconstructed as K_S0 pi+ pi-
alg_kspipi.tag_side(:D0) do |t|
  t.modes :D0toKsPiPi
  t.charm 1
end

# tag side: the other D (charm -1) reconstructed in flavor-specific, CP-even,
# CP-odd and self-conjugate modes; double-tag candidates ranked by invariant mass
alg_kspipi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKK, :D0toPiPi,
          :D0toKsPi0Pi0, :D0toPiPiPi0, :D0toKsPi0, :D0toKsEta,
          :D0toKsEtaPrime, :D0toKsOmega, :D0toKsPiPi, :D0toKsKK
  t.charm  -1
  t.rank_by :inv
end

# kinematic fit: four-momentum conservation + D mass constraints, chi2 < 200
alg_kspipi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_kspipi
  .note(:final_state_selection, "Final-state particle selection follows the previous BESIII binned strong-phase measurement, JHEP 06 (2025) 086.")
  .note(:ks_kl_mass_constraint, "The kinematic fit also constrains the K_S0 and K_L0 masses to their nominal values; these constraints act on the tag sub-decays and are not expressible at the tag-group level of the DSL.")
  .note(:kl_missing_mass, "For the K_L0 h+ h- signal and for K_L0 tags the K_L0 is treated as missing; the missing-mass-squared is used at the ROOT stage.")
  .note(:background_veto, "Dominant peaking backgrounds (D -> pi+ pi- h+ h- for K_S0 h+ h-, and D -> K_S0(->pi0 pi0) h+ h- for K_L0 h+ h-) are handled in the ROOT fit.")
  .apply

alg_kspipi.execute_on([data_3773, incMC_3773, exMC_kspipi])

# =====================================================================
# Double-tag selection II : signal D0 -> K_S0 K+ K-
# =====================================================================
alg_kskk = TagAnalysis.new("D0ToKsKKTag")
alg_kskk.set_header(["D0ToKsKKTagAlg/D0ToKsKKTag.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .with_decay_card(decay_card_kskk)

# signal side: D0 (charm +1) reconstructed as K_S0 K+ K-
alg_kskk.tag_side(:D0) do |t|
  t.modes :D0toKsKK
  t.charm 1
end

# tag side: the other D (charm -1) reconstructed in the same set of
# flavor-specific, CP-even, CP-odd and self-conjugate modes
alg_kskk.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKK, :D0toPiPi,
          :D0toKsPi0Pi0, :D0toPiPiPi0, :D0toKsPi0, :D0toKsEta,
          :D0toKsEtaPrime, :D0toKsOmega, :D0toKsPiPi, :D0toKsKK
  t.charm  -1
  t.rank_by :inv
end

# kinematic fit: four-momentum conservation + D mass constraints, chi2 < 200
alg_kskk.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_kskk
  .note(:final_state_selection, "Final-state particle selection follows the previous BESIII binned strong-phase measurement, JHEP 06 (2025) 086.")
  .note(:ks_kl_mass_constraint, "The kinematic fit also constrains the K_S0 and K_L0 masses to their nominal values; these constraints act on the tag sub-decays and are not expressible at the tag-group level of the DSL.")
  .note(:kl_missing_mass, "For the K_L0 h+ h- signal and for K_L0 tags the K_L0 is treated as missing; the missing-mass-squared is used at the ROOT stage.")
  .note(:background_veto, "Dominant peaking backgrounds (D -> pi+ pi- h+ h- for K_S0 h+ h-, and D -> K_S0(->pi0 pi0) h+ h- for K_L0 h+ h-) are handled in the ROOT fit.")
  .apply

alg_kskk.execute_on([data_3773, incMC_3773, exMC_kskk])