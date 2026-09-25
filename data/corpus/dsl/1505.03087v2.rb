# ============================================================================
# arXiv:1505.03087v2 (BESIII)
# Search for the FCNC decay D0 -> gamma gamma using a double-tag technique at
# sqrt(s) = 3.773 GeV (psi(3770) -> D0 D0bar), with an improved measurement of
# B(D0 -> pi0 pi0) from the same sample.
#
# Data: 2.92 fb^-1 collected at the psi(3770) peak.
# Tag side: hadronically decaying D0bar in five modes,
#   K+ pi-, K+ pi- pi0, K+ pi- pi+ pi-, K+ pi- pi+ pi- pi0, K+ pi- pi0 pi0.
# Signal side: D0 -> gamma gamma (or D0 -> pi0 pi0 for the by-product
# measurement) from the showers not used by the tag.
#
# Tag-based analysis: TagAnalysis (DTagAlg / DTagTool), no Selection object.
# BOSS part only: decay cards, exclusive MC, tag declarations and the
# kinematic fit.
# ============================================================================

### Dataset description ###
# psi(3770) data and inclusive MC (tag-reformed generic sample) at 3.773 GeV
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ----------------------------------------------------------------------------
# Signal decay cards: psi(3770) -> D0 D0bar, one D decaying to a hadronic tag
# mode and the other to the signal final state.
# ----------------------------------------------------------------------------
# D0 -> gamma gamma (the signal mode); the tag side decays generically
decay_card_gammagamma = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# D0 -> pi0 pi0 (the by-product measurement) with pi0 -> gamma gamma
decay_card_pi0pi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# ----------------------------------------------------------------------------
# Exclusive MC. The signal samples consist of e+e- -> psi(3770) -> D0 D0bar
# events in which one D decays into a hadronic tag mode or into the signal
# final state, while the other D decays without restriction.
# ----------------------------------------------------------------------------
exMC_gammagamma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3770_D0togammagamma"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_gammagamma
  config.cross_section   = :default
end

exMC_pi0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3770_D0topi0pi0"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_pi0pi0
  config.cross_section   = :default
end

### Tag analysis (BOSS) ###
# Signal D0 -> gamma gamma (the tag side is the anti-D0)
alg_name_gg = "D0ToGammaGammaDTag"
alg_gg = TagAnalysis.new(alg_name_gg)
alg_gg.set_header(["#{alg_name_gg}Alg/#{alg_name_gg}.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .set_alias({"std::vector<double>" => "VDouble"})
      .with_decay_card(decay_card_gammagamma)

# Tag side 1: the anti-D0 in the signal-side flavour combination
alg_gg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPiPiPi0, :D0toKPiPi0Pi0
  t.charm -1
  t.rank_by :inv
end

# Tag side 2: the D0 recoiling against the signal D0 (double tag)
alg_gg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPiPiPi0, :D0toKPiPi0Pi0
  t.charm 1
  t.rank_by :inv
end

# Signal side: D0 -> gamma gamma from the two most energetic showers that are
# not used by either tag
alg_gg.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0        # > 10 deg from any charged track (QED Bhabha veto)
  s.min_photon_energy 0.025      # barrel shower energy floor (GeV)
end

alg_gg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:D0)   # signal D0 mass constraint
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)            # tag-side 6C lines
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200                 # loose BOSS-pass cut; the tight selection is applied in ROOT
end

# Signal D0 -> pi0 pi0 (by-product measurement)
alg_name_pipi = "D0ToPi0Pi0DTag"
alg_pipi = TagAnalysis.new(alg_name_pipi)
alg_pipi.set_header(["#{alg_name_pipi}Alg/#{alg_name_pipi}.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .set_alias({"std::vector<double>" => "VDouble"})
        .with_decay_card(decay_card_pi0pi0)

alg_pipi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPiPiPi0, :D0toKPiPi0Pi0
  t.charm -1
  t.rank_by :inv
end

alg_pipi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPiPiPi0, :D0toKPiPi0Pi0
  t.charm 1
  t.rank_by :inv
end

# Signal side: D0 -> pi0 pi0 -> four photons from the showers not used by the tags
alg_pipi.signal_side do |s|
  s.photons 4
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_pipi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # first pi0
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # second pi0
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# Notes for BOSS-side criteria without a corresponding DSL construct
alg_gg
  .note(:background_veto, "For any DT combination including the D0 -> K+ pi- tag the " \
        "dominant background is the doubly-radiative Bhabha process " \
        "e+e- -> e+e- gamma gamma. It is suppressed by requiring the angle between " \
        "each signal photon candidate and any charged track to be greater than " \
        "10 degrees, which removes 93% of the QED background.")
  .note(:background_veto, "The dominant peaking background in the DeltaE(gamma gamma) " \
        "signal region is D0 -> pi0 pi0. A pi0 veto rejects events in which one of " \
        "the two signal photons can be combined with any other photon in the event " \
        "to form a pi0; this rejects 82% of the D0 -> pi0 pi0 background while " \
        "keeping 88% of the signal.")
  .note(:tag_side_selection, "The tag-side selection is applied on the stored tag " \
        "observables: 1.847 < mBC(tag) < 1.883 GeV/c^2 and |DeltaE(tag)| < 0.1 GeV. " \
        "In events with multiple tag candidates the one candidate per mode with " \
        "reconstructed energy closest to the beam energy is chosen, and among " \
        "multiple DT candidates the combination whose average " \
        "mBC = (mBC(tag) + mBC(gamma gamma))/2 is closest to the known D0 mass. " \
        "The tight signal-side requirement is |DeltaE(gamma gamma)| < 0.25 GeV " \
        "(mBC(gamma gamma) > 1.85 GeV/c^2); these windows are applied at the ROOT " \
        "stage.")
  .note(:efficiency_curve, "The DT-to-ST ratio method provides an absolute " \
        "normalisation independent of the integrated luminosity and of the " \
        "D0 D0bar production cross section; most tag-reconstruction systematic " \
        "uncertainties cancel in the ratio. The remaining independent sources are " \
        "photon reconstruction (2.0%), the signal-side mBC requirement (3.1%) and " \
        "the ST yields (1.0%), giving a total of 3.8%.")
  .note(:fit_model_variation, "The two-dimensional unbinned maximum-likelihood fit " \
        "of DeltaE(gamma gamma) versus DeltaE(tag) uses PDFs extracted from MC: " \
        "signal and D0 -> pi0 pi0 shapes from MC, a flat continuum component with " \
        "floating normalisation, and a Crystal Ball line plus Gaussian in " \
        "DeltaE(tag) with a second-order exponential polynomial in DeltaE(gamma gamma) " \
        "for other D0 D0bar decays. The systematic variations (fit ranges, " \
        "background shapes, D0 -> pi0 pi0 normalisation and shape, signal shape) " \
        "are part of the ROOT fitting stage.")

alg_pipi
  .note(:background_veto, "For any DT combination including the D0 -> K+ pi- tag the " \
        "dominant background is the doubly-radiative Bhabha process " \
        "e+e- -> e+e- gamma gamma, suppressed by requiring the angle between each " \
        "signal photon candidate and any charged track to be greater than 10 degrees.")
  .note(:tag_side_selection, "The tag-side selection is applied on the stored tag " \
        "observables: 1.847 < mBC(tag) < 1.883 GeV/c^2 and |DeltaE(tag)| < 0.1 GeV. " \
        "For the D0 -> pi0 pi0 measurement the two pi0 candidates are formed from " \
        "the photons not used by the tag, the pair giving the smallest " \
        "|DeltaE(pi0 pi0)| is selected, and -0.070 < DeltaE(pi0 pi0) < +0.075 GeV " \
        "is required before fitting mBC(pi0 pi0).")
  .note(:efficiency_curve, "The D0 -> pi0 pi0 yield is normalised to the number of " \
        "D0 D0bar pairs, N_D0D0bar = L * sigma(e+e- -> psi(3770) -> D0 D0bar) with " \
        "L the integrated luminosity and sigma = (3.607 +/- 0.017 +/- 0.056) nb; " \
        "the uncertainty on N_D0D0bar contributes 1.9% to the total 3.6% " \
        "systematic uncertainty of the branching fraction.")

# Render the tag specs (no Selection argument: the tag declarations are the selection)
alg_gg.apply
alg_pipi.apply

# Execute on real data, inclusive MC and the signal exclusive MC
root_files_gg   = alg_gg.execute_on([psi3770_data, psi3770_incMC, exMC_gammagamma])
root_files_pipi = alg_pipi.execute_on([psi3770_data, psi3770_incMC, exMC_pi0pi0])
