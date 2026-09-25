# =============================================================================
# BOSS-side HepScript DSL
# psi(3770) -> D0 D0bar, one D tagged hadronically, the other reconstructed as
# the signal:  D0 -> gamma gamma  (FCNC search)   and   D0 -> pi0 pi0  (BF)
# Tag-based analysis -> TagAnalysis (no Selection; the tag carries its own
# selected / PID'd tracks and showers, the signal side comes from the tag's
# unused showers).
# =============================================================================

### Datasets ###
data_3773  = DatasetManager.real_data.find("712_3773")      # psi(3770) real data, 2.92 fb^-1
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # inclusive MC at 3.773 GeV

### Decay cards (EvtGen syntax) ###
# Signal mode I : D0 -> gamma gamma, tag D0bar -> K+ pi-
decay_card_gg = <<~DECAYCARD
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

# Signal mode II : D0 -> pi0 pi0 (pi0 -> gamma gamma), tag D0bar -> K+ pi-
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

### Exclusive signal MC (500k events per signal mode) ###
exMC_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_D0toGG"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_gg
  config.cross_section   = :default
end

exMC_pi0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_D0topi0pi0"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_pi0pi0
  config.cross_section   = :default
end

### ---------------------------------------------------------------------------
### Signal mode I : D0 -> gamma gamma
### ---------------------------------------------------------------------------
alg_gg = TagAnalysis.new("D0ToGG")
alg_gg.set_header(["D0ToGGAlg/D0ToGG.h"])                        # algorithm header
      .set_constant({ "ECMS" => [:double, 3.773] })              # sqrt(s) = 3.773 GeV
      .note(:tag_candidate_selection,
            "per tag mode one candidate is kept: the one whose reconstructed D0 " \
            "energy is closest to the beam energy (DTagAlg internal ranking, not " \
            "expressible through the tag-side declarations)")
      .note(:dt_combination_ranking,
            "among multiple tag+signal combinations the one whose average of tag " \
            "and signal mBC is closest to the nominal D0 mass is retained; handled " \
            "internally by DTagAlg/DTagTool (the signal side is photon-based, i.e. " \
            "an ST + signal pattern rather than a two-tag DT)")
      .note(:background_veto,
            "pi0 veto for D0 -> gamma gamma: a signal photon combined with any " \
            "other photon must not form a pi0; events whose M(gamma gamma) falls " \
            "in the pi0 window are rejected (suppresses D0 -> pi0 pi0 and " \
            "continuum pi0 backgrounds)")
      .with_decay_card(decay_card_gg)

# Hadronic tag side: five D0 tag modes (charge-conjugate included; charm left
# free so both the D0 and the D0bar recoil are tagged)
alg_gg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPiPiPi0, :D0toKPiPi0Pi0
  t.window :mBC,    min: 1.847, max: 1.883   # 1.847 < mBC(tag) < 1.883 GeV/c^2
  t.window :deltaE, abs: 0.1                 # |DeltaE(tag)| < 0.1 GeV
end

# Signal side: the two most energetic showers not used by the tag
alg_gg.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0     # each photon > 10 deg from any charged track (Bhabha veto)
  s.min_photon_energy 0.025   # 25 MeV barrel energy floor
end

# Kinematic fit: 4C + D0 mass constraint on the signal gamma gamma + D0 mass
# constraint on the tag  ->  6 constraints
alg_gg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:D0)   # signal
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)            # tag
  f.chi2_cut 200
end

alg_gg.apply
alg_gg.execute_on([data_3773, incMC_3773, exMC_gg])

### ---------------------------------------------------------------------------
### Signal mode II : D0 -> pi0 pi0
### ---------------------------------------------------------------------------
alg_pi0pi0 = TagAnalysis.new("D0ToPi0Pi0")
alg_pi0pi0.set_header(["D0ToPi0Pi0Alg/D0ToPi0Pi0.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })
         .note(:tag_candidate_selection,
               "per tag mode one candidate is kept: the one whose reconstructed D0 " \
               "energy is closest to the beam energy (DTagAlg internal ranking)")
         .note(:dt_combination_ranking,
               "among multiple tag+signal combinations the one whose average of tag " \
               "and signal mBC is closest to the nominal D0 mass is retained " \
               "(internal to DTagAlg/DTagTool)")
         .note(:pi0_mass_constraints,
               "the two pi0 mass constraints on the signal pi0 pi0 cannot be " \
               "expressed in the tag-fit participant vocabulary (only :gamma with " \
               "the declared arity); they are applied as one additional 1-C " \
               "constraint per pi0 in the fit")
         .note(:photon_pair_selection,
               "the two pi0's are formed from the four most energetic signal showers; " \
               "the pairing with the smallest |DeltaE(pi0 pi0)| is selected and the " \
               "window -0.070 < DeltaE(pi0 pi0) < +0.075 GeV together with the " \
               "mBC(pi0 pi0) fit is applied afterwards at the ROOT stage")
         .with_decay_card(decay_card_pi0pi0)

# Same hadronic tag side as mode I
alg_pi0pi0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPiPiPi0, :D0toKPiPi0Pi0
  t.window :mBC,    min: 1.847, max: 1.883
  t.window :deltaE, abs: 0.1
end

# Signal side: the four most energetic showers not used by the tag
alg_pi0pi0.signal_side do |s|
  s.photons 4
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# Kinematic fit: 4C + D0 mass constraint on the signal pi0 pi0 + D0 mass
# constraint on the tag + two pi0 mass constraints
alg_pi0pi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma, :gamma, :gamma).constrain_to_nominal_mass_of(:D0)  # signal pi0 pi0 -> D0
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)                           # tag -> D0
  f.chi2_cut 200
end

alg_pi0pi0.apply
alg_pi0pi0.execute_on([data_3773, incMC_3773, exMC_pi0pi0])