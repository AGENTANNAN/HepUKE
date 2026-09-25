### Dataset description ###
# ψ(3770) data and corresponding inclusive MC at 3.773 GeV (BOSS 7.1.2 sample name 712_3773)
data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

### Decay cards (EvtGen format) ###
# Signal mode I: ψ(3770) → D0(→ K_S π+π−) anti-D0(→ K+ π−), with K_S → π+π−
decay_card_modeI = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_S0 pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Signal mode II: ψ(3770) → D0(→ K_S K+ K−) anti-D0(→ K+ π−), with K_S → π+π−
decay_card_modeII = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_S0 K+ K- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC (500,000 events per signal mode) ###
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKsPiPi_DTag"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKsKK_DTag"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection (tag-based / double-tag) ###
# DTagAlg channel names for the tag anti-D0 (charm −1)
tag_modes = [
  :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi,                      # flavor tags  (Kπ, Kππ0, Kπππ)
  :D0toKK, :D0toPiPi, :D0toKsPi0Pi0, :D0toPiPiPi0,          # CP-even tags (KK, ππ, Ksπ0π0, πππ0)
  :D0toKsPi0, :D0toKsEta, :D0toKsEtaPrime, :D0toKsOmega,    # CP-odd tags  (Ksπ0, Ksη, Ksη′, Ksω)
  :D0toKsPiPi, :D0toKsKK                                    # self-conjugate K_S h+h− modes
]

# ---------- Mode I: signal D0 → K_S π+π− ----------
alg_modeI = TagAnalysis.new("D0ToKsPiPiDTag")
alg_modeI.set_header(["D0ToKsPiPiDTagAlg/D0ToKsPiPiDTag.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })
         .with_decay_card(decay_card_modeI)

# Signal side of the double tag: D0 (charm +1) → K_S π+π−
alg_modeI.tag_side(:D0) do |t|
  t.modes :D0toKsPiPi
  t.charm 1
  t.rank_by :inv
end

# Tag side of the double tag: anti-D0 (charm −1), flavor / CP / self-conjugate modes
alg_modeI.tag_side(:D0) do |t|
  t.modes(*tag_modes)
  t.charm -1
  t.rank_by :inv
end

# Four-momentum-constrained kinematic fit; D masses fixed to their PDG values
alg_modeI.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
end

# Inexpressible BOSS-side procedures
alg_modeI
  .note(:missing_kl_treatment, "the missing K_L is treated via the missing-mass-squared M_miss^2 or the U_miss variable, with the K_L mass fixed to its PDG value when forming these observables")
  .note(:efficiency_curve, "the Dalitz-plot efficiency for D0 -> K_S pi+pi- is taken from phase-space exclusive MC (efficiency is not flat over the Dalitz plane)")
  .note(:yield_extraction, "signal yields are extracted from M_BC fits; tag mBC/deltaE are stored and windowed downstream, not cut in BOSS")
  .note(:fourier_weights, "per-event optimal Fourier weights are applied to the signal")
  .note(:background_veto, "peaking backgrounds for the K_S K+ K- channel are handled at the fit stage")

alg_modeI.apply
alg_modeI.execute_on([data_3773, incMC_3773, exMC_modeI])

# ---------- Mode II: signal D0 → K_S K+K− ----------
alg_modeII = TagAnalysis.new("D0ToKsKKDTag")
alg_modeII.set_header(["D0ToKsKKDTagAlg/D0ToKsKKDTag.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
          .with_decay_card(decay_card_modeII)

# Signal side of the double tag: D0 (charm +1) → K_S K+K−
alg_modeII.tag_side(:D0) do |t|
  t.modes :D0toKsKK
  t.charm 1
  t.rank_by :inv
end

# Tag side of the double tag: anti-D0 (charm −1), flavor / CP / self-conjugate modes
alg_modeII.tag_side(:D0) do |t|
  t.modes(*tag_modes)
  t.charm -1
  t.rank_by :inv
end

# Four-momentum-constrained kinematic fit; D masses fixed to their PDG values
alg_modeII.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
end

# Inexpressible BOSS-side procedures
alg_modeII
  .note(:missing_kl_treatment, "the missing K_L is treated via the missing-mass-squared M_miss^2 or the U_miss variable, with the K_L mass fixed to its PDG value when forming these observables")
  .note(:yield_extraction, "signal yields are extracted from M_BC fits; tag mBC/deltaE are stored and windowed downstream, not cut in BOSS")
  .note(:fourier_weights, "per-event optimal Fourier weights are applied to the signal")
  .note(:background_veto, "peaking backgrounds for the K_S K+ K- channel are handled at the fit stage")

alg_modeII.apply
alg_modeII.execute_on([data_3773, incMC_3773, exMC_modeII])