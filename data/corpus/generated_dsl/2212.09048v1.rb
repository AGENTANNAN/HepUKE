# BOSS-stage DSL for the ψ(3770) → D⁰ D̄⁰ double-tag analysis at 3.773 GeV.
# Two independent signal modes (different final states / selection):
#   Mode 1 : D⁰ → K_L⁰ π⁺π⁻   (K_L⁰ missing)          tag D̄⁰ → K⁺π⁻, K⁺π⁻π⁺π⁻, K⁺π⁻π⁰
#   Mode 2 : D⁰ → K_S⁰ π⁺π⁻   (K_S⁰ → π⁺π⁻ reco.)      tag D̄⁰ → K⁺π⁻, K⁺π⁻π⁺π⁻, K⁺π⁻π⁰
# Tag-based analyses use TagAnalysis (no Selection, no select_track/select_photon/pid).

### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")       # ψ(3770) real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")    # matching inclusive MC

# Signal decay card — Mode 1: ψ(3770) → D⁰ D̄⁰, D⁰ → K_L⁰ π⁺π⁻, D̄⁰ → K⁺π⁻
decay_card_kl = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_L0 pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Signal decay card — Mode 2: ψ(3770) → D⁰ D̄⁰, D⁰ → K_S⁰ π⁺π⁻, K_S⁰ → π⁺π⁻, D̄⁰ → K⁺π⁻
decay_card_ks = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_S0 pi+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC — 1,000,000 events per signal mode
exMC_kl = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_D0toKLpipi_DT"
  config.related_dataset = data_3773
  config.events = 1_000_000
  config.decay_card = decay_card_kl
  config.cross_section = :default
end

exMC_ks = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_D0toKSpipi_DT"
  config.related_dataset = data_3773
  config.events = 1_000_000
  config.decay_card = decay_card_ks
  config.cross_section = :default
end

### Mode 1 — tag D̄⁰ (flavor tag) + signal D⁰ → K_L⁰ π⁺π⁻ with K_L⁰ missing ###
alg_kl = TagAnalysis.new("D0TagKLpipi")
alg_kl.set_header(["D0TagKLpipiAlg/D0TagKLpipi.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .note(:background_veto, "π⁰ veto: reject any γγ pair with M(γγ) in [0.095, 0.165] GeV/c²; " \
                              "η veto: reject any γγ pair with M(γγ) in [0.48, 0.58] GeV/c²")
      .note(:tag_side_windows, "tag-side ΔE and M_BC windows are ±3σ (resolution-based, not a fixed " \
                               "BOSS cut) -> stored unconditionally and applied in ROOT; the tag " \
                               "candidate with the smallest |ΔE| is selected in ROOT")
      .with_decay_card(decay_card_kl)

# Tag side: hadronic flavor tags of the D̄⁰
alg_kl.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPiPi, :D0toKPiPi0   # K⁺π⁻, K⁺π⁻π⁺π⁻, K⁺π⁻π⁰
  t.charm -1                                    # pinned to the D̄⁰ tag
end

# Signal side: one π⁺ and one π⁻ (net charge 0); K_L⁰ left missing
alg_kl.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.require_charge 0
  s.missing :K_L0                               # K_L⁰ unreconstructed -> missing-mass observables auto-stored
end

# Fit: total four-momentum conservation (tag + signal + missing), χ² < 200
alg_kl.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_kl.apply
alg_kl.execute_on([data_3773, incMC_3773, exMC_kl])

### Mode 2 — tag D̄⁰ (flavor tag) + signal D⁰ → K_S⁰ π⁺π⁻ with K_S⁰ → π⁺π⁻ fully reconstructed ###
alg_ks = TagAnalysis.new("D0TagKSpipi")
alg_ks.set_header(["D0TagKSpipiAlg/D0TagKSpipi.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .note(:ks0_secondary_vertex_fit, "K_S⁰ reconstructed from π⁺π⁻ via a secondary-vertex fit, with " \
                                       "mass window [0.485, 0.510] GeV/c² and flight significance L/σ_L > 2; " \
                                       "the resulting K_S⁰ candidate is combined with the two remaining " \
                                       "signal pions into the D⁰")
      .note(:background_veto, "π⁰ veto: reject any γγ pair with M(γγ) in [0.095, 0.165] GeV/c²; " \
                              "η veto: reject any γγ pair with M(γγ) in [0.48, 0.58] GeV/c²")
      .note(:tag_side_windows, "tag-side ΔE and M_BC windows are ±3σ (resolution-based, not a fixed " \
                               "BOSS cut) -> stored unconditionally and applied in ROOT; the tag " \
                               "candidate with the smallest |ΔE| is selected in ROOT")
      .with_decay_card(decay_card_ks)

# Tag side: hadronic flavor tags of the D̄⁰
alg_ks.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPiPi, :D0toKPiPi0   # K⁺π⁻, K⁺π⁻π⁺π⁻, K⁺π⁻π⁰
  t.charm -1                                    # pinned to the D̄⁰ tag
end

# Signal side: exactly four charged tracks (K_S⁰ → π⁺π⁻ and D⁰ → K_S⁰ π⁺π⁻) and zero photons
alg_ks.signal_side do |s|
  s.charged(pip: 2, pim: 2)   # exactly 4 charged tracks; no photons declared -> zero photons
end

# 6C fit: 4C four-momentum conservation + nominal D⁰ mass + nominal K_S⁰ mass, χ² < 200
alg_ks.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut 200
end

alg_ks.apply
alg_ks.execute_on([data_3773, incMC_3773, exMC_ks])