# =====================================================================
# BOSS / DSL specification for the ψ(3770) → D0/D+ → ω/η + pions analysis
# using the double-tag technique.  Tag-based analysis → the TagAnalysis
# surface (no Selection, no select_track/select_photon/pid): the tag side
# carries its own selected, PID'd tracks/showers; the signal side is built
# from what the tag did not use.
# =====================================================================

### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")      # ψ(3770) real data
incmc_3773 = DatasetManager.inclusive_mc.find("712_3773")   # matching inclusive MC

# Decay card for the 200k-event exclusive signal MC:
# ψ(3770) → D0 anti-D0, D0 → ω π+π−, anti-D0 → K− π+, ω → π+π−π0, π0 → γγ.
# The description states this single sample is shared by both the D0 and D+
# analyses, so the same card is attached to both TagAnalysis objects.
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 omega pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K- pi+ PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0barD0_omega_pipi"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# =====================================================================
# Analysis 1 — anti-D0 tag + D0 → ω/η π+π− and D0 → ω/η π0π0
# =====================================================================
alg_d0 = TagAnalysis.new("D0TagOmegaEta")
alg_d0.set_header(["D0TagOmegaEtaAlg/D0TagOmegaEta.h"])
      .set_constant({ "ECMS" => [:double, 3.773] })
      .with_decay_card(decay_card_signal)

# Tag side: anti-D0 reconstructed in hadronic modes (K+π−, K+π−π0, K+π−π−π+).
# charm −1 pins the anti-D0 (c̄) side opposite the charm-+1 D0 signal.
alg_d0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: remaining content of the D0 decay.
# D0 modes → at least one π+ and one π− (net charge 0) with 2–6 photons
# (γγ from π0 → γγ, up to the ω/η π0π0 channel).
alg_d0.signal_side do |s|
  s.charged(pip: 1, pim: 1, at_least: true)
  s.require_charge(0)
  s.photons 2..6
end

# Kinematic fit: 4C conservation, γγ constrained to the π0 nominal mass,
# loose χ² < 200 (tight/mode-dependent cut applied later in ROOT).
alg_d0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)          # π0 mass window
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:pip, :pim).out_of(0.475, 0.520)               # K_S0 → π+π− veto
  f.chi2_cut 200
end

# BOSS-side steps that are not expressible as DSL constructs.
alg_d0
  .note(:background_veto,
        "Tag-side background vetoes: K_S0 veto |M(pi+pi-) - M(K_S0)| > 30 MeV/c^2 " \
        "for the four-track anti-D0 tag mode (K+pi-pi-pi+).")
  .note(:signal_ks0_veto,
        "Signal-side K_S0 rejection in the pi0pi0 channel: 0.448 < M(pi0pi0) < 0.548 GeV/c^2.")
  .note(:signal_resonance_selection,
        "ω/η reconstruction from π+π−π0 combinations: M(π+π−π0) < 0.9 GeV/c^2; " \
        "ω signal window 0.74–0.82 GeV/c^2 and η signal window 0.52–0.57 GeV/c^2 " \
        "(resonance windows applied on the reconstructed ω/η candidates); " \
        "π0 candidates require at least one photon in the EMC barrel.")
  .note(:tag_deltaE_cut,
        "Mode-dependent tag-side ΔE windows: 3.0σ for D0 → ω π+π− and 3.5σ for " \
        "D0 → ω π0π0; the candidate with minimum |ΔE| is retained. These σ-based " \
        "windows are resolution-dependent and are kept out of the BOSS-level DSL " \
        "(tag ΔE is stored and windowed in ROOT).")

alg_d0.apply           # takes NO Selection argument for a TagAnalysis
alg_d0.execute_on([data_3773, incmc_3773, exMC_signal])

# =====================================================================
# Analysis 2 — D− tag + D+ → ω/η π+π0
# =====================================================================
alg_dp = TagAnalysis.new("DpTagOmegaEta")
alg_dp.set_header(["DpTagOmegaEtaAlg/DpTagOmegaEta.h"])
      .set_constant({ "ECMS" => [:double, 3.773] })
      .with_decay_card(decay_card_signal)

# Tag side: D− (charm −1) reconstructed in the listed hadronic modes.
alg_dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0,
          :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi,
          :DptoKKPi
  t.charm -1
end

# Signal side: D+ → ω/η π+π0 → at least two π+ and one π− (net charge +1)
# with at least four photons (π0 → γγ from both the D+ and the ω/η π0).
alg_dp.signal_side do |s|
  s.charged(pip: 2, pim: 1, at_least: true)
  s.require_charge(1)
  s.photons 4
end

# Kinematic fit: 4C conservation with the γγ invariant mass constrained to
# the π0 nominal mass (each π0 constrained separately — see note), χ² < 200.
alg_dp.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)          # π0 mass window
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:pip, :pim).out_of(0.475, 0.520)               # K_S0 → π+π− veto
  f.chi2_cut 200
end

alg_dp
  .note(:background_veto,
        "Tag-side background vetoes for the D− modes: Λ (M(pbar π+) ∈ [1.110,1.120]), " \
        "K_S0 (M(π+π−) ∈ [0.480,0.520]) and Σ− (M(pbar π0) ∈ [1.170,1.200]) GeV/c^2.")
  .note(:signal_resonance_selection,
        "ω/η reconstruction from π+π−π0 combinations: M(π+π−π0) < 0.9 GeV/c^2; " \
        "ω window 0.74–0.82 GeV/c^2 and η window 0.52–0.57 GeV/c^2; π0 candidates " \
        "require at least one photon in the EMC barrel.")
  .note(:tag_deltaE_cut,
        "Mode-dependent tag-side ΔE window at 3.5σ, choosing the minimum |ΔE| " \
        "combination; the σ-based window is stored and cut in ROOT rather than " \
        "imposed at BOSS level.")
  .note(:pi0_multiplicity_in_fit,
        "Each π0 of the D+ final state must be mass-constrained separately in the " \
        "kinematic fit; the DSL TagFit declares one γγ → π0 constraint, the " \
        "per-π0 iteration being carried out in the generated fit step.")

alg_dp.apply
alg_dp.execute_on([data_3773, incmc_3773, exMC_signal])