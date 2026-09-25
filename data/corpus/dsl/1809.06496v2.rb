# 1809.06496v2: D0(+) → π-π0(+) e+ ν_e semileptonic decays and PWA at BESIII
# ψ(3770) at √s = 3.773 GeV with 2.93 fb⁻¹ at BESIII
# Double-tag (ST+missing) technique — TagAnalysis
# Two independent signal modes per Rule T1

# === Datasets ===
psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# =========================================================
# Mode I: D0 → π- π0 e+ ν_e (ρ- P-wave dominant)
# Tag: anti-D0 (charm -1) hadronic, Signal: D0 → π⁻π⁰ e⁺ ν_e
# =========================================================
decay_card_I = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0000 pi- pi0 e+ nu_e PHSP;
  Enddecay
  Decay anti-D0
  1.0000 K+ pi- PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_I = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_D0_pipi0_enu"
  config.related_dataset = psi3770_data
  config.events          = 200_000
  config.decay_card      = decay_card_I
  config.cross_section   = :default
end

alg_I = TagAnalysis.new("D0topipi0enu")
alg_I.set_header(["D0topipi0enuAlg/D0topipi0enu.h"])
     .set_constant({ "ECMS" => [:double, 3.773] })

alg_I.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPi0Pi0, :D0toKPiPiPi, :D0toKPiPiPiPi0
  t.charm -1
end

alg_I.signal_side do |s|
  s.photons 2
  s.charged(pim: 1, ep: 1)
  s.require_charge 0
  s.missing :nu_e
end

alg_I.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_I.with_decay_card(decay_card_I)
     .note(:double_tag_technique,
       "ST+missing pattern: hadronic tag (anti-D0) + semileptonic D0 → π⁻π⁰ e⁺ ν_e; " \
       "Umiss ≡ Emiss − |Pmiss| peaks at zero for signal neutrino")
     .note(:tag_modes,
       "Neutral D tags (anti-D0): K⁺π⁻, K⁺π⁻π⁰, K⁺π⁻π⁰π⁰, K⁺π⁻π⁻π⁺, K⁺π⁻π⁻π⁺π⁰; " \
       "tag yields from M_BC fits per Ref. [13]")
     .note(:pwa_fit,
       "Simultaneous unbinned PWA fit on both D0 and D+ modes in 5D phase space " \
       "(m_ππ, q², θ_e, θ_π, χ); D0 mode: ρ⁻ P-wave only; " \
       "FF ratios rV = V(0)/A1(0), r2 = A2(0)/A1(0) extracted")
     .note(:extra_photon_veto,
       "E_γ,max < 0.25 GeV for extra photons not used by tag or signal reconstruction; " \
       "applied in ROOT")
     .note(:bremsstrahlung_recovery,
       "e⁺ momentum improved by recovering FSR/bremsstrahlung energy in inner detector region")
     .note(:umiss_fit,
       "Signal yield from unbinned ML fit to Umiss distribution; " \
       "signal shape: signal MC convolved with Gaussian; " \
       "background shape: inclusive MC convolved with same Gaussian; " \
       "|Umiss| < 0.06 GeV applied for PWA sample selection in ROOT")
     .note(:pi0_selection,
       "Multiple π⁰ candidates resolved by choosing γγ pair with invariant mass closest to nominal π⁰ mass")
     .note(:bf_measurement,
       "B(D0 → ρ⁻ e⁺ ν_e) = (1.445 ± 0.058 ± 0.039) × 10⁻³; " \
       "systematic uncertainty 2.5% total")
     .apply

alg_I.execute_on([psi3770_data, psi3770_incMC, exMC_I])

# =========================================================
# Mode II: D+ → π- π+ e+ ν_e (ρ0 P-wave + ω + f0(500) S-wave)
# Tag: D- (charm -1) hadronic, Signal: D+ → π⁻π⁺ e⁺ ν_e
# =========================================================
decay_card_II = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D- PHSP;
  Enddecay
  Decay D+
  1.0000 pi- pi+ e+ nu_e PHSP;
  Enddecay
  Decay D-
  1.0000 K+ pi- pi- PHSP;
  Enddecay
  End
DECAYCARD

exMC_II = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Dp_pipienu"
  config.related_dataset = psi3770_data
  config.events          = 200_000
  config.decay_card      = decay_card_II
  config.cross_section   = :default
end

alg_II = TagAnalysis.new("Dptopipienu")
alg_II.set_header(["DptopipienuAlg/Dptopipienu.h"])
      .set_constant({ "ECMS" => [:double, 3.773] })

alg_II.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

alg_II.signal_side do |s|
  s.charged(pim: 1, pip: 1, ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

alg_II.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_II.with_decay_card(decay_card_II)
      .note(:double_tag_technique,
        "ST+missing pattern: hadronic tag (D⁻) + semileptonic D⁺ → π⁻π⁺ e⁺ ν_e; " \
        "Umiss ≡ Emiss − |Pmiss| peaks at zero for signal neutrino")
      .note(:tag_modes,
        "Charged D tags (D⁻): K⁺π⁻π⁻, K⁺π⁻π⁻π⁰, K_S⁰π⁻, K_S⁰π⁻π⁰, K_S⁰π⁻π⁺π⁻, K⁺K⁻π⁻; " \
        "tag yields from M_BC fits per Ref. [13]")
      .note(:pwa_fit,
        "Simultaneous unbinned PWA fit on both D0 and D+ modes in 5D phase space " \
        "(m_ππ, q², θ_e, θ_π, χ); D+ mode: ρ⁰ P-wave (GS) + ρ⁰−ω interference + " \
        "f0(500) S-wave; rV, r2, a_S, φ_S extracted from simultaneous fit")
      .note(:ks0_veto,
        "Veto π⁺π⁻ combinations with |M(π⁺π⁻) − M(K_S⁰)| < 70 MeV/c² " \
        "to suppress D⁺ → K_S⁰ e⁺ ν_e background (~98.3% rejection); applied in ROOT")
      .note(:extra_photon_veto,
        "E_γ,max < 0.25 GeV for extra photons not used by tag or signal reconstruction; " \
        "applied in ROOT")
      .note(:bremsstrahlung_recovery,
        "e⁺ momentum improved by recovering FSR/bremsstrahlung energy in inner detector region")
      .note(:umiss_fit,
        "Signal yield from unbinned ML fit to Umiss distribution; " \
        "signal shape: signal MC convolved with Gaussian; " \
        "background shape: inclusive MC convolved with same Gaussian; " \
        "|Umiss| < 0.06 GeV applied for PWA sample selection in ROOT")
      .note(:s_wave_observation,
        "First observation of S-wave (f0(500)) contribution in D⁺ → π⁻π⁺ e⁺ ν_e at >10σ; " \
        "fraction (25.7 ± 1.6 ± 1.1)%")
      .note(:ff_measurement,
        "Form factor ratios from simultaneous PWA: " \
        "rV = V(0)/A1(0) = 1.695 ± 0.083 ± 0.051, " \
        "r2 = A2(0)/A1(0) = 0.845 ± 0.056 ± 0.039")
      .note(:bf_measurement,
        "B(D⁺ → π⁻π⁺ e⁺ ν_e) = (2.449 ± 0.074 ± 0.073) × 10⁻³; " \
        "B(D⁺ → ρ⁰ e⁺ ν_e) = (1.860 ± 0.070 ± 0.061) × 10⁻³; " \
        "B(D⁺ → f0(500) e⁺ ν_e, f0(500)→π⁺π⁻) = (6.30 ± 0.43 ± 0.32) × 10⁻⁴; " \
        "systematic uncertainty 3.0% total")
      .apply

alg_II.execute_on([psi3770_data, psi3770_incMC, exMC_II])