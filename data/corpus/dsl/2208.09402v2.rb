# 2208.09402v2: Improved measurement of strong-phase difference δ_D^{Kπ}
# in quantum-correlated D⁰D⁰bar decays at ψ(3770)
# Tag-based double-tag analysis:
# Signal mode: D⁰ → K⁻π⁺
# Tags: CP-even eigenstates (K⁺K⁻, π⁺π⁻, π⁰π⁰, K_S⁰π⁰π⁰) + K_L⁰X modes
#       CP-odd eigenstates (K_S⁰π⁰, K_S⁰η, K_S⁰η', K_S⁰ω, K_S⁰φ) + K_L⁰π⁰π⁰
#       Quasi CP-even: π⁺π⁻π⁰
#       Mixed CP: K_S,L⁰π⁺π⁻ (binned in phase space)
# Uses DTagAlg pre-stored tag candidates

decay_card = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 VSS;
  Enddecay
  Decay D0
  1.000 K- pi+ PHSP;
  Enddecay
  Decay anti-D0
  1.000 K+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# CP-even tag modes available in DTagAlg vocabulary
# :D0toKK (K⁺K⁻), :D0toPiPi (π⁺π⁻), :D0toPi0Pi0 (π⁰π⁰), :D0toKsPi0Pi0 (K_S⁰π⁰π⁰)
# CP-odd tag modes available:
# :D0toKsPi0 (K_S⁰π⁰), :D0toKsEta (K_S⁰η), :D0toKsPiPiEta (K_S⁰η'→π⁺π⁻η~mode)

alg_cp_even = TagAnalysis.new("D0toKPi_CPevenTag")
alg_cp_even.set_header(["D0toKPiCPevenTag/D0toKPiCPevenTag.h"])
alg_cp_even.set_constant({ "ECMS" => [:double, 3.773] })
alg_cp_even.with_decay_card(decay_card)

# Tag side: D⁰ CP-even eigenstates
alg_cp_even.tag_side(:D0) do |t|
  t.modes :D0toKK, :D0toPiPi, :D0toPi0Pi0, :D0toKsPi0Pi0
  t.rank_by :inv
end

# Signal side: K⁻π⁺ (remaining tracks not used by tag)
alg_cp_even.signal_side do |s|
  s.photons 0
  s.charged(km: 1, pip: 1)
  s.require_charge 0
end

alg_cp_even.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_cp_even.note(:tag_mode_unavailable, "K_L⁰π⁰, K_L⁰ω (CP-even) and K_L⁰π⁰π⁰ (CP-odd) not in DTagAlg vocabulary; K_S⁰η' modes, K_S⁰ω, K_S⁰φ, π⁺π⁻π⁰, K_S,L⁰π⁺π⁻ tags handled in separate selections or ROOT analysis")
alg_cp_even.note(:cp_even_analysis, "CP-even tag yields for D⁰→K⁻π⁺ used to extract B(D_- → K⁻π⁺) for asymmetry A_Kπ measurement")
alg_cp_even.note(:kl_modes, "K_L⁰X modes require missing-mass technique M_miss² = (√s/2 - E_X)² - |p_T + p_X|², handled via partial_rec in separate algorithm")
alg_cp_even.note(:yield_extraction, "ST yields from M_BC fits, DT yields from M_BC (fully reco) or M_miss² (K_L modes) fits; asymmetry A_Kπ = (B(D_-→Kπ)-B(D_+→Kπ))/(B(D_-→Kπ)+B(D_+→Kπ)) extracted in ROOT")

alg_cp_even.apply
alg_cp_even.execute_on([DatasetManager.load_real_data.find("712_3773"), DatasetManager.load_inclusive_mc.find("712_3773")])

# CP-odd tags
alg_cp_odd = TagAnalysis.new("D0toKPi_CPoddTag")
alg_cp_odd.set_header(["D0toKPiCPoddTag/D0toKPiCPoddTag.h"])
alg_cp_odd.set_constant({ "ECMS" => [:double, 3.773] })
alg_cp_odd.with_decay_card(decay_card)

alg_cp_odd.tag_side(:D0) do |t|
  t.modes :D0toKsPi0, :D0toKsEta, :D0toKsPiPiEta
  t.rank_by :inv
end

alg_cp_odd.signal_side do |s|
  s.photons 0
  s.charged(km: 1, pip: 1)
  s.require_charge 0
end

alg_cp_odd.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_cp_odd.note(:cp_odd_analysis, "CP-odd tag yields for D⁰→K⁻π⁺ used to extract B(D_+ → K⁻π⁺) for asymmetry A_Kπ measurement; K_S⁰ω, K_S⁰φ, K_L⁰π⁰π⁰ tags processed separately")
alg_cp_odd.note(:ks_phi_and_omega, "K_S⁰φ (→K⁺K⁻) and K_S⁰ω (→π⁺π⁻π⁰) tag reconstructions not directly in DTagAlg vocabulary; handled via dedicated selection or ROOT")

alg_cp_odd.apply
alg_cp_odd.execute_on([DatasetManager.load_real_data.find("712_3773"), DatasetManager.load_inclusive_mc.find("712_3773")])