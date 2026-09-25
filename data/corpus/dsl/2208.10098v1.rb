# 2208.10098v1: Measurement of CP-even fraction of D⁰ → π⁺π⁻π⁺π⁻
# in quantum-correlated D⁰D⁰bar decays at ψ(3770)
# Tag-based double-tag analysis:
# Signal mode: D⁰ → π⁺π⁻π⁺π⁻
# Tags: CP-even eigenstates (K⁺K⁻, K_S⁰π⁰π⁰, K_L⁰π⁰, K_L⁰ω)
#       CP-odd eigenstates (K_S⁰π⁰, K_S⁰η, K_S⁰η', K_S⁰ω, K_S⁰φ, K_L⁰π⁰π⁰)
#       quasi CP-even: π⁺π⁻π⁰
#       mixed CP: K_S,L⁰π⁺π⁻ (binned in phase space)
# Uses DTagAlg pre-stored tag candidates

decay_card = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 VSS;
  Enddecay
  Decay D0
  1.000 pi+ pi- pi+ pi- PHSP;
  Enddecay
  Decay anti-D0
  1.000 pi+ pi- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# CP-even tags available in DTagAlg vocabulary:
# :D0toKK (K⁺K⁻), :D0toKsPi0Pi0 (K_S⁰π⁰π⁰)
# CP-odd tags available:
# :D0toKsPi0 (K_S⁰π⁰), :D0toKsEta (K_S⁰η), :D0toKsPiPiEta (K_S⁰η')

# CP-even tag analysis
alg_cp_even = TagAnalysis.new("D0to4Pi_CPevenTag")
alg_cp_even.set_header(["D0to4PiCPevenTag/D0to4PiCPevenTag.h"])
alg_cp_even.set_constant({ "ECMS" => [:double, 3.773] })
alg_cp_even.with_decay_card(decay_card)

alg_cp_even.tag_side(:D0) do |t|
  t.modes :D0toKK, :D0toKsPi0Pi0
  t.rank_by :inv
end

# Signal side: π⁺π⁻π⁺π⁻ (4 charged pions)
alg_cp_even.signal_side do |s|
  s.photons 0
  s.charged(pip: 2, pim: 2)
  s.require_charge 0
end

alg_cp_even.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_cp_even.note(:tag_mode_unavailable, "CP-even tags: K_L⁰π⁰, K_L⁰ω, π⁰π⁰ not in DTagAlg vocabulary; CP-odd tags: K_S⁰ω, K_S⁰φ, K_L⁰π⁰π⁰ not in DTagAlg vocabulary")
alg_cp_even.note(:kl_modes, "K_L⁰X modes require missing-mass technique for DT selection; quasi-CP-even π⁺π⁻π⁰ tag and K_S,L⁰π⁺π⁻ binned analysis handled in separate algorithms or ROOT")
alg_cp_even.note(:ks_veto, "K_S⁰ veto applied to signal π⁺π⁻π⁺π⁻: events rejected if any π⁺π⁻ pair has invariant mass in [0.481, 0.514] GeV/c² and L/σ_L > 2")
alg_cp_even.note(:cp_fraction, "CP-even fraction F_+^{4π} determined from combination of CP-eigenstate tags, D→π⁺π⁻π⁰ tag, and D→K_S,L⁰π⁺π⁻ binned analysis; final fit in ROOT")
alg_cp_even.note(:deltae_signal, "DT signal-side ΔE required within [-0.026, 0.023] GeV")
alg_cp_even.note(:four_track_requirement, "DT selection: exactly 4 extra charged tracks, all identified as pions, net charge zero")

alg_cp_even.apply
alg_cp_even.execute_on([DatasetManager.load_real_data.find("712_3773"), DatasetManager.load_inclusive_mc.find("712_3773")])

# CP-odd tag analysis
alg_cp_odd = TagAnalysis.new("D0to4Pi_CPoddTag")
alg_cp_odd.set_header(["D0to4PiCPoddTag/D0to4PiCPoddTag.h"])
alg_cp_odd.set_constant({ "ECMS" => [:double, 3.773] })
alg_cp_odd.with_decay_card(decay_card)

alg_cp_odd.tag_side(:D0) do |t|
  t.modes :D0toKsPi0, :D0toKsEta, :D0toKsPiPiEta
  t.rank_by :inv
end

alg_cp_odd.signal_side do |s|
  s.photons 0
  s.charged(pip: 2, pim: 2)
  s.require_charge 0
end

alg_cp_odd.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_cp_odd.note(:yield_determination, "ST yields from M_BC fits with ARGUS background; DT yields from M_BC fits (fully reco) or M_miss² fits (K_L⁰X); peaking backgrounds estimated from MC with quantum-correlation corrections")

alg_cp_odd.apply
alg_cp_odd.execute_on([DatasetManager.load_real_data.find("712_3773"), DatasetManager.load_inclusive_mc.find("712_3773")])