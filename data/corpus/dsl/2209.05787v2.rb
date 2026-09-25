# 2209.05787v2: Search for baryon and lepton number violating decays
# D± → n(anti-n) e± at ψ(3770) (√s = 3.773 GeV)
# Tag-based double-tag analysis:
# ST: D± in 6 hadronic modes
# DT: D± → n(anti-n) e± with Δ|B-L|=0 or Δ|B-L|=2
# Signal side: electron + missing (anti-)neutron, 2C kinematic fit
# Uses DTagAlg pre-stored tag candidates

decay_card = <<~DECAYCARD
  Decay psi(3770)
  1.000 D+ D- VSS;
  Enddecay
  # Signal decays
  Decay D+
  1.000 anti-n e+ PHSP;
  Enddecay
  Decay D-
  1.000 n e- PHSP;
  Enddecay
  End
DECAYCARD

dataset_3773 = DatasetManager.load_real_data.find("712_3773")
inc_mc_3773 = DatasetManager.load_inclusive_mc.find("712_3773")

# ============================================================
# D⁺ → anti-n e⁺ (Δ|B-L| = 0) and D⁺ → n e⁺ (Δ|B-L| = 2)
# Single tag side: D⁺ reconstructed in 6 hadronic modes
# ============================================================

# All 6 ST D⁺ modes are available in the DTagAlg vocabulary:
# :DptoKPiPi (K⁻π⁺π⁺), :DptoKPiPiPi0 (K⁻π⁺π⁺π⁰),
# :DptoKsPi (K_S⁰π⁺), :DptoKsPiPi0 (K_S⁰π⁺π⁰),
# :DptoKsPiPiPi (K_S⁰π⁺π⁺π⁻), :DptoKKPi (K⁺K⁻π⁺)

alg_dplus = TagAnalysis.new("DplustoNeutronE")
alg_dplus.set_header(["DplustoNeutronE/DplustoNeutronE.h"])
alg_dplus.set_constant({ "ECMS" => [:double, 3.773] })
alg_dplus.with_decay_card(decay_card)

# Tag side: D⁺ in 6 hadronic ST modes (charm +1)
alg_dplus.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm 1
  t.rank_by :inv
end

# Signal side: electron + missing neutron (massive, unknown mass)
# The neutron (anti-neutron) is treated as a missing particle
# with unknown mass in the 2C kinematic fit
alg_dplus.signal_side do |s|
  s.photons 0
  s.charged(ep: 1)
  s.missing :X0, mass: nil   # neutron/anti-neutron with unknown mass
  s.require_charge 1
end

# 2C kinematic fit: energy-momentum conservation +
# constrain ST D⁺ mass + constrain (e⁺ + missing) to D⁻ mass
alg_dplus.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Dplus)
  f.invariant_mass_of(:ep, :X0).constrain_to_nominal_mass_of(:"D-")
  f.chi2_cut 9999   # no chi2 cut applied; fit only required to converge
end

alg_dplus.note(:dt_selection, "DT: electron candidate + no additional charged tracks; ST ΔE within (-25,+25) MeV (or (-55,+40) MeV for modes with π⁰); ST M_BC ∈ (1.863, 1.877) GeV/c² from ARGUS fit")
alg_dplus.note(:neutron_emc_matching, "Post-fit: EMC shower within 30° of fitted neutron direction; further 10°(15°) opening angle cut; GBDT MVA on shower shape variables (E_tot, N_hit, A20, A42 Zernike moments) in neutron momentum bins")
alg_dplus.note(:st_yield, "ST yields from binned maximum-likelihood fits to M_BC with MC-convolved double-Gaussian signal + ARGUS background; N_ST^tot(D⁺) = 758.2k")
alg_dplus.note(:dt_fit, "DT signal yield from unbinned ML fit to neutron mass M_n/anti-n from kinematic fit; signal and background shapes from MC; no significant signal observed")
alg_dplus.note(:upper_limits, "Upper limits at 90% CL set via likelihood integration with systematic uncertainties convolved; B(D⁺→anti-n e⁺) < 1.43×10⁻⁵ (Δ|B-L|=0), B(D⁺→n e⁺) < 2.91×10⁻⁵ (Δ|B-L|=2)")

alg_dplus.apply
alg_dplus.execute_on([dataset_3773, inc_mc_3773])

# ============================================================
# D⁻ → n e⁻ (Δ|B-L| = 0) and D⁻ → anti-n e⁻ (Δ|B-L| = 2)
# Same ST modes for D⁻ (charge conjugate)
# ============================================================

alg_dminus = TagAnalysis.new("DminustoNeutronE")
alg_dminus.set_header(["DminustoNeutronE/DminustoNeutronE.h"])
alg_dminus.set_constant({ "ECMS" => [:double, 3.773] })
alg_dminus.with_decay_card(decay_card)

alg_dminus.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.rank_by :inv
end

alg_dminus.signal_side do |s|
  s.photons 0
  s.charged(em: 1)
  s.missing :X0, mass: nil   # neutron/anti-neutron with unknown mass
  s.require_charge -1
end

alg_dminus.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Dplus)
  f.invariant_mass_of(:em, :X0).constrain_to_nominal_mass_of(:"D+")
  f.chi2_cut 9999
end

alg_dminus.note(:st_yield, "ST yields from M_BC fits; N_ST^tot(D⁻) = 763.9k; same selection criteria as D⁺")
alg_dminus.note(:dt_fit, "Same DT selection as D⁺; no signal observed; upper limits set")

alg_dminus.apply
alg_dminus.execute_on([dataset_3773, inc_mc_3773])