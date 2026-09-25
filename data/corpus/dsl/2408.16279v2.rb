# D0 → π+π-π+π- strong-phase difference measurement
# Tag-based double-tag analysis at ψ(3770), 2.93 fb⁻¹
# arXiv:2408.16279v2

psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

decay_card = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay
  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "d0_4pi_strongphase"
  config.related_dataset = psi3770_data
  config.events = 200_000
  config.decay_card = decay_card
  config.cross_section = :default
end

alg = TagAnalysis.new("D0To4PiStrongPhase")
alg.set_header(["D0To4PiStrongPhaseAlg/D0To4PiStrongPhase.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })
   .with_decay_card(decay_card)

# Tag side 1: D0 hadronic flavor-tag modes (D⁰ → K⁻π⁺, K⁻π⁺π⁰, K⁻π⁺π⁻π⁺)
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1
end

# Signal side: D⁰ → π⁺π⁻π⁺π⁻ from remaining tracks
alg.signal_side do |s|
  s.charged(pip: 2, pim: 2)
  s.require_charge 0
end

# 4C kinematic fit: tag + signal = ecms_lab with D⁰ mass constraint on tag
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# Tag mBC and ΔE are stored unconditionally; windowed in ROOT.
# Semileptonic K⁻e⁺ν_e and K_L⁰ modes use partial reconstruction (ROOT-level UMmiss/Mmiss² fits).
alg.note(:partial_reco_tags, "semileptonic D→K⁻e⁺ν_e and K_L⁰ modes use partial reconstruction with U_miss/M²_miss fits in ROOT; these tag modes are not available in DTagAlg")

# DCS correction factors applied to hadronic flavor-tag yields in ROOT
alg.note(:dcs_correction, "doubly-Cabibbo-suppressed corrections from r_D^F, δ_D^F, R_F for K⁻π⁺π⁰ and K⁻π⁺π⁻π⁺ tags applied in ROOT")

# Phase-space binning determined from amplitude model; bin migration efficiency matrix from MC
alg.note(:phase_space_binning, "two binning schemes (equal-Δδ_D, optimal); 2×5 bins with KS⁰ veto [0.481,0.514] GeV; efficiency+migration matrix from MC")

alg.apply
alg.execute_on([psi3770_data, psi3770_incMC, exMC_signal])