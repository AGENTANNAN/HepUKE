# DSL for paper 2403.14998v1: Precise measurement of e+e- → Ds+Ds- cross sections
# at √s from threshold to 4.95 GeV using a single tag method
# Tag-based analysis: ST + missing Ds

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# Decay card for signal MC: e+e- → Ds+Ds- via KKMC
# Top mother psi(4260) per BESIII convention for KKMC generator;
# ISR and beam energy spread handled by KKMC
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s+ D_s- PHSP;
    Enddecay

    Decay D_s+
    1.000 K+ K- pi+ PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Query real data and inclusive MC in 3.94-4.95 GeV range
# 138 energy points (XYZ + R-scan), BOSS 703/705/706/707, 22.9 fb-1 total
scan_data = DatasetManager.real_data.where(cms_energy: {value: 3940..4950})
scan_incMC = DatasetManager.inclusive_mc.where(cms_energy: {value: 3940..4950})

# Exclusive signal MC for each energy scan point
# Flat Born cross section as initial input; iteratively updated via MC-weighting
scan_exMC = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name = "sig_DsDs"
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

# TagAnalysis: single tag Ds- → K+K-π-, Ds+ identified by recoil mass
alg = TagAnalysis.new("DsDsTag")
alg.set_header(["DsDsTagAlg/DsDsTag.h"])
   .set_constant({ "ECMS" => [:double, 4.400] })
   .note(:submode_selection, "Ds- candidates selected via two intermediate decay modes: Ds- → φπ- with φ → K+K- (1.005 < M(K+K-) < 1.035 GeV/c²) or Ds- → K*(892)0 K- with K*(892)0 → K+π- (0.832 < M(K+π-) < 0.928 GeV/c², helicity angle |cos θ(K+)| < 0.52)")
   .note(:recoil_mass_window, "Recoil mass window 1.945 < RM(Ds-) < 1.990 GeV/c² applied to suppress background; recoil mass resolution improved using mass-constrained formula: M_recoil(K+K-π-) + M(K+K-π-) - m(Ds-)")
   .note(:peaking_background_subtraction, "Peaking background from e+e- → Ds+Ds*- at Ecms > 4.6 GeV estimated from exclusive MC, normalized to luminosity and cross sections, and subtracted from signal yields. e+e- → Ds*+Ds*- fails the recoil mass requirement.")
   .note(:iterative_isr_correction, "MC-weighting method with 5 iterations: flat Born cross section as initial input; smoothed lineshape from variable span smoother updates ISR correction factors (1+δ), efficiencies, and measured cross sections iteratively until convergence")
   .note(:cross_section_formula, "Born cross section σ_Born = (N_fit - N_bkg) / (2 × B(Ds→KKπ) × ε × (1+δ) × (1/|1-Π|²) × L) where factor 1/2 accounts for both Ds+ and Ds- single tag reconstruction")
   .with_decay_card(decay_card)

alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

alg.signal_side do |s|
  s.missing :Ds
end

alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.apply
alg.execute_on(scan_data + scan_incMC + scan_exMC)