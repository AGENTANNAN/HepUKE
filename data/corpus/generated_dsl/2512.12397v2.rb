### Dataset preparation ###
# psi(3770) real data and the matching tag-based inclusive MC (712_3773, ~20.3 fb^-1)
data_3770  = DatasetManager.real_data.find("712_3773")
incMC_3770 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card: psi(3770) -> D+ D-, signal D+ -> pi+ pi0 pi0,
# tag D- decays into the six hadronic tag modes used by the analysis.
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-                        PHSP;
    Enddecay

    Decay D+
    1.0000 pi+ pi0 pi0                  PHSP;
    Enddecay

    Decay D-
    0.1667 K+ pi- pi-                   PHSP;
    0.1667 K+ pi- pi- pi0               PHSP;
    0.1666 K_S0 pi-                     PHSP;
    0.1667 K_S0 pi- pi0                 PHSP;
    0.1666 K_S0 pi- pi- pi+             PHSP;
    0.1667 K+ K- pi-                    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                  PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                      PHSP;
    Enddecay

    End
DECAYCARD

# 2M-event exclusive signal MC of psi(3770) -> D+ D- (signal + tag decays included)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_DpToPipPi0Pi0_tag"
  config.related_dataset = data_3770
  config.events          = 2_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (tag-based: D- tag + signal D+ -> pi+ pi0 pi0) ###
alg_name = "DpTagPipPi0Pi0"
tag_alg  = TagAnalysis.new(alg_name)
tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.773]})
       .with_decay_card(decay_card_signal)

# Tag side: D- reconstructed in six hadronic modes
tag_alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1                                             # pin the tagged side to D-
  t.window :mBC, min: 1.865, max: 1.875                  # tag-side M_BC requirement
end

# Signal side: everything the tag did not use -> one pi+ (net charge +1) and good showers
tag_alg.signal_side do |s|
  s.charged(pip: 1)                                      # exactly one pi+
  s.photons 4..48                                        # 4-48 signal-side photons (event-level cap 48)
  s.require_charge 1                                     # net charge +1
  s.min_photon_angle 10.0                                # lab angle vs beam > 10 degrees
  s.min_photon_energy 0.025                              # photon energy > 25 MeV
end

# Kinematic fit: 4-momentum constraint plus pi0 mass constraints from gamma-gamma pairs
tag_alg.fit do |f|
  f.constrain_four_momentum                                                       # 4C
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)                       # pi0 mass window
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)          # 1C pi0 mass constraint
  f.chi2_cut 200                                                                  # loose chi2 < 200
end

tag_alg
  .note(:pi0_reconstruction, "Each pi0 formed from a gamma-gamma pair with 0.115 < M(gamma gamma) < 0.150 GeV/c^2, requiring at least one photon in the barrel EMC; 1C mass-constrained pi0 fit with chi2 < 50.")
  .note(:background_veto_ks_pi, "Veto D+ -> K_S0 pi+ by rejecting 0.428 < M(pi0 pi0) < 0.548 GeV/c^2 (K_S0 -> pi0 pi0).")
  .note(:background_veto_mispartition, "Reject D+ -> pi+ pi0 pi0 / D- -> K+ K- pi- mispartition when both |M(K- pi+ pi0) - 1.865| < 0.05 and |M(K+ pi- pi0) - 1.865| < 0.05 GeV/c^2.")
  .note(:deltaE_window_per_mode, "Per-mode DeltaE windows, applied in ROOT (store-not-cut): D+ -> pi+ pi0 pi0 (-0.100, 0.045); D- -> K+ pi- pi- (-0.025, 0.024); D- -> K+ pi- pi- pi0 (-0.057, 0.046); D- -> K_S0 pi- (-0.025, 0.026); D- -> K_S0 pi- pi0 (-0.062, 0.049); D- -> K_S0 pi- pi- pi+ (-0.028, 0.027); D- -> K+ K- pi- (-0.024, 0.023) GeV.")
  .note(:signal_mbc_window, "Signal-side M_BC required in 1.865-1.875 GeV/c^2 (tag-side window declared above).")
  .note(:best_candidate, "Keep the single best candidate minimizing the sum of squared DeltaE on the tag and signal sides.")
  .note(:seven_c_constraint, "Amplitude analysis uses a 7C fit (initial e+e- four-momentum, D+ mass, and both pi0 masses); applied after the BOSS selection.")
  .apply

# Execute on real data, inclusive MC and the exclusive signal MC
root_files = tag_alg.execute_on([data_3770, incMC_3770, exMC_signal])