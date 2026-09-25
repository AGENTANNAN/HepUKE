# =====================================================================
# psi(3770) @ sqrt(s) = 3.773 GeV : D0 -> K- pi+ pi0 pi0
# Double-tag analysis, tag side = anti-D0 -> K+ pi-
# =====================================================================

### Dataset preparation ###
d3770_data  = DatasetManager.real_data.find("712_3773")      # psi(3770) real data (2.93 fb^-1)
d3770_incMC = DatasetManager.inclusive_mc.find("712_3773")    # matching inclusive MC sample

# Signal decay card: full psi(3770) -> D0 anti-D0 chain (EvtGen syntax)
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0                   PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ pi0 pi0               PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi-                       PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                  PHSP;
    Enddecay

    End
DECAYCARD

# 1M exclusive MC events for the full signal chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKmPipPi0Pi0_Dbar0toKpPim"
  config.related_dataset = d3770_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: "temp_for_test")

### Event selection (BOSS) — tag-based (double tag) ###
alg_name = "D0toKmPipPi0Pi0DTAg"
my_alg = TagAnalysis.new(alg_name)
my_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
      .set_constant({"ECMS" => [:double, 3.773]})

# --- Tag side: anti-D0 -> K+ pi- (charm = -1 pins the anti-D0) ---
# Explicit opt-in tag-side windows: Delta E in [-0.03, 0.02] GeV,
# m_BC in [1.8575, 1.8775] GeV/c^2.  All other tag observables stay stored-not-cut.
my_alg.tag_side(:D0) do |t|
  t.modes  :D0toKPi
  t.charm  -1
  t.window :deltaE, min: -0.03,  max: 0.02
  t.window :mBC,    min: 1.8575, max: 1.8775
end

# --- Signal side: D0 -> K- pi+ pi0 pi0  (two pi0 -> 4 photons) ---
# Charged multiset fixed to exactly one K- and one pi+.
my_alg.signal_side do |s|
  s.photons 4                      # 4 signal-side photons (2 x pi0 -> gamma gamma)
  s.charged(km: 1, pip: 1)         # exactly one K- and one pi+
  s.min_photon_angle  10.0         # photon angle to nearest charged track > 10 deg
  s.min_photon_energy 0.025        # photon energy > 25 MeV
end

# --- Combined tag + signal kinematic fit: 4C + 2 x m(gamma gamma)=m(pi0) = 6C ---
my_alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # pi0 (1)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # pi0 (2)
  f.chi2_cut 80
end

# --- Paper-level signal-side items: recorded here, applied in the ROOT analysis ---
my_alg
  .note(:signal_vertex_fit,
        "Paper applies a signal-side K- pi+ secondary vertex fit with chi2 < 100; the " \
        "implemented BOSS step combines tag and signal in a single 6C kinematic fit " \
        "(4C four-momentum + two pi0 mass constraints) and does not contain this " \
        "signal-side vertex fit.")
  .note(:mbc_deltae_window,
        "Paper applies signal-side m_BC in [1.8600, 1.8730] GeV/c^2 and Delta E in " \
        "[-0.04, 0.02] GeV; these are stored by the tag/signal reconstruction and " \
        "windowed in the ROOT analysis, not cut at BOSS level.")
  .note(:pi0_mass_window,
        "Paper applies a diphoton invariant-mass window [0.115, 0.150] GeV/c^2 with at " \
        "least one photon in the EMC barrel; applied as a ROOT-level cut.")
  .note(:background_veto,
        "K_S0 veto: events with M(pi+ pi-) in [0.458, 0.520] GeV/c^2 are rejected to " \
        "suppress K_S0 -> pi+ pi- contamination.")
  .note(:photon_timing,
        "Signal-side photon EMC timing required below 700 ns; not tunable through the " \
        "signal-side shower declarations (DTagTool isGoodShower defaults used).")
  .note(:efficiency_curve,
        "Signal-side photon energy thresholds: 25 MeV (barrel) / 50 MeV (endcap); " \
        "minimum photon angle to charged tracks 10 degrees.")
  .note(:fit_hypothesis,
        "Implemented fit is the combined tag + signal 6C fit (4C four-momentum plus two " \
        "pi0 mass constraints) rather than the paper's signal-side 3C fit.")
  .note(:background_veto_charge,
        "Description requests a total signal-side charge of -1, which is inconsistent " \
        "with the declared K- pi+ multiset (net charge 0); the exact charged multiset " \
        "{K-: 1, pi+: 1} is used and no explicit charge pin is emitted.")

# Render the tag specification (apply takes no Selection) and run on the datasets
my_alg.with_decay_card(decay_card_signal).apply
root_files = my_alg.execute_on([d3770_data, d3770_incMC, exMC_signal])