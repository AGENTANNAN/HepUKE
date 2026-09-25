# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# ψ(3770) at √s = 3.773 GeV
psi3770_data  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

# Signal decay card: e+e- → D⁰D̄⁰; tag D̄⁰ → K⁺π⁻; signal D⁰ → K⁻π⁺π⁺π⁻ (phase space)
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0  D0  anti-D0            PHSP;
    Enddecay

    Decay anti-D0
    1.0  K+  pi-                PHSP;
    Enddecay

    Decay D0
    1.0  K-  pi+  pi+  pi-      PHSP;
    Enddecay

    End
DECAYCARD

# Peaking-background decay card: same tag D̄⁰ → K⁺π⁻; signal D⁰ → K_S⁰ K⁻π⁺
decay_card_peaking = <<~DECAYCARD
    Decay psi(3770)
    1.0  D0  anti-D0            PHSP;
    Enddecay

    Decay anti-D0
    1.0  K+  pi-                PHSP;
    Enddecay

    Decay D0
    1.0  K_S0  K-  pi+          PHSP;
    Enddecay

    Decay K_S0
    1.0  pi+  pi-               PHSP;
    Enddecay

    End
DECAYCARD

# 1,000,000 phase-space signal D⁰ → K⁻π⁺π⁺π⁻ exclusive MC events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_D0toKPiPiPi_signal"
  config.related_dataset = psi3770_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# 500,000 D⁰ → K_S⁰ K⁻π⁺ peaking-background exclusive MC events (same tag)
exMC_peaking = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_D0toKsKPi_peaking"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_peaking
  config.cross_section   = :default
end

### Event selection (tag-based D⁰D̄⁰ double tag) ###
alg = TagAnalysis.new("D0D0barDT")
alg.set_header(["D0D0barDTAlg/D0D0barDT.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })
   .with_decay_card(decay_card_signal)

# Tag side: D̄⁰ → K⁺π⁻, pinned to the D̄⁰ charm.
# Explicit opt-in tag-side windows (mBC, ΔE); all other tag observables are stored
# unconditionally and windowed in ROOT.
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi
  t.charm -1
  t.window :mBC,    min: 1.8575, max: 1.8775   # 1.8575 < mBC < 1.8775 GeV/c²
  t.window :deltaE, abs: 0.03                  # |ΔE| < 0.03 GeV
end

# Signal side: the D⁰ built from the tag's unused tracks — K⁻ π⁺ π⁺ π⁻, net charge zero.
# No photon cuts (fully charged final state); PID per species is handled by the signal side.
alg.signal_side do |s|
  s.charged(km: 1, pip: 2, pim: 1)
  s.require_charge 0
end

# 5C kinematic fit: four-momentum conservation to the initial e⁺e⁻ system +
# m(K⁻π⁺π⁺π⁻) constrained to the nominal D⁰ mass.
# Loose χ² < 200 in BOSS (the paper's χ² < 40 is imposed later in ROOT).
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:km, :pip, :pip, :pim).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# BOSS-side procedures that have no direct tag-DSL construct
alg.note(:d0d0bar_vertex_fit, "a common vertex fit on the D⁰ and D̄⁰ candidate tracks is required, with vertex χ² < 200")
alg.note(:background_veto, "signal-side π⁺π⁻ combinations with |m(π⁺π⁻) − m(K_S⁰)| < 0.03 GeV/c² and decay-length significance > 2σ are rejected to suppress the D⁰ → K_S⁰ K⁻π⁺ peaking background")

alg.apply
alg.execute_on([psi3770_data, psi3770_incMC, exMC_signal, exMC_peaking])