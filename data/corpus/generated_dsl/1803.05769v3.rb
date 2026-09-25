# =====================================================================
# ψ(3770) → D0 D̄0 double-tag analysis
#   tag side   : anti-D0 → K+π−, K+π−π0, K+π−π−π+  (DTagAlg candidates)
#   signal side: all-neutral D0 → π0π0π0, π0π0η, π0ηη, ηηη  (π0, η → γγ)
# =====================================================================

### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")      # 3.773 GeV real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # matching inclusive MC

# --- decay cards: four all-neutral D0 signal modes in phase space ---
# a representative charged tag (D̄0 → K+π−) is included; the tag modes
# K+π−, K+π−π0, K+π−π−π+ are all taken from DTagAlg tag candidates.
dcard_3pi0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi0 pi0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

dcard_2pi0eta = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi0 pi0 eta PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

dcard_pi0etaeta = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi0 eta eta PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

dcard_3eta = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 eta eta eta PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- 100k-event exclusive MC for each signal mode ---
exMC_3pi0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3773_D0to3pi0"
    config.related_dataset = data_3773
    config.events          = 100_000
    config.decay_card      = dcard_3pi0
    config.cross_section   = :default
end

exMC_2pi0eta = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3773_D0to2pi0eta"
    config.related_dataset = data_3773
    config.events          = 100_000
    config.decay_card      = dcard_2pi0eta
    config.cross_section   = :default
end

exMC_pi0etaeta = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3773_D0topi0etaeta"
    config.related_dataset = data_3773
    config.events          = 100_000
    config.decay_card      = dcard_pi0etaeta
    config.cross_section   = :default
end

exMC_3eta = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3773_D0to3eta"
    config.related_dataset = data_3773
    config.events          = 100_000
    config.decay_card      = dcard_3eta
    config.cross_section   = :default
end

### Event selection (tag-based, BOSS) ###
alg_name = "D0AllNeutralDT"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(dcard_3pi0)

# --- tag side: charged anti-D0 (K+π−, K+π−π0, K+π−π−π+) ---
alg.tag_side(:D0) do |t|
    t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
    t.charm -1                    # tag the anti-D0, opposite the D0 signal
    t.window :mBC, min: 1.83      # tag mBC > 1.83 GeV
end

# --- signal side: all-neutral D0 (six photons, no charged track) ---
alg.signal_side do |s|
    s.photons 6
    s.min_photon_angle  10.0
    s.min_photon_energy 0.025     # barrel E > 25 MeV (endcap 50 MeV, see note)
end

# --- kinematic fit: 4C plus one γγ mass constraint (π0, or η for ηηη) ---
alg.fit do |f|
    f.constrain_four_momentum
    f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    f.chi2_cut 200
end

# --- inexpressible BOSS-side procedures ---
alg.note(:photon_multiplicity,
         "The all-neutral signal D0 → π0π0π0 / π0π0η / π0ηη / ηηη requires six " \
         "photons (π0, η → γγ). The generated code consumes only up to two signal " \
         "photons in the kinematic fit; the full six-photon requirement is applied " \
         "at the ROOT level.")
   .note(:mass_constraints,
         "Three 1C π0/η mass constraints (one per γγ pair) are needed to fully " \
         "reconstruct the all-neutral D0; only one γγ mass constraint is " \
         "expressible in the kinematic-fit block.")
   .note(:background_veto,
         "Vetoes: D0 → 4π0 (χ²_4π > 20 when at least four independent π0 " \
         "candidates); K_S0 → π0π0 (M(π0π0) in [445, 535] MeV); cross-feed " \
         "χ²_3π0 > 20 and χ²_π0π0η > 20 with no other combination having χ² < 20.")
   .note(:selection_windows,
         "Tag ΔE windows: Kπ (−0.027, 0.025), Kππ0 (−0.071, 0.041), " \
         "Kπππ (−0.025, 0.022) GeV. Signal ΔE windows: 3π0 (−0.115, 0.059), " \
         "2π0η (−0.088, 0.053), π0ηη (−0.061, 0.045), 3η (−0.030, 0.028) GeV; " \
         "best candidate chosen by minimum |ΔE|. The per-mode windows are stored " \
         "unconditionally and applied at the ROOT level.")
   .note(:photon_selection,
         "Signal photons: E > 25 MeV in the barrel (|cosθ| < 0.8) and E > 50 MeV in " \
         "the endcap (0.86 < |cosθ| < 0.92); EMC TDC in [0, 700] ns and a minimum " \
         "photon angle of 10°. Only the single energy floor and the 10° angle are " \
         "expressible in the signal-side block.")
   .note(:tag_track_selection,
         "Tag-side charged tracks require |cosθ| < 0.93, Vr < 1 cm and |Vz| < 10 cm, " \
         "with PID from TOF and dE/dx selecting the highest confidence level; these " \
         "are applied internally by DTagTool and are not re-expressed in the DSL.")
   .note(:decay_model,
         "The π0ηη signal mode is dominated by the a0(980)0η intermediate state; a " \
         "phase-space card is used for this exclusive MC sample.")

alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_3pi0, exMC_2pi0eta, exMC_pi0etaeta, exMC_3eta])