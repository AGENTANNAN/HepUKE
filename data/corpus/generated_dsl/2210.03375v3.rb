# =============================================================================
# Cabibbo-suppressed Λc+ → nπ+π0, nπ+π-π+, nK-π+π+
# Double-tag analysis: the event is tagged by a single-tag (ST) Λc- built from
# the pre-stored DTag candidates (DTagTool); the Λc+ is reconstructed in the
# system recoiling against the tag and the signal neutron is identified through
# the missing mass.  Tag-based spec → TagAnalysis (no Selection object).
# =============================================================================

### ------------------------- Dataset preparation ------------------------- ###
# 7 c.m. energy points, 4.5995–4.6988 GeV (BOSS 703 for the 4.600 GeV point,
# BOSS 706 for the six 4.610–4.700 GeV points).
data_points = [
  DatasetManager.real_data.find("703_4600"),   # 4599.53 MeV
  DatasetManager.real_data.find("706_4610"),   # 4611.86 MeV
  DatasetManager.real_data.find("706_4620"),   # 4628.00 MeV
  DatasetManager.real_data.find("706_4640"),   # 4640.91 MeV
  DatasetManager.real_data.find("706_4660"),   # 4661.24 MeV
  DatasetManager.real_data.find("706_4680"),   # 4681.92 MeV
  DatasetManager.real_data.find("706_4700")    # 4698.82 MeV
]

# Matching inclusive MC samples (background estimation / efficiency study)
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700")
]

datasets = data_points + incMC_points

### ----------------------------- Decay cards ----------------------------- ###
# Mode 1: Λc+ → n π+ π0
decay_card_mode1 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ Lambda_c-        PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 n0 pi+ pi0                PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma               PHSP;
    Enddecay

    End
DECAYCARD

# Mode 2: Λc+ → n π+ π- π+
decay_card_mode2 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ Lambda_c-        PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 n0 pi+ pi- pi+            PHSP;
    Enddecay

    End
DECAYCARD

# Mode 3: Λc+ → n K- π+ π+
decay_card_mode3 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ Lambda_c-        PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 n0 K- pi+ pi+             PHSP;
    Enddecay

    End
DECAYCARD

### ------------------------ Λc- ST tag mode lists ------------------------ ###
# Full 12-mode Λc- tag list declared by the analysis:
#   pK-π0, pπ+π-, pK-π0π0, pK-π0π+π-, pK-π0π+π-π0, Λπ+, Λπ+π0,
#   Λπ+π-π+, Σ+π0, Σ+π0π0, Σ+π+π-, pπ-π+
lambdac_tag_modes = [
  :LambdacPtoPKPi0,            # p K- π0
  :LambdacPtoPPiPi,            # p π+ π-
  :LambdacPtoPKPi0Pi0,         # p K- π0 π0
  :LambdacPtoPKPi0PiPi,        # p K- π0 π+ π-
  :LambdacPtoPKPi0PiPiPi0,     # p K- π0 π+ π- π0
  :LambdacPtoLambdaPi,         # Λ π+
  :LambdacPtoLambdaPiPi0,      # Λ π+ π0
  :LambdacPtoLambdaPiPiPi,     # Λ π+ π- π+
  :LambdacPtoSigmaPPi0,        # Σ+ π0
  :LambdacPtoSigmaPPi0Pi0,     # Σ+ π0 π0
  :LambdacPtoSigmaPPiPi,       # Σ+ π+ π-
  :LambdacPtoPPiPiCC           # p π- π+
]

# Mode 1 (nπ+π0): 10 tags — pK+π-π0 and pπ-π+ excluded
modes_mode1 = lambdac_tag_modes - [:LambdacPtoPPiPi, :LambdacPtoPPiPiCC]
# Mode 2 (nπ+π-π+): 11 tags — pπ-π+ excluded
modes_mode2 = lambdac_tag_modes - [:LambdacPtoPPiPiCC]
# Mode 3 (nK-π+π+): all 12 tags
modes_mode3 = lambdac_tag_modes

# =============================================================================
# Mode 1 : Λc+ → n π+ π0
# =============================================================================
alg_mode1 = TagAnalysis.new("LambdacTagNpiPi0")
alg_mode1.set_header(["LambdacTagNpiPi0Alg/LambdacTagNpiPi0.h"])
         .set_constant({ "ECMS" => [:double, 4.600] })   # reference point; per-run measured beam energy is used by the fit
         .with_decay_card(decay_card_mode1)

# One tag side → single tag (ST) of Λc-, charm = -1
alg_mode1.tag_side(:Lambdac) do |t|
  t.modes(*modes_mode1)   # 10 Cabibbo-allowed / -suppressed Λc- tag modes
  t.charm(-1)             # pin the tagged side to Λc-
end

# Signal side = everything the tag did not use
alg_mode1.signal_side do |s|
  s.photons 2                      # two photons to form the π0
  s.min_photon_angle 10.0          # minimum photon–track separation angle (deg)
  s.charged(pip: 1)                # π+ (total signal-side charge +1)
  s.require_charge 1
  s.missing :n                     # neutron identified via missing mass
end

# Kinematic fit: 4-momentum constraint, χ² < 200; the π0 is required via its mass
alg_mode1.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_mode1
  .note(:background_veto,
        "Peaking background from Λc+ → Λπ+ / Σ+π0 / Σ0π+ feed-down into " \
        "Λc+ → nπ+π0 removed by missing-mass cuts Mmiss(π+) > 1300 MeV/c² and " \
        "Mmiss(π0) > 1370 MeV/c²; the missing mass is stored per event by the " \
        "tag fit (m_P4_miss_fit) and the veto is applied in the ROOT analysis.")
  .note(:signal_window,
        "Signal extraction windows are applied in ROOT (store-not-cut): " \
        "MBC ∈ [2.275, 2.300] GeV/c² at 4.600 GeV, [2.275, 2.306] GeV/c² at " \
        "4.612–4.641 GeV and [2.275, 2.310] GeV/c² at 4.661–4.699 GeV, " \
        "together with the mode-dependent (~3σ) ΔE tag windows.")

alg_mode1.apply

# =============================================================================
# Mode 2 : Λc+ → n π+ π- π+
# =============================================================================
alg_mode2 = TagAnalysis.new("LambdacTagNpiPiPi")
alg_mode2.set_header(["LambdacTagNpiPiPiAlg/LambdacTagNpiPiPi.h"])
         .set_constant({ "ECMS" => [:double, 4.600] })
         .with_decay_card(decay_card_mode2)

alg_mode2.tag_side(:Lambdac) do |t|
  t.modes(*modes_mode2)   # 11 Λc- tag modes
  t.charm(-1)
end

alg_mode2.signal_side do |s|
  s.charged(pip: 2, pim: 1)   # π+ π- π+  (total charge +1)
  s.require_charge 1
  s.missing :n                # neutron identified via missing mass
end

alg_mode2.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mode2
  .note(:background_veto,
        "Peaking backgrounds in Λc+ → nπ+π-π+ removed by vetos applied on the " \
        "stored observables in the ROOT analysis: K_S0 → π+π- (M(π+π-) outside " \
        "[487, 511] MeV/c²), Σ+ → nπ+ (Mmiss(π+π-) outside [1150, 1250] MeV/c²) " \
        "and Σ- → nπ- (Mmiss(π+π+) outside [1150, 1250] MeV/c²).")
  .note(:signal_window,
        "Signal extraction windows applied in ROOT: MBC ∈ [2.275, 2.300] GeV/c² " \
        "at 4.600 GeV, [2.275, 2.306] GeV/c² at 4.612–4.641 GeV and " \
        "[2.275, 2.310] GeV/c² at 4.661–4.699 GeV, with mode-dependent (~3σ) " \
        "ΔE tag windows.")

alg_mode2.apply

# =============================================================================
# Mode 3 : Λc+ → n K- π+ π+
# =============================================================================
alg_mode3 = TagAnalysis.new("LambdacTagNKPiPi")
alg_mode3.set_header(["LambdacTagNKPiPiAlg/LambdacTagNKPiPi.h"])
         .set_constant({ "ECMS" => [:double, 4.600] })
         .with_decay_card(decay_card_mode3)

alg_mode3.tag_side(:Lambdac) do |t|
  t.modes(*modes_mode3)   # the full 12-mode Λc- tag list
  t.charm(-1)
end

alg_mode3.signal_side do |s|
  s.charged(km: 1, pip: 2)    # K- π+ π+  (total charge +1)
  s.require_charge 1
  s.missing :n                # neutron identified via missing mass
end

alg_mode3.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mode3
  .note(:signal_window,
        "Signal extraction windows applied in ROOT: MBC ∈ [2.275, 2.300] GeV/c² " \
        "at 4.600 GeV, [2.275, 2.306] GeV/c² at 4.612–4.641 GeV and " \
        "[2.275, 2.310] GeV/c² at 4.661–4.699 GeV, with the mode-dependent " \
        "(~3σ) ΔE tag windows.")

alg_mode3.apply

# =============================================================================
# Execution — 7 real-data points + the matching inclusive MC samples
# =============================================================================
root_files_mode1 = alg_mode1.execute_on(datasets)
root_files_mode2 = alg_mode2.execute_on(datasets)
root_files_mode3 = alg_mode3.execute_on(datasets)

# NOTE on ECMS: a single ECMS constant is declared per algorithm (4.600 GeV),
# while the spec is executed over the 7 c.m. energy points 4.5995-4.6988 GeV.
# The 4-momentum constraint of the tag fit uses the per-run measured CMS
# four-vector (MeasuredEcmsSvc) on data, so the running energy is handled
# correctly; ECMS only supplies the default beam-energy spread.