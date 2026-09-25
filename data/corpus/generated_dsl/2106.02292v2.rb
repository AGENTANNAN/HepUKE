# frozen_string_literal: true

### Dataset description ###
# ψ(3770) real data and the corresponding inclusive MC sample
# (sample-name convention: [BOSS version]_[CMS energy in MeV])
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the signal process ψ(3770) → D0 D0bar,
# with D0 → ρ− μ+ ν_μ, ρ− → π− π0 → γγ; the anti-D0 (tag side) decays into the
# three hadronic tag modes of the analysis.
decay_card_signal = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0000 rho- mu+ nu_mu PHSP;
  Enddecay

  Decay rho-
  1.0000 pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay anti-D0
  0.3333 K+ pi- PHSP;
  0.3333 K+ pi- pi0 PHSP;
  0.3334 K+ pi- pi- pi+ PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive signal MC: 1,000,000 ψ(3770) → D0 D0bar events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0_rhomunu"
  config.related_dataset = data_3773            # generated for the 712_3773 data taking
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
alg_name = "D0TagRhoMuNu"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})   # √s = 3.773 GeV (ψ(3770))

# BOSS-side procedures that have no formal DSL construct are kept as notes.
alg
  .note(:deltaE_mode_dependence,
        "Tag-side ΔE windows are mode dependent: (−0.055, 0.040) GeV for " \
        "D0bar → K+π−π0 and (−0.025, 0.025) GeV for D0bar → K+π− and " \
        "D0bar → K+π−π−π+. The DSL emits the lossless union window " \
        "(−0.055, 0.040) GeV; the per-mode tightening is applied on the stored " \
        "tag candidates in the ROOT analysis.")
  .note(:tag_candidate_selection,
        "For each tag mode and each tag charge the candidate with the smallest " \
        "|ΔE| is retained (single tag per mode per charge); the ranking is done " \
        "on the stored tag candidates in the ROOT analysis.")
  .note(:pid_correction_method,
        "Signal μ+ identified from combined dE/dx + TOF + EMC probabilities with " \
        "the MUC not used: CL_μ > 0.001, CL_μ > CL_e, CL_μ > CL_K and no extra " \
        "μ/π separation requirement; the generated tag code applies its fixed " \
        "lepton-PID recipe, retuning implies editing the generated code.")
  .note(:background_veto,
        "K_S0 veto: events with M(π−π+) ∈ 0.458–0.538 GeV/c² or recoil " \
        "M(π+π−) ∈ 0.468–0.528 GeV/c² are rejected; K→π misidentification veto: " \
        "|M²_miss(π→K)| > 0.05 GeV²/c⁴; no extra charged track and no extra π0 " \
        "are allowed.")
  .note(:signal_side_windows,
        "Signal-side windows applied on the measured (pre-fit) quantities: " \
        "E_μ,EMC ∈ (0.125, 0.275) GeV, M(ρ−μ+) < 1.5 GeV/c² and extra-photon " \
        "energy < 0.25 GeV; the final signal is extracted from the " \
        "missing-mass-squared distribution M²_miss (stored by the fit).")

# One tag side → single tag; the missing ν_μ on the signal side makes it ST + missing.
# The tagged particle is the hadronic anti-D0 (charm = −1), reconstructed in three
# DTagAlg modes. Tag mBC / ΔE are stored unconditionally; only the explicitly
# requested windows are emitted.
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi   # K+π−, K+π−π0, K+π−π−π+ (charm −1)
  t.charm -1
  t.window :mBC,    min: 1.859, max: 1.873     # M_BC ∈ (1.859, 1.873) GeV/c²
  t.window :deltaE, min: -0.055, max: 0.040    # lossless union of the per-mode ΔE windows
end

# Signal side: two photons (π0 → γγ), one π−, one μ+, net charge zero,
# and one massless missing ν_μ.
alg.signal_side do |s|
  s.photons 2                 # two showers forming the π0
  s.charged(pim: 1, mup: 1)   # exactly one π− and one μ+ (μ+ via the lepton PID recipe)
  s.require_charge 0          # net charge of the signal side
  s.missing :nu_mu            # massless missing muon neutrino
end

# Kinematic fit: tag ⊗ signal(π−, μ+, γγ) ⊗ ν_μ constrained to the ψ(3770) four-momentum,
# with M(γγ) constrained to the nominal π0 mass.
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # π0 → γγ
  f.invariant_mass_of(:pim, :gamma, :gamma).between(0.625, 0.925)        # M(π−π0) window
  f.chi2_cut 200                                                          # loose BOSS-level χ²
end

# Render the tag analysis (no Selection argument) and run on data, inclusive MC and signal MC
alg.apply
root_files = alg.execute_on([data_3773, incMC_3773, exMC_signal])