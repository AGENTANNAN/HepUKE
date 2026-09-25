# frozen_string_literal: true

### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")     # psi(3770) real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC

# Decay card: psi(3770) -> D0 anti-D0; signal D0 -> K+ K- pi0 pi0; generator card
# decays the anti-D0 to the tag mode K+ pi-.
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K+ K- pi0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive MC events for the signal decay D0 -> K+ K- pi0 pi0
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKKPi0Pi0"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
alg = TagAnalysis.new("D0KKPi0Pi0Tag")
alg.set_header(["D0KKPi0Pi0TagAlg/D0KKPi0Pi0Tag.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(decay_card_signal)

# Tag side: anti-D0 reconstructed from the pre-stored DTag candidates in three modes
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi  # K+pi-, K+pi-pi0, K+pi-pi-pi+ (charge conj.)
  t.charm -1                                   # pin the tagged side to the anti-D0
  # DeltaE is stored unconditionally; the loose BOSS-level envelope below covers all
  # tags, the exact (mode-dependent) windows are applied in ROOT (see note).
  t.window :deltaE, min: -0.055, max: 0.040
end

# Signal side: D0 -> K+ K- pi0 pi0, built from the showers/tracks not used by the tag
alg.signal_side do |s|
  s.charged(kp: 1, km: 1)   # exactly one K+ and one K-
  s.require_charge(0)       # net charge zero
  s.photons 4               # four photons (two pi0 candidates)
  s.min_photon_energy 0.025 # photon energy above 25 MeV
  s.min_photon_angle 10.0   # photon opening angle above 10 degrees
end

# Kinematic fit: 4-momentum conservation + two pi0 mass constraints to the nominal pi0 mass
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # first pi0 (1C)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # second pi0 (1C)
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)               # pi0 candidate mass window
  f.chi2_cut 200
end

# BOSS-side procedures that cannot be expressed in the tag DSL (applied in ROOT)
alg.note(:tag_deltaE_windows,
         "Tag-side DeltaE windows are mode dependent: K+pi- (-0.025, 0.025) GeV, " \
         "K+pi-pi0 (-0.055, 0.040) GeV, K+pi-pi-pi+ (-0.025, 0.025) GeV; a single loose " \
         "envelope is encoded in BOSS, the exact per-mode windows are applied in ROOT. " \
         "Signal-side DeltaE window (-0.020, 0.020) GeV is also applied in ROOT.")
   .note(:amplitude_analysis_region,
         "Amplitude-analysis region 1.861 < MBC < 1.870 GeV/c2 for both signal and tag; " \
         "MBC is stored unconditionally in BOSS and the region is applied in ROOT.")
   .note(:pi0_selection,
         "pi0 candidates must satisfy 0.115 < M(gamma gamma) < 0.150 GeV/c2, contain at " \
         "least one barrel photon, and pass a 1C Kalman pi0 mass-constrained fit with " \
         "chi2 < 50 (mass constraint handled by the Kalman fit, not the tag fit).")
   .note(:background_veto,
         "K_S0 veto: reject events with M(pi0 pi0) in [0.460, 0.525] GeV/c2; additionally " \
         "require the D0 - anti-D0 opening angle > 167 degrees.")
   .note(:best_candidate_selection,
         "When several candidates survive, the best candidate is chosen by minimum |DeltaE|.")

alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_signal])