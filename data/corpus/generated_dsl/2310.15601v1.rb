# =============================================================================
# BESIII  e+e- -> Ds*+ Ds-  (eight points: 4.128, 4.157, 4.178, 4.189, 4.199,
#                                  4.209, 4.219, 4.226 GeV ; BOSS 705 / 703)
# Doubly Cabibbo-suppressed decays
#   Ds+ -> K+ K+ pi-        and      Ds+ -> K+ K+ pi- pi0
# Hadronic tag on the Ds-, with the Ds*+ -> gamma Ds+ / pi0 Ds+ transition.
# Tag-based analysis -> TagAnalysis surface (no Selection / track / PID blocks).
# =============================================================================

### Datasets — eight real-data points and the matching inclusive MC ###
data_points = [
  DatasetManager.real_data.find("705_4130"),   # ~4.128 GeV
  DatasetManager.real_data.find("705_4160"),   # ~4.157 GeV
  DatasetManager.real_data.find("703_4180"),   # ~4.178 GeV
  DatasetManager.real_data.find("703_4190"),   # ~4.189 GeV
  DatasetManager.real_data.find("703_4200"),   # ~4.199 GeV
  DatasetManager.real_data.find("703_4210"),   # ~4.209 GeV
  DatasetManager.real_data.find("703_4220"),   # ~4.219 GeV
  DatasetManager.real_data.find("703_4230")    # ~4.226 GeV
]

incMC_points = [
  DatasetManager.inclusive_mc.find("705_4130"),
  DatasetManager.inclusive_mc.find("705_4160"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230")
]

### Decay cards (EvtGen) — one per Ds+ signal mode, generated in phase space ###
# Mode I : Ds+ -> K+ K+ pi-
decay_card_kkpi = <<~DECAYCARD
    Decay psi(4260)
    1.0 Ds*+ Ds- PHSP;
    Enddecay

    Decay Ds*+
    1.0 gamma Ds+ PHSP;
    Enddecay

    Decay Ds+
    1.0 K+ K+ pi- PHSP;
    Enddecay

    Decay Ds-
    1.0 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode II : Ds+ -> K+ K+ pi- pi0, pi0 -> gamma gamma
decay_card_kkpipi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0 Ds*+ Ds- PHSP;
    Enddecay

    Decay Ds*+
    1.0 gamma Ds+ PHSP;
    Enddecay

    Decay Ds+
    1.0 K+ K+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    Decay Ds-
    1.0 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive signal MC — same signal generated at every one of the eight energies ###
# (no explicit event count was specified; the generator default is used)
exMCs_kkpi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_dcs_ds_kkpi"
  config.decay_card    = decay_card_kkpi
  config.cross_section = :default
end

exMCs_kkpipi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_dcs_ds_kkpipi0"
  config.decay_card    = decay_card_kkpipi0
  config.cross_section = :default
end

# =============================================================================
# Mode I : Ds+ -> K+ K+ pi-
# =============================================================================
alg_name_kkpi = "DsDcsKKPi"
alg_kkpi = TagAnalysis.new(alg_name_kkpi)
alg_kkpi.set_header(["#{alg_name_kkpi}Alg/#{alg_name_kkpi}.h"])
        # ECMS seed value only — the 4-momentum constraint builds the CMS
        # four-vector per run from the measured beam energy (beam_energy :db).
        .set_constant({"ECMS" => [:double, 4.178]})

# Tag Ds-: hadronic Ds -> K+K-pi mode, covering the K+K-pi- and K+K-pi+
# charge-conjugate final states.
alg_kkpi.tag_side(:Ds) do |t|
  t.modes :DstoKKPi                      # K+ K- pi- and K+ K- pi+ final states
  t.window :mBC, min: 2.010, max: 2.088  # union of the eight per-energy M_BC windows
end

# Signal side: the tracks not used by the tag -> DCS Ds+ candidate.
alg_kkpi.signal_side do |s|
  s.charged(kp: 2, pim: 1)               # exactly two K+ and one pi-
  s.photons 0..48                        # 0 photons required for the K+K+pi- mode
end

# 6C kinematic fit: 4-momentum conservation + both Ds masses to the nominal Ds mass.
alg_kkpi.fit do |f|
  f.constrain_four_momentum                                            # 4C
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)         # tag Ds- mass (1C)
  f.invariant_mass_of(:kp, :kp, :pim).constrain_to_nominal_mass_of(:Ds) # signal Ds+ mass (1C)
  f.chi2_cut 200
end

# Tagging defaults: local re-reconstruction on, beam energy from the conditions DB,
# DTagAlg mode lists trimmed to the declared modes.
alg_kkpi.dtag_reconstruction do |r|
  r.local true
  r.beam_energy :db
  r.trim_mode_lists true
end

alg_kkpi
  .note(:tag_mode_availability,
        "Only the Ds -> K K pi hadronic tag mode (:DstoKKPi) is available in the frozen "
        "BOSS releases; the other paper tag modes (K+K-pi-pi0, pi-pi+pi-, K_S0 K-, "
        "eta / eta' / rho modes) are unavailable and are not reconstructed.")
  .note(:dsstar_transition,
        "Ds*+ -> gamma Ds+ and Ds*+ -> pi0 Ds+ hypotheses are resolved by looping over the "
        "unused gamma and pi0 candidates and choosing the combination with minimum |deltaE|.")
  .note(:tag_candidate_ranking,
        "Tag / double-tag candidates are ranked by invariant-mass difference (DTagTool :inv).")
  .note(:mbc_windows,
        "Per-energy M_BC windows applied in ROOT: 2.010-2.061 (4.128 GeV), 2.010-2.070 (4.157), "
        "2.010-2.073 (4.178), 2.010-2.076 (4.189), 2.010-2.079 (4.199), 2.010-2.082 (4.209), "
        "2.010-2.085 (4.219), 2.010-2.088 GeV (4.226); the BOSS-level window is their union.")

alg_kkpi.with_decay_card(decay_card_kkpi).apply
alg_kkpi.execute_on(data_points + incMC_points + exMCs_kkpi)

# =============================================================================
# Mode II : Ds+ -> K+ K+ pi- pi0  (pi0 -> gamma gamma)
# =============================================================================
alg_name_kkpipi0 = "DsDcsKKPiPi0"
alg_kkpipi0 = TagAnalysis.new(alg_name_kkpipi0)
alg_kkpipi0.set_header(["#{alg_name_kkpipi0}Alg/#{alg_name_kkpipi0}.h"])
           .set_constant({"ECMS" => [:double, 4.178]})   # per-run beam energy read from DB

alg_kkpipi0.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.window :mBC, min: 2.010, max: 2.088
end

alg_kkpipi0.signal_side do |s|
  s.charged(kp: 2, pim: 1)               # exactly two K+ and one pi-
  s.photons 2..48                        # at least two photons (pi0 -> gamma gamma)
end

alg_kkpipi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)                        # tag Ds- (1C)
  f.invariant_mass_of(:kp, :kp, :pim, :gamma, :gamma).constrain_to_nominal_mass_of(:Ds) # signal Ds+ (1C)
  f.chi2_cut 200
end

alg_kkpipi0.dtag_reconstruction do |r|
  r.local true
  r.beam_energy :db
  r.trim_mode_lists true
end

alg_kkpipi0
  .note(:tag_mode_availability,
        "Only the Ds -> K K pi hadronic tag mode (:DstoKKPi) is available in the frozen BOSS "
        "releases; the remaining paper tag modes are not reconstructed.")
  .note(:dsstar_transition,
        "Ds*+ -> gamma Ds+ and Ds*+ -> pi0 Ds+ hypotheses resolved by looping over the unused "
        "gamma / pi0 candidates and choosing minimum |deltaE|.")
  .note(:mbc_windows,
        "Per-energy M_BC windows (2.010-2.061 ... 2.010-2.088 GeV) applied in ROOT.")

alg_kkpipi0.with_decay_card(decay_card_kkpipi0).apply
alg_kkpipi0.execute_on(data_points + incMC_points + exMCs_kkpipi0)