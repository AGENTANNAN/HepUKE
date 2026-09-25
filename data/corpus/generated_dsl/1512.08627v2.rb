### Dataset description ###
# psi(3770) real data (2.93 fb^-1) and the matching inclusive MC sample.
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the exclusive signal MC: psi(3770) -> D+ D-, with the signal
# D+ -> K- pi+ e+ nu_e generated with PHOTOS and phase space.
# The opposite D- (tag side) is generated with a generic hadronic decay; the
# DTagTool reconstructs whichever of the declared tag modes applies.
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 K- pi+ e+ nu_e PHOTOS PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive signal MC events for the efficiency.
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToKPiENu"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) — tagged-D analysis ###
alg_name = "DpToKPiENuTag"
tag_alg  = TagAnalysis.new(alg_name)
tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.773]})       # sqrt(s) = 3.773 GeV
       .set_alias({"std::vector<double>" => "Vdouble"})

# --- Tag side: single tag of the D- in six hadronic modes ---
tag_alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,        # D- -> K+ pi- pi-
          :DptoKPiPiPi0,     # D- -> K+ pi- pi- pi0
          :DptoKsPi,         # D- -> K_S0 pi-
          :DptoKsPiPi0,      # D- -> K_S0 pi- pi0
          :DptoKsPiPiPi,     # D- -> K_S0 pi- pi- pi+
          :DptoKKPi          # D- -> K+ K- pi-
  t.charm -1                 # pin the tagged side to the D-
  # Only the M_BC window is applied here; the DeltaE window (~2 sigma around
  # zero) is applied downstream in ROOT (store-not-cut).
  t.window :mBC, min: 1.863, max: 1.877
end

# --- Signal side: K- pi+ e+ together with a missing massless nu_e ---
tag_alg.signal_side do |s|
  s.charged(km: 1, pip: 1, ep: 1)  # exactly one K-, one pi+ and one e+
  s.require_charge 1               # net charge of the signal side is +1
  s.missing :nu_e                  # massless missing neutrino (semileptonic tag)
  s.min_photon_angle 10.0          # minimum photon angle to the closest track
end

# --- 4C kinematic fit ---
tag_alg.fit do |f|
  f.constrain_four_momentum        # tag + K- pi+ e+ + nu constrained to the measured CMS
  f.chi2_cut 200                   # loose chi2 here; tighter cut applied downstream in ROOT
end

# --- BOSS-side procedures with no formal DSL construct ---
tag_alg
  .note(:tag_track_photon_quality,
        "Tag-side good tracks (|Vz| < 10 cm, |Vr| < 1 cm, |cos(theta)| < 0.93) and good photons "
        "(E > 25 MeV in the barrel |cos(theta)| < 0.80, E > 50 MeV in the endcap "
        "0.86 < |cos(theta)| < 0.92, shower time within 700 ns) are the criteria applied inside "
        "DTagAlg's tag reconstruction (findSTag/findDTag); the generated selection does not "
        "re-apply them.")
  .note(:tag_pid,
        "Tag-side pi/K separation from combined dE/dx and TOF likelihoods (P(K) > P(pi) for "
        "kaons, P(pi) > P(K) for pions) is performed inside DTagAlg's tag reconstruction and is "
        "not re-applied by the generated selection.")
  .note(:tag_pi0_ks_reconstruction,
        "Tag-side pi0 candidates (0.115 < M(gamma gamma) < 0.150 GeV/c^2 with a 1C "
        "mass-constrained fit chi2 < 200, events with both photons in the endcap vetoed) and "
        "K_S0 candidates (0.487 < M(pi+ pi-) < 0.511 GeV/c^2 with the secondary vertex within "
        "20 cm of the IP along the beam and no good-track or PID requirement on the daughters) "
        "are built inside the DTagAlg tag reconstruction.")
  .note(:tag_candidate_selection,
        "At most one tag candidate per tag mode and charge is kept, namely the one with the "
        "smallest |DeltaE|; this per-mode/per-charge best-candidate choice has no DSL primitive "
        "for a single-tag pattern. DTagAlg's own mBC/DeltaE selector cuts stay disabled in the "
        "generated jobOptions: the M_BC window 1.863-1.877 GeV/c^2 declared above is the only "
        "tag-side window applied, and |DeltaE| is windowed downstream in ROOT.")
  .note(:background_veto,
        "D0 cross-feed veto for the K_S0 pi- pi- pi+, K_S0 pi- pi0 and K+ pi- pi- pi0 tag modes: "
        "events with 1.860 < M_BC(D0) < 1.875 GeV/c^2 and |DeltaE(D0)| < 0.01 GeV are rejected; "
        "no DSL primitive expresses a tag-side cross-feed veto.")
  .note(:signal_electron_pid,
        "The signal-side positron is identified by the charged(ep: 1) key using the fixed "
        "SimplePIDSvc lepton recipe. The analysis-specific electron requirement - combined "
        "dE/dx + TOF + EMC likelihoods with P^2(e)/(P^2(K)+P^2(pi)+P^2(e)) > 0.8, P^2(e) > 0.001, "
        "and EMC energy above 80% of the MDC momentum - is not DSL-tunable and has to be imposed "
        "by retuning the generated tag code.")
  .note(:extra_shower_veto,
        "Signal-side veto on unassociated EMC showers above 0.25 GeV, where unassociated clusters "
        "are defined as more than 15 degrees from the closest charged track; this veto has no DSL "
        "primitive. The declared min_photon_angle 10.0 is the good-shower isolation angle passed "
        "to DTagTool, while the 15-degree cluster definition is applied downstream.")
  .note(:missing_mass_selection,
        "The signal-side selections |U_miss| = |E_miss - |p_miss|| < 0.04 GeV and "
        "E_miss > 0.04 GeV use the auto-stored semileptonic diagnostics (m_Umiss, m_Umiss2, "
        "m_P4_miss_fit) and are applied downstream in ROOT rather than at BOSS level.")

# Render the tag analysis (takes NO Selection argument) and run it.
tag_alg.with_decay_card(decay_card_signal).apply
root_files = tag_alg.execute_on([data_3773, incMC_3773, exMC_signal])