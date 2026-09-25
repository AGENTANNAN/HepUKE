# =============================================================================
# BOSS DSL — psi(3770) double-tag measurement of Cabibbo-suppressed hadronic
#             D0 and D+ decays into multi-pion final states
#
#   * data / incMC : dataset 712_3773 (sqrt(s) = 3.773 GeV)
#   * tag side     : hadronic D0bar / D- tags taken from DTagAlg pre-stored
#                    candidates
#   * signal side  : recoil D0 / D+ Cabibbo-suppressed multi-pion decays,
#                    reconstructed from the showers/tracks the tag did NOT use
#
#   Because a D0bar tag always recoils against a D0 and a D- tag always
#   recoils against a D+, the two event topologies (D0D0bar and D+D-) are
#   encoded as two independent TagAnalysis objects.
# =============================================================================

### ------------------------------ Datasets ------------------------------- ###
data_3773  = DatasetManager.real_data.find("712_3773")       # psi(3770) real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")    # matching inclusive MC

### --------------------- Exclusive signal MC decay cards ------------------ ###
# Three-body signal modes use the BODY3 generator; every multi-body mode with
# >= 4 bodies uses phase space (PHSP).  No per-mode event count was given, so a
# common default is used.  The signal D is generated inside
# psi(3770) -> D Dbar so that the correct beam-constrained kinematics is kept;
# the opposite D is left to a generic decay (it is only the flavour-tag side and
# is not used for the signal efficiency).
def build_cs_signal_card(d_mother, d_other, final_state, model)
  <<~DECAYCARD
    Decay psi(3770)
    1.0 #{d_mother} #{d_other} PHSP;
    Enddecay

    Decay #{d_mother}
    1.0 #{final_state} #{model};
    Enddecay

    Decay #{d_other}
    1.0 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.0 gamma gamma PHSP;
    Enddecay

    End
  DECAYCARD
end

# --- 9 Cabibbo-suppressed D0 signal modes ---
d0_signal_modes = {
  "pip_pim_pi0"        => ["pi+ pi- pi0",             "BODY3"],  # three-body
  "pip_pim_pi0pi0"     => ["pi+ pi- pi0 pi0",         "PHSP"],
  "pip_pim_etaeta"     => ["pi+ pi- eta eta",         "PHSP"],
  "pi0pi0pi0pi0"       => ["pi0 pi0 pi0 pi0",         "PHSP"],
  "pi0pi0pi0eta"       => ["pi0 pi0 pi0 eta",         "PHSP"],
  "pipipimpimp_pi0"    => ["pi+ pi+ pi- pi- pi0",     "PHSP"],
  "pipipimpimp_eta"    => ["pi+ pi+ pi- pi- eta",     "PHSP"],
  "pipim_pi0pi0pi0"    => ["pi+ pi- pi0 pi0 pi0",     "PHSP"],
  "pipipimpimp_pi0pi0" => ["pi+ pi+ pi- pi- pi0 pi0", "PHSP"]
}

# --- 11 Cabibbo-suppressed D+ signal modes ---
dplus_signal_modes = {
  "pip_pip_pim"             => ["pi+ pi+ pi-",               "BODY3"],  # three-body
  "pip_pi0pi0"              => ["pi+ pi0 pi0",               "BODY3"],  # three-body
  "pip_pip_pim_pi0"         => ["pi+ pi+ pi- pi0",           "PHSP"],
  "pip_pi0pi0pi0"           => ["pi+ pi0 pi0 pi0",           "PHSP"],
  "pip_pip_pip_pim_pim"     => ["pi+ pi+ pi+ pi- pi-",       "PHSP"],
  "pip_pip_pim_pi0pi0"      => ["pi+ pi+ pi- pi0 pi0",       "PHSP"],
  "pip_pip_pim_pi0eta"      => ["pi+ pi+ pi- pi0 eta",       "PHSP"],
  "pip_pi0pi0pi0pi0"        => ["pi+ pi0 pi0 pi0 pi0",       "PHSP"],
  "pip_pi0pi0pi0eta"        => ["pi+ pi0 pi0 pi0 eta",       "PHSP"],
  "pip_pip_pip_pim_pim_pi0" => ["pi+ pi+ pi+ pi- pi- pi0",   "PHSP"],
  "pip_pip_pim_pi0pi0pi0"   => ["pi+ pi+ pi- pi0 pi0 pi0",   "PHSP"]
}

# --- Create one exclusive signal MC sample per signal mode ---
d0_signal_mcs = d0_signal_modes.map do |name, (fs, model)|
  DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3773_D0_#{name}"
    config.related_dataset = data_3773
    config.events          = 100_000
    config.decay_card      = build_cs_signal_card("D0", "anti-D0", fs, model)
    config.cross_section   = :default
  end
end

dplus_signal_mcs = dplus_signal_modes.map do |name, (fs, model)|
  DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3773_Dplus_#{name}"
    config.related_dataset = data_3773
    config.events          = 100_000
    config.decay_card      = build_cs_signal_card("D+", "D-", fs, model)
    config.cross_section   = :default
  end
end

### ---------------------------- Event selection -------------------------- ###
# ============================================================================
# Channel 1 : e+e- -> D0 D0bar  —  hadronic D0bar tag  +  recoil D0 signal
# ============================================================================
alg_d0 = TagAnalysis.new("DT_D0bar_CS_D0")
alg_d0.set_header(["DT_D0bar_CS_D0Alg/DT_D0bar_CS_D0.h"])
      .set_constant({ "ECMS" => [:double, 3.773] })
      .with_decay_card(build_cs_signal_card("D0", "anti-D0", "pi+ pi- pi0", "BODY3"))

# --- Tag side: hadronic D0bar (charm = -1) ---
alg_d0.tag_side(:D0) do |t|
  t.modes :D0toKPi,            # K- pi+
          :D0toKPiPi0,         # K- pi+ pi0
          :D0toKPiPiPi,        # K- pi+ pi+ pi-
          :D0toKsPiPi,         # K_S0 pi+ pi-
          :D0toKsPiPiPi0,      # K_S0 pi+ pi- pi0
          :D0toKsPiPiPiPi,     # K_S0 pi+ pi- pi+ pi-
          :D0toKKPi            # K+ K- pi+
  t.charm -1
  # Opt-in tag-side |DeltaE| window.  The reference analysis uses a
  # mode-dependent window; a single looser union window is declared here (see
  # the tag_deltae_window note).
  t.window :deltaE, min: -0.055, max: 0.040
end

# --- Signal side: recoil D0 (representative D0 -> pi+ pi- pi0) ---
alg_d0.signal_side do |s|
  s.charged(pip: 1, pim: 1)     # one pi+ and one pi-
  s.photons 2                   # at least two photons (for pi0 -> gamma gamma)
  s.min_photon_angle 10.0       # photon-track isolation, degrees
  s.min_photon_energy 0.025     # shower energy floor, GeV
end

# --- Kinematic fit over the derived participants (tag + signal) ---
alg_d0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200                # loose; optimal cut determined in ROOT
end

alg_d0.note(:signal_mode_list,
            "9 Cabibbo-suppressed D0 signal modes form the recoil side: " \
            "pi+pi-pi0, pi+pi-2pi0, pi+pi-2eta, 4pi0, 3pi0eta, 2pi+2pi-pi0, " \
            "2pi+2pi-eta, pi+pi-3pi0, 2pi+2pi-2pi0. The DSL signal_side declares " \
            "only the representative pi+pi-pi0 configuration (1 pi+, 1 pi-, >=2 photons); " \
            "all modes share the same shower criteria (>=10 deg, >=25 MeV).")
alg_d0.note(:tag_deltae_window,
            "Tag DeltaE window is mode dependent: (-55,+40) MeV for tag modes " \
            "containing a pi0 (K-pi+pi0, K_S0pi+pi-pi0) and (-25,+25) MeV otherwise; " \
            "the candidate with the best |DeltaE_tag| is kept per tag mode. A single " \
            "looser union window is declared at BOSS level.")
alg_d0.note(:background_veto,
            "K_S0 background vetoed on the signal side by excluding M(pi+pi-) in " \
            "[0.468, 0.528] GeV/c2 and M(pi0pi0) in [0.428, 0.548] GeV/c2.")
alg_d0.note(:no_signal_side_kinematic_fit,
            "The reference analysis applies NO kinematic fit on the signal side; the " \
            "signal D0 is identified from DeltaE_sig and M_BC_sig with mode-dependent " \
            "DeltaE_sig windows, and yields are extracted in ROOT from 2D unbinned ML " \
            "fits to M_BC^tag vs M_BC^sig (including quantum-correlation corrections).")
alg_d0.note(:track_quality_selection,
            "No explicit |cos(theta)|, |Vz|, Vr or dedicated PID algorithm is imposed " \
            "beyond the charged multiplicities, net charge, photon requirements and " \
            "vetoes listed above.")

alg_d0.apply                                   # no Selection argument
alg_d0.execute_on([data_3773, incMC_3773] + d0_signal_mcs)

# ============================================================================
# Channel 2 : e+e- -> D+ D-  —  hadronic D- tag  +  recoil D+ signal
# ============================================================================
alg_dplus = TagAnalysis.new("DT_Dminus_CS_Dplus")
alg_dplus.set_header(["DT_Dminus_CS_DplusAlg/DT_Dminus_CS_Dplus.h"])
        .set_constant({ "ECMS" => [:double, 3.773] })
        .with_decay_card(build_cs_signal_card("D+", "D-", "pi+ pi+ pi-", "BODY3"))

# --- Tag side: hadronic D- (charm = -1) ---
alg_dplus.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,          # K- pi+ pi-
          :DptoKsPi,           # K_S0 pi-
          :DptoKsPiPiPi,       # K_S0 pi+ pi- pi-
          :DptoKKPi            # K+ K- pi-
  t.charm -1
  # None of the declared D- tag modes contains a pi0 -> narrower (-25,+25) MeV window
  t.window :deltaE, min: -0.025, max: 0.025
end

# --- Signal side: recoil D+ (representative D+ -> 2pi+ pi-) ---
alg_dplus.signal_side do |s|
  s.charged(pip: 2, pim: 1)     # two pi+ and one pi-
  s.require_charge(1)           # net charge +1
  s.min_photon_angle 10.0       # same shower criteria are set; no photons required
  s.min_photon_energy 0.025
end

# --- Kinematic fit over the derived participants (tag + signal) ---
alg_dplus.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_dplus.note(:signal_mode_list,
               "11 Cabibbo-suppressed D+ signal modes form the recoil side: " \
               "2pi+pi-, pi+2pi0, 2pi+pi-pi0, pi+3pi0, 3pi+2pi-, 2pi+pi-2pi0, " \
               "2pi+pi-pi0eta, pi+4pi0, pi+3pi0eta, 3pi+2pi-pi0, 2pi+pi-3pi0. " \
               "The DSL signal_side declares only the representative 2pi+pi- " \
               "configuration (2 pi+, 1 pi-, net charge +1, no required photons); " \
               "all modes share the same shower criteria (>=10 deg, >=25 MeV).")
alg_dplus.note(:tag_deltae_window,
               "Tag DeltaE window (-25,+25) MeV (no declared D- tag mode contains a " \
               "pi0); the candidate with the best |DeltaE_tag| is kept per tag mode.")
alg_dplus.note(:background_veto,
               "K_S0 background vetoed on the signal side by excluding M(pi+pi-) in " \
               "[0.468, 0.528] GeV/c2 and M(pi0pi0) in [0.428, 0.548] GeV/c2; the D+ " \
               "modes use the 4-sigma K_S0 mass windows.")
alg_dplus.note(:no_signal_side_kinematic_fit,
               "The reference analysis applies NO kinematic fit on the signal side; the " \
               "signal D+ is identified from DeltaE_sig and M_BC_sig with mode-dependent " \
               "DeltaE_sig windows, and yields are extracted in ROOT from 2D unbinned ML " \
               "fits to M_BC^tag vs M_BC^sig.")
alg_dplus.note(:track_quality_selection,
               "No explicit |cos(theta)|, |Vz|, Vr or dedicated PID algorithm is imposed " \
               "beyond the charged multiplicities, net charge, photon requirements and " \
               "vetoes listed above.")

alg_dplus.apply                                # no Selection argument
alg_dplus.execute_on([data_3773, incMC_3773] + dplus_signal_mcs)