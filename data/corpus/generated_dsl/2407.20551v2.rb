# =============================================================================
#  BOSS part — ψ(3770) semileptonic double-tag analysis
#    tag side  : D0bar (6 hadronic modes) / D- (6 hadronic modes)   [DTag]
#    signal    : D0 -> b1(1235)- e+ nu_e   and   D+ -> b1(1235)0 e+ nu_e
#               (b1 -> omega pi, omega -> pi+pi-pi0, pi0 -> gamma gamma)
#  Tag-based analysis  =>  TagAnalysis layer (no Selection object, no explicit
#  select_track / select_photon / pid: the tag carries its own selected, PID'ed
#  tracks and showers, and the signal side is built from what the tag did not use)
# =============================================================================

### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data (≈7.9 fb⁻¹)
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

# -----------------------------------------------------------------------------
# Decay cards (EvtGen) — one per signal mode
# -----------------------------------------------------------------------------
# ψ(3770) -> D0 D0bar, D0 -> b1(1235)- e+ nu_e, b1- -> omega pi-, omega -> pi+pi-pi0,
# pi0 -> gamma gamma, D0bar -> K+ pi-  (the hadronic tag channel)
decay_card_D0_b1_enu = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0            PHSP;
    Enddecay

    Decay D0
    1.0000 b_1(1235)- e+ nu_e    PHSP;
    Enddecay

    Decay b_1(1235)-
    1.0000 omega pi-             PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0           OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma           PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi-                PHSP;
    Enddecay

    End
DECAYCARD

# ψ(3770) -> D+ D-, D+ -> b1(1235)0 e+ nu_e, b10 -> omega pi0, omega -> pi+pi-pi0,
# pi0 -> gamma gamma, D- -> K+ pi- pi-  (the hadronic tag channel)
decay_card_Dp_b1_enu = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-                 PHSP;
    Enddecay

    Decay D+
    1.0000 b_1(1235)0 e+ nu_e    PHSP;
    Enddecay

    Decay b_1(1235)0
    1.0000 omega pi0             PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0           OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma           PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi-            PHSP;
    Enddecay

    End
DECAYCARD

# -----------------------------------------------------------------------------
# Exclusive MC (500k events per signal mode)
# -----------------------------------------------------------------------------
exMC_D0_b1enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0_b1enu_D0barKPi"   # ψ(3770) -> (D0 -> b1- e+ nu_e) (D0bar -> K+ pi-)
  config.related_dataset = psi3770_data
  config.events          = 500000
  config.decay_card      = decay_card_D0_b1_enu
  config.cross_section   = :default
end

exMC_Dp_b1enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dp_b1enu_DmKPiPi"    # ψ(3770) -> (D+ -> b10 e+ nu_e) (D- -> K+ pi- pi-)
  config.related_dataset = psi3770_data
  config.events          = 500000
  config.decay_card      = decay_card_Dp_b1_enu
  config.cross_section   = :default
end

# =============================================================================
# Algorithm I — D0bar hadronic tag + D0 -> b1(1235)- e+ nu_e (missing nu_e)
# =============================================================================
alg_name_D0 = "D0barTagB1Eenu"
alg_D0 = TagAnalysis.new(alg_name_D0)
alg_D0.set_header(["#{alg_name_D0}Alg/#{alg_name_D0}.h"])
      .set_constant({"ECMS" => [:double, 3.773]})          # √s = 3.773 GeV (fit reference)
      .set_alias({"std::vector<double>" => "Vdouble"})
      .with_decay_card(decay_card_D0_b1_enu)
      .note(:tag_deltaE_windows,
            "per-mode ΔE windows (Table 1) required for the D0bar tag modes; the tag layer " \
            "offers only a single opt-in window per observable per side and no numeric " \
            "per-mode values were provided — the per-mode ΔE selection is applied on the " \
            "stored tag candidates")
      .note(:tag_candidate_selection,
            "when several tag candidates of the same mode and charge survive in one event, " \
            "only the one with the minimum |ΔE| is kept; this per-candidate ranking is not " \
            "expressible in the tag layer and is performed on the stored tag quantities")
      .note(:electron_pid_criteria,
            "signal-side positron identified with CL_e > 0.001, " \
            "CL_e/(CL_e+CL_pi+CL_K) > 0.8 and E/p > 0.8; the tag layer applies fixed " \
            "SimplePIDSvc lepton criteria for the signal-side electron which are not " \
            "DSL-tunable")
      .note(:background_veto,
            "K_S0 veto |M(pi+pi-) - m(K_S0)| > 0.008 GeV/c^2; a0(980) veto " \
            "M(pi+pi-pi0) > 0.6 GeV/c^2; M(b1 e+) < 1.82 GeV/c^2")
      .note(:omega_mass_window,
            "omega signal region M(pi+pi-pi0) = 0.757-0.807 GeV/c^2, sidebands " \
            "0.697-0.742 and 0.822-0.867 GeV/c^2 used for background estimation; the " \
            "3-particle mass window is not encodable as a fit constraint over the derived " \
            "tag-fit participants")
      .note(:extra_shower_veto,
            "maximum energy of any extra photon < 0.30 GeV, no extra pi0 allowed, and " \
            "cos(theta(gamma, missing)) < 0.3 for the D0 signal")

# --- tag side: D0bar, six hadronic modes, pinned by charm = -1 -----------------
alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKsPiPi,
          :D0toKPiPi0Pi0, :D0toKPiPiPiPi0
  t.charm  -1                                  # tag the D0bar (K+ pi- ...) combination
  t.window :mBC, min: 1.859, max: 1.873        # D0bar MBC window (GeV/c^2)
end

# --- signal side: D0 -> b1- e+ nu_e (2 photons, 1 e+, 1 pi+, 2 pi-, net charge 0,
#     one missing nu_e, photon-track angle > 10 deg) -----------------------------
alg_D0.signal_side do |s|
  s.photons          2                         # pi0 -> gamma gamma
  s.charged(ep: 1, pip: 1, pim: 2)             # e+ ; omega -> pi+ pi- ; b1- -> omega pi-
  s.require_charge   0                         # net charge of the signal side
  s.missing          :nu_e                     # massless missing neutrino
  s.min_photon_angle 10.0                      # deg
end

# --- kinematic fit: 4-momentum conservation + gamma gamma mass constraint to pi0,
#     loose chi2 cut (tight cut is determined later in ROOT) ---------------------
alg_D0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_D0.apply                                     # no Selection argument for tag analyses
alg_D0.execute_on([psi3770_data, psi3770_incMC, exMC_D0_b1enu])

# =============================================================================
# Algorithm II — D- hadronic tag + D+ -> b1(1235)0 e+ nu_e (missing nu_e)
# =============================================================================
alg_name_Dm = "DmTagB1Eenu"
alg_Dm = TagAnalysis.new(alg_name_Dm)
alg_Dm.set_header(["#{alg_name_Dm}Alg/#{alg_name_Dm}.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      .with_decay_card(decay_card_Dp_b1_enu)
      .note(:tag_deltaE_windows,
            "per-mode ΔE windows (Table 1) required for the D- tag modes; only a single " \
            "opt-in window per observable per side exists and no numeric per-mode values " \
            "were provided — the per-mode ΔE selection is applied on the stored tag candidates")
      .note(:tag_candidate_selection,
            "when several tag candidates of the same mode and charge survive in one event, " \
            "only the one with the minimum |ΔE| is kept; this per-candidate ranking is not " \
            "expressible in the tag layer")
      .note(:electron_pid_criteria,
            "signal-side positron identified with CL_e > 0.001, " \
            "CL_e/(CL_e+CL_pi+CL_K) > 0.8 and E/p > 0.8; the tag layer applies fixed " \
            "SimplePIDSvc lepton criteria for the signal-side electron which are not " \
            "DSL-tunable")
      .note(:background_veto,
            "K_S0 veto |M(pi+pi-) - m(K_S0)| > 0.008 GeV/c^2; a0(980) veto " \
            "M(pi+pi-pi0) > 0.6 GeV/c^2; M(b1 e+) < 1.82 GeV/c^2")
      .note(:omega_mass_window,
            "omega signal region M(pi+pi-pi0) = 0.757-0.807 GeV/c^2, sidebands " \
            "0.697-0.742 and 0.822-0.867 GeV/c^2 used for background estimation; not " \
            "encodable as a fit constraint over the derived tag-fit participants")
      .note(:extra_shower_veto,
            "maximum energy of any extra photon < 0.30 GeV, no extra pi0 allowed, and " \
            "cos(theta(gamma, missing)) < 0.4 for the D+ signal")

# --- tag side: D-, six hadronic modes (charge fixed by the species alias) ------
alg_Dm.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0,
          :DptoKsPiPiPi, :DptoKKPi
  t.window :mBC, min: 1.863, max: 1.877        # D- MBC window (GeV/c^2)
end

# --- signal side: D+ -> b10 e+ nu_e (4 photons from pi0pi0, 1 e+, 1 pi+, 1 pi-,
#     net charge +1, one missing nu_e, photon-track angle > 10 deg) --------------
alg_Dm.signal_side do |s|
  s.photons          4                         # omega -> pi+pi-pi0 and b10 -> omega pi0, both pi0 -> gamma gamma
  s.charged(ep: 1, pip: 1, pim: 1)             # e+ ; omega -> pi+ pi-
  s.require_charge   1                         # net charge of the signal side
  s.missing          :nu_e                     # massless missing neutrino
  s.min_photon_angle 10.0                      # deg
end

# --- kinematic fit: 4-momentum conservation + gamma gamma mass constraint to pi0,
#     loose chi2 cut ------------------------------------------------------------
alg_Dm.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_Dm.apply
alg_Dm.execute_on([psi3770_data, psi3770_incMC, exMC_Dp_b1enu])