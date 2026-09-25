# =============================================================================
# BESIII arXiv:1509.01398v2
# Confirmation of the charged charmoniumlike state Zc(3885)-/+ in
# e+e- -> pi± (D D*)∓ with a DOUBLE D TAG.
#
# Data: 1092 pb^-1 at sqrt(s) = 4.23 GeV (sample 703_4230) and
#        826 pb^-1 at sqrt(s) = 4.26 GeV (sample 703_4260) + inclusive MC.
#
# Method: the bachelor pi+ and the D meson pair are fully reconstructed
# ("double D tag", DT). The soft pion from the D*- / D*0 decay is too low in
# momentum to be reconstructed; it is called the "missing pion" and is inferred
# from energy-momentum conservation. Two isospin processes are combined:
#   (a) e+e- -> pi+ D0 D*-    (the D0 D0bar pair is tagged)
#   (b) e+e- -> pi+ D- D*0    (the D- D0 pair is tagged)
# Charge-conjugate processes are included throughout.
#
# -----------------------------------------------------------------------------
# DSL mapping notes (v1 tag vocabulary limitations, recorded per algorithm):
#   * Process (b) pairs two *different* tag species (D- and D0). The DT
#     vocabulary requires both tag_side calls to name the same species
#     ([DSL:tag_dt_species_mismatch]) and `missing` is ST-only
#     ([DSL:tag_dt_missing_unsupported]). The mixed-species DT is therefore
#     expressed as a single D- tag (tag_side(:Dplus), charm -1) with the D0
#     fully reconstructed on the signal side -- one Algorithm per D0 mode
#     (Rule T1: different charged multisets and photon content).
#   * Process (a) is a genuine same-species DT (D0 / anti-D0), so it is written
#     with two tag_side(:D0) calls; the missing pi- of D*- -> anti-D0 pi- cannot
#     enter the DT fit (DT + missing is unsupported) and is recorded in notes.
#   * ECMS is not declared: the analysis spans two centre-of-mass energies and
#     the generated code reads the measured per-run beam energy and boost from
#     MeasuredEcmsSvc (DtagReconstructionBuilder default `beam_energy :db`).
#     The declared energy is only the MC fallback.
# =============================================================================

### Dataset description ###
data_4230  = DatasetManager.real_data.find("703_4230")      # sqrt(s) = 4.226 GeV, 1092 pb^-1
data_4260  = DatasetManager.real_data.find("703_4260")      # sqrt(s) = 4.258 GeV,  826 pb^-1
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")   # generic inclusive MC (KKMC + EvtGen + Lundcharm)
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

data_points = [data_4230, data_4260]
incMCs      = [incMC_4230, incMC_4260]

# -----------------------------------------------------------------------------
# Decay cards (EvtGen syntax).
# The D mesons decay inclusively over the reconstructed modes; the fractional
# weights are the PDG branching fractions normalised over the modes used in the
# analysis (the paper's own MC uses the full PDG values).
# The top mother is psi(4260) (the BESIII KKMC convention for e+e- annihilation
# in this energy region).
# -----------------------------------------------------------------------------

# Signal card (a): e+e- -> pi+ Zc(3885)- , Zc(3885)- -> D0 anti-D*- ,
# anti-D*- -> anti-D0 pi-  (the pi- is the undetected "missing pion").
decay_card_zc_a = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ Zc(3885)-                      PHSP;
  Enddecay

  Decay Zc(3885)-
  1.0000 D0 anti-D*-                        PHSP;
  Enddecay

  Decay anti-D*-
  1.0000 anti-D0 pi-                        PHSP;
  Enddecay

  Decay D0
  0.1250 K- pi+                             PHSP;
  0.4490 K- pi+ pi0                         PHSP;
  0.2580 K- pi+ pi+ pi-                     PHSP;
  0.1680 K- pi+ pi+ pi- pi0                 PHSP;
  Enddecay

  Decay anti-D0
  0.1250 K+ pi-                             PHSP;
  0.4490 K+ pi- pi0                         PHSP;
  0.2580 K+ pi- pi+ pi-                     PHSP;
  0.1680 K+ pi- pi+ pi- pi0                 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                        PHSP;
  Enddecay

  End
DECAYCARD

# Signal card (b): e+e- -> pi+ Zc(3885)- , Zc(3885)- -> D- D*0 ,
# D*0 -> D0 pi0  (the pi0 is the undetected "missing pion").
decay_card_zc_b = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ Zc(3885)-                      PHSP;
  Enddecay

  Decay Zc(3885)-
  1.0000 D- D*0                             PHSP;
  Enddecay

  Decay D*0
  1.0000 D0 pi0                             PHSP;
  Enddecay

  Decay D-
  0.3300 K+ pi- pi-                         PHSP;
  0.2100 K+ pi- pi- pi0                     PHSP;
  0.0540 K_S0 pi-                           PHSP;
  0.1480 K_S0 pi- pi0                       PHSP;
  0.1130 K_S0 pi+ pi- pi-                   PHSP;
  0.0360 K+ K- pi-                          PHSP;
  0.1090 pi- pi0                            PHSP;
  Enddecay

  Decay D0
  0.1250 K- pi+                             PHSP;
  0.4490 K- pi+ pi0                         PHSP;
  0.2580 K- pi+ pi+ pi-                     PHSP;
  0.1680 K- pi+ pi+ pi- pi0                 PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-                            PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                        PHSP;
  Enddecay

  End
DECAYCARD

# Phase-space card (a): e+e- -> pi+ (D D*)- with no Zc resonance.
decay_card_phsp_a = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ D0 anti-D*-                    PHSP;
  Enddecay

  Decay anti-D*-
  1.0000 anti-D0 pi-                        PHSP;
  Enddecay

  Decay D0
  0.1250 K- pi+                             PHSP;
  0.4490 K- pi+ pi0                         PHSP;
  0.2580 K- pi+ pi+ pi-                     PHSP;
  0.1680 K- pi+ pi+ pi- pi0                 PHSP;
  Enddecay

  Decay anti-D0
  0.1250 K+ pi-                             PHSP;
  0.4490 K+ pi- pi0                         PHSP;
  0.2580 K+ pi- pi+ pi-                     PHSP;
  0.1680 K+ pi- pi+ pi- pi0                 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                        PHSP;
  Enddecay

  End
DECAYCARD

# Phase-space card (b): e+e- -> pi+ (D D*)- with no Zc resonance.
decay_card_phsp_b = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ D- D*0                         PHSP;
  Enddecay

  Decay D*0
  1.0000 D0 pi0                             PHSP;
  Enddecay

  Decay D-
  0.3300 K+ pi- pi-                         PHSP;
  0.2100 K+ pi- pi- pi0                     PHSP;
  0.0540 K_S0 pi-                           PHSP;
  0.1480 K_S0 pi- pi0                       PHSP;
  0.1130 K_S0 pi+ pi- pi-                   PHSP;
  0.0360 K+ K- pi-                          PHSP;
  0.1090 pi- pi0                            PHSP;
  Enddecay

  Decay D0
  0.1250 K- pi+                             PHSP;
  0.4490 K- pi+ pi0                         PHSP;
  0.2580 K- pi+ pi+ pi-                     PHSP;
  0.1680 K- pi+ pi+ pi- pi0                 PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-                            PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                        PHSP;
  Enddecay

  End
DECAYCARD

# Background card: e+e- -> D1(2420) anti-D0 with D1(2420) -> pi- D*+ ,
# D*+ -> D0 pi+ . This is the only process able to produce a peak at the
# D Dbar* threshold; its |cos(theta_piD)| asymmetry in generic MC (0.43 +- 0.01)
# is compared with data (0.11 +- 0.07) to bound its contribution.
decay_card_bkg_d1 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D1(2420)0 anti-D0                  PHSP;
  Enddecay

  Decay D1(2420)0
  1.0000 pi- D*+                            PHSP;
  Enddecay

  Decay D*+
  1.0000 D0 pi+                             PHSP;
  Enddecay

  Decay D0
  1.0000 K- pi+                             PHSP;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi-                             PHSP;
  Enddecay

  End
DECAYCARD

exMC_zc_a = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Zc3885_pi_D0_Dstar_sig"
  config.events        = 200_000
  config.decay_card    = decay_card_zc_a
  config.cross_section = :default
end

exMC_zc_b = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Zc3885_pi_Dm_Dstar0_sig"
  config.events        = 200_000
  config.decay_card    = decay_card_zc_b
  config.cross_section = :default
end

exMC_phsp_a = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Zc3885_pi_D0_Dstar_phsp"
  config.events        = 200_000
  config.decay_card    = decay_card_phsp_a
  config.cross_section = :default
end

exMC_phsp_b = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Zc3885_pi_Dm_Dstar0_phsp"
  config.events        = 200_000
  config.decay_card    = decay_card_phsp_b
  config.cross_section = :default
end

exMC_bkg_d1 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "D1_2420_D0bar_bkg"
  config.events        = 200_000
  config.decay_card    = decay_card_bkg_d1
  config.cross_section = :default
end

### ---------------------------------------------------------------------------
### Process (a): e+e- -> pi+ D0 anti-D*-  (pi+ D0 D0bar tagged, DT)
### ---------------------------------------------------------------------------
alg_d0d0 = TagAnalysis.new("Zc3885PiD0D0barTag")
alg_d0d0.set_header(["Zc3885PiD0D0barTagAlg/Zc3885PiD0D0barTag.h"])
        # Two centre-of-mass energies (4.23 and 4.26 GeV): no fixed ECMS
        # constant, the measured per-run beam energy/boost is read from the DB.
        .with_decay_card(decay_card_zc_a)
        .note(:track_selection,
              "Good charged tracks: |cos(theta)| < 0.93 and the point of closest " \
              "approach to the interaction point within 10 cm along the beam " \
              "direction and 1 cm in the plane perpendicular to it.")
        .note(:pid,
              "dE/dx and TOF information are combined into Prob(K) and Prob(pi); a " \
              "track is a K (pi) when Prob(K) > Prob(pi) (Prob(pi) > Prob(K)). " \
              "Tracks used in K_S0 reconstruction are exempt from these " \
              "requirements.")
        .note(:photon_selection,
              "EMC clusters: E > 25 MeV in the barrel (|cos(theta)| < 0.8) and " \
              "E > 50 MeV in the endcap (0.84 < |cos(theta)| < 0.92); angle to the " \
              "nearest charged track > 20 deg; EMC timing used to suppress " \
              "electronic noise and unrelated energy deposits.")
        .note(:pi0_selection,
              "pi0 candidates from photon pairs with 0.115 < M(gamma gamma) < " \
              "0.150 GeV/c^2; a 1C kinematic fit constrains M(gamma gamma) to the " \
              "PDG pi0 mass. The modes carrying a pi0 enter the DT through the " \
              "DTagAlg channel names (:D0toKPiPi0, :D0toKPiPiPiPi0).")
        .note(:ks0_selection,
              "K_S0 candidates from oppositely charged track pairs with " \
              "|cos(theta)| < 0.93 and distance to the IP along the beam within " \
              "20 cm; vertex fit to a common decay vertex with " \
              "0.487 < M(pi+ pi-) < 0.511 GeV/c^2, plus a secondary-vertex " \
              "constraint between the production and decay vertices.")
        .note(:dt_candidate_selection,
              "If there is more than one candidate per possible DT mode, the " \
              "candidate with the minimum Delta Mhat is kept, where " \
              "Delta Mhat = [M(D) + M(Dbar)]/2 - [M_PDG(D) + M_PDG(Dbar)]/2. The " \
              "D0 D0bar signal region is -20 < Delta Mhat < 15 MeV/c^2 together " \
              "with |Delta M| = |M(D) - M(Dbar)| < 40 MeV/c^2. Both observables are " \
              "stored and windowed downstream in ROOT. The DTagAlg selector cuts " \
              "stay disabled; the DT ranking across tag modes uses findDTag's own " \
              "ordering.")
        .note(:missing_pi,
              "The soft pi- from anti-D*- -> anti-D0 pi- is not reconstructed; it " \
              "is the 'missing pion' inferred from energy-momentum conservation. " \
              "Because a double tag leaves nothing unreconstructed, the v1 tag " \
              "vocabulary cannot declare `missing` together with two tag_side " \
              "calls ([DSL:tag_dt_missing_unsupported]); the D*- mass constraint " \
              "(missing pi- + anti-D0 -> M_PDG(D*-)) of the paper's fit is therefore " \
              "not emitted here and is applied on the stored four-momenta in ROOT.")
        .note(:bachelor_pion,
              "The bachelor pi+ is required as at least one good charged track not " \
              "among the decay products of the two D candidates (declared on the " \
              "signal side with at_least: true). The charge-conjugate configuration " \
              "(bachelor pi-) is covered by the corresponding charge-conjugate " \
              "sample and selection at the ROOT stage.")
        .note(:kinematic_fit,
              "The paper performs a 4C fit imposing momentum and energy " \
              "conservation, constraining the invariant mass of both D candidates " \
              "to M_PDG(D) and constraining the missing pi plus the corresponding D " \
              "candidate to M_PDG(D*) -- 7 constraints in total, of which four " \
              "remain free because the missing pion three-momentum is unknown. " \
              "chi2_4C < 100 is required and the candidate with the minimum " \
              "chi2_4C is kept; the loose BOSS value 200 is used here and the " \
              "published cut is applied in ROOT (Rule T3).")
        .note(:background_veto,
              "e+e- -> D* D* background is suppressed by M(pi+ D0) >= 2.03 GeV/c^2 " \
              "(process a) and M(pi+ D-) > 2.08 GeV/c^2 (process b). The " \
              "reconstructed D pi recoil mass, " \
              "M_recoil(D pi)^2 c^4 = (E_cm - E_D - E_pi)^2 - |p_cm - p_D - p_pi|^2 " \
              "c^2, must satisfy |M_recoil(D pi) - M_PDG(D*)| < 30 MeV/c^2. Both are " \
              "computed from the fitted four-momenta and applied in ROOT.")
        .note(:angular_distribution,
              "The |cos(theta_pi)| distribution, where theta_pi is the polar angle " \
              "of the bachelor pi+ relative to the beam in the CMS frame, is flat " \
              "as expected for J^P = 1+ (chi2/NDF = 16.5/9) and disfavours " \
              "J^P = 0- (sin^2 theta_pi, 103.1/9) and J^P = 1- " \
              "(1 + cos^2 theta_pi, 106.3/9). The |cos(theta_piD)| asymmetry " \
              "A = (n_>0.5 - n_<0.5)/(n_>0.5 + n_<0.5) is used in ROOT.")
        .note(:results,
              "Simultaneous unbinned maximum-likelihood fit of the M(D Dbar*) " \
              "distributions of both isospin processes at both energies with a " \
              "mass-dependent Breit-Wigner: statistical significance > 10 sigma; " \
              "M_pole = (3881.7 +- 1.6 (stat.) +- 1.6 (syst.)) MeV/c^2 and " \
              "Gamma_pole = (26.6 +- 2.0 (stat.) +- 2.1 (syst.)) MeV. " \
              "sigma x Br = (141.6 +- 7.9 +- 12.3) pb at 4.23 GeV and " \
              "(108.4 +- 6.9 +- 8.8) pb at 4.26 GeV. The mass resolution is " \
              "1.1 +- 0.1 MeV/c^2 for the pi+ D0 D0bar-tagged process.")

# Tag side 1: D0 (the directly produced neutral D)
alg_d0d0.tag_side(:D0) do |t|
  t.modes :D0toKPi,           # K- pi+
          :D0toKPiPi0,        # K- pi+ pi0
          :D0toKPiPiPi,       # K- pi+ pi+ pi-
          :D0toKPiPiPiPi0     # K- pi+ pi+ pi- pi0
  t.charm 1
end

# Tag side 2: anti-D0 (recoiling D from the D*- decay)
alg_d0d0.tag_side(:D0) do |t|
  t.modes :D0toKPi,
          :D0toKPiPi0,
          :D0toKPiPiPi,
          :D0toKPiPiPiPi0
  t.charm -1
end

# Signal side: the bachelor pi+ not used by either tag.
alg_d0d0.signal_side do |s|
  s.charged(pip: 1, at_least: true)   # at least one extra good charged track
  s.require_charge 1                  # +1 (bachelor pi+)
  s.min_photon_angle 20.0
end

alg_d0d0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)   # m(D0) constraint
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)   # m(anti-D0) constraint
  f.chi2_cut 200                                                 # loose BOSS cut; chi2 < 100 in ROOT
end

alg_d0d0.apply
alg_d0d0.execute_on(data_points + incMCs + exMC_zc_a + exMC_phsp_a + exMC_bkg_d1)

### ---------------------------------------------------------------------------
### Process (b): e+e- -> pi+ D- D*0  (pi+ D- D0 tagged)
### Mixed-species DT (D- and D0) -> one D- tag side + the D0 reconstructed on
### the signal side. One Algorithm per D0 decay mode (Rule T1).
### ---------------------------------------------------------------------------

# D- tag modes shared by all four algorithms below.
dm_tag_modes = [:DptoKPiPi,      # D- -> K+ pi- pi-
                :DptoKPiPiPi0,   # D- -> K+ pi- pi- pi0
                :DptoKsPi,       # D- -> K_S0 pi-
                :DptoKsPiPi0,    # D- -> K_S0 pi- pi0
                :DptoKsPiPiPi,   # D- -> K_S0 pi+ pi- pi-
                :DptoKKPi]       # D- -> K+ K- pi-

# --- (b1) D0 -> K- pi+ -------------------------------------------------------
alg_b_kpi = TagAnalysis.new("Zc3885PiDmD0TagKPi")
alg_b_kpi.set_header(["Zc3885PiDmD0TagKPiAlg/Zc3885PiDmD0TagKPi.h"])
         .with_decay_card(decay_card_zc_b)
         .note(:mixed_species_dt,
               "The D- D0 pair is a mixed-species double tag. The v1 tag " \
               "vocabulary pairs two tag_side calls of the same species only " \
               "([DSL:tag_dt_species_mismatch]), so the D- is taken from the " \
               "pre-stored tag collection (tag_side(:Dplus), charm -1) while the " \
               "D0 is fully reconstructed on the signal side. One Algorithm per D0 " \
               "mode is needed because the D0 modes carry different charged " \
               "multisets and photon content (Rule T1).")
         .note(:missing_pi,
               "The soft pi0 from D*0 -> D0 pi0 (and from the isospin feed " \
               "D*- -> D- pi0) is not reconstructed; it is the 'missing pion' " \
               "inferred from energy-momentum conservation. It is declared via " \
               "s.missing so that a four-momentum-conserving fit is performed; the " \
               "paper's additional constraint (missing pi0 + D0 -> M_PDG(D*0)) " \
               "cannot be expressed as a resonance over the derived participants " \
               "and is applied on the stored four-momenta in ROOT.")
         .note(:dt_candidate_selection,
               "The D- D0 signal region is -17 < Delta Mhat < 14 MeV/c^2 together " \
               "with |Delta M| = |M(D-) - M(D0)| < 35 MeV/c^2; multiple candidates " \
               "are ranked by the minimum Delta Mhat. Both observables are stored " \
               "and windowed downstream in ROOT.")
         .note(:isospin_feed,
               "Events from the isospin partner channel e+e- -> pi+ D0 D*- with " \
               "D*- -> D- pi0 can satisfy the same requirements with a different " \
               "efficiency and mass resolution; they are treated as signal and " \
               "combined with the pi+ D- D0-tagged process. This is why the mass " \
               "resolution is described by a sum of two Crystal Ball functions " \
               "(2.2 +- 0.1 MeV/c^2).")
         .note(:track_pid_photon,
               "Good tracks |cos(theta)| < 0.93, PCA within 10 cm along the beam " \
               "and 1 cm transverse; K/pi separation from combined dE/dx and TOF " \
               "likelihoods (Prob(K) > Prob(pi) for kaons, Prob(pi) > Prob(K) for " \
               "pions), K_S0 daughters exempt; photons E > 25 MeV barrel " \
               "(|cos(theta)| < 0.8) / E > 50 MeV endcap (0.84 < |cos(theta)| < " \
               "0.92) and more than 20 deg from any charged track; pi0 from " \
               "0.115 < M(gamma gamma) < 0.150 GeV/c^2 with a 1C mass-constrained " \
               "fit; K_S0 with a vertex fit and 0.487 < M(pi+ pi-) < " \
               "0.511 GeV/c^2 plus a secondary-vertex constraint.")
         .note(:background_veto,
               "M(pi+ D-) > 2.08 GeV/c^2 suppresses e+e- -> D* D*; the D pi recoil " \
               "mass must satisfy |M_recoil(D pi) - M_PDG(D*)| < 30 MeV/c^2. Both " \
               "are applied in ROOT on the fitted four-momenta.")
         .note(:results,
               "M_pole = (3881.7 +- 1.6 +- 1.6) MeV/c^2, " \
               "Gamma_pole = (26.6 +- 2.0 +- 2.1) MeV; sigma x Br = " \
               "(141.6 +- 7.9 +- 12.3) pb at 4.23 GeV and (108.4 +- 6.9 +- 8.8) pb " \
               "at 4.26 GeV (simultaneous fit of both isospin processes).")

alg_b_kpi.tag_side(:Dplus) do |t|
  t.modes(*dm_tag_modes)
  t.charm -1                       # tagged D-
end

# D0 -> K- pi+ plus the bachelor pi+
alg_b_kpi.signal_side do |s|
  s.charged(km: 1, pip: 2, at_least: true)
  s.require_charge 1               # -1 (K-) + 2 (pi+) = +1
  s.min_photon_angle 20.0
  s.missing :pi0                   # missing soft pi0 from D*0 -> D0 pi0
end

alg_b_kpi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)  # m(D0) constraint
  f.chi2_cut 200                   # loose BOSS cut; chi2 < 100 in ROOT
end

alg_b_kpi.apply
alg_b_kpi.execute_on(data_points + incMCs + exMC_zc_b + exMC_phsp_b + exMC_bkg_d1)

# --- (b2) D0 -> K- pi+ pi0 ---------------------------------------------------
alg_b_kpipi0 = TagAnalysis.new("Zc3885PiDmD0TagKPiPi0")
alg_b_kpipi0.set_header(["Zc3885PiDmD0TagKPiPi0Alg/Zc3885PiDmD0TagKPiPi0.h"])
            .with_decay_card(decay_card_zc_b)
            .note(:mixed_species_dt,
                  "Mixed-species D- D0 double tag expressed as a D- tag side plus " \
                  "the fully reconstructed D0 on the signal side (see " \
                  "Zc3885PiDmD0TagKPi for the full rationale).")
            .note(:missing_pi,
                  "The soft pi0 from D*0 -> D0 pi0 is the missing pion, inferred " \
                  "from energy-momentum conservation (s.missing :pi0). The " \
                  "D*0 mass constraint of the paper's fit is applied in ROOT.")
            .note(:pi0_on_signal_side,
                  "The D0 -> K- pi+ pi0 candidate is completed by the pi0 built " \
                  "from the two signal-side photons; the 1C mass-constrained fit " \
                  "to the PDG pi0 mass uses 0.115 < M(gamma gamma) < " \
                  "0.150 GeV/c^2 with both photons in the barrel-favoured " \
                  "(|cos(theta)| < 0.8) / endcap (0.84 < |cos(theta)| < 0.92) " \
                  "acceptance and more than 20 deg from any charged track.")
            .note(:background_veto,
                  "M(pi+ D-) > 2.08 GeV/c^2 and " \
                  "|M_recoil(D pi) - M_PDG(D*)| < 30 MeV/c^2, applied in ROOT.")

alg_b_kpipi0.tag_side(:Dplus) do |t|
  t.modes(*dm_tag_modes)
  t.charm -1
end

# D0 -> K- pi+ pi0 plus the bachelor pi+; the pi0 is the two signal photons.
alg_b_kpipi0.signal_side do |s|
  s.charged(km: 1, pip: 2, at_least: true)
  s.require_charge 1
  s.photons 2                      # pi0 -> gamma gamma of the D0
  s.min_photon_angle 20.0
  s.min_photon_energy 0.025
  s.missing :pi0
end

alg_b_kpipi0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_b_kpipi0.apply
alg_b_kpipi0.execute_on(data_points + incMCs + exMC_zc_b + exMC_phsp_b + exMC_bkg_d1)

# --- (b3) D0 -> K- pi+ pi+ pi- ----------------------------------------------
alg_b_kpipipi = TagAnalysis.new("Zc3885PiDmD0TagKPiPiPi")
alg_b_kpipipi.set_header(["Zc3885PiDmD0TagKPiPiPiAlg/Zc3885PiDmD0TagKPiPiPi.h"])
             .with_decay_card(decay_card_zc_b)
             .note(:mixed_species_dt,
                   "Mixed-species D- D0 double tag expressed as a D- tag side plus " \
                   "the fully reconstructed D0 on the signal side (see " \
                   "Zc3885PiDmD0TagKPi for the full rationale).")
             .note(:missing_pi,
                   "The soft pi0 from D*0 -> D0 pi0 is the missing pion, inferred " \
                   "from energy-momentum conservation (s.missing :pi0). The " \
                   "D*0 mass constraint of the paper's fit is applied in ROOT.")
             .note(:d0_mass_constraint_omitted,
                   "The D0 -> K- pi+ pi+ pi- mass constraint is not emitted: the " \
                   "signal side carries three pi+ (two from the D0 and one bachelor) " \
                   "so the invariant_mass_of(:km, :pip, :pip, :pim) sub-system is " \
                   "ambiguous under the derived-participant rules. The D0 mass " \
                   "window is applied in ROOT.")
             .note(:background_veto,
                   "M(pi+ D-) > 2.08 GeV/c^2 and " \
                   "|M_recoil(D pi) - M_PDG(D*)| < 30 MeV/c^2, applied in ROOT.")

alg_b_kpipipi.tag_side(:Dplus) do |t|
  t.modes(*dm_tag_modes)
  t.charm -1
end

# D0 -> K- pi+ pi+ pi- plus the bachelor pi+
alg_b_kpipipi.signal_side do |s|
  s.charged(km: 1, pip: 3, pim: 1, at_least: true)
  s.require_charge 1               # -1 + 3 - 1 = +1
  s.min_photon_angle 20.0
  s.missing :pi0
end

alg_b_kpipipi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_b_kpipipi.apply
alg_b_kpipipi.execute_on(data_points + incMCs + exMC_zc_b + exMC_phsp_b + exMC_bkg_d1)

# --- (b4) D0 -> K- pi+ pi+ pi- pi0 ------------------------------------------
alg_b_kpipipipi0 = TagAnalysis.new("Zc3885PiDmD0TagKPiPiPiPi0")
alg_b_kpipipipi0.set_header(["Zc3885PiDmD0TagKPiPiPiPi0Alg/Zc3885PiDmD0TagKPiPiPiPi0.h"])
                .with_decay_card(decay_card_zc_b)
                .note(:mixed_species_dt,
                      "Mixed-species D- D0 double tag expressed as a D- tag side " \
                      "plus the fully reconstructed D0 on the signal side (see " \
                      "Zc3885PiDmD0TagKPi for the full rationale).")
                .note(:missing_pi,
                      "The soft pi0 from D*0 -> D0 pi0 is the missing pion " \
                      "(s.missing :pi0); the D*0 mass constraint of the paper's fit " \
                      "is applied in ROOT. This mode also carries a real pi0 from " \
                      "the D0 decay, reconstructed from the two signal photons.")
                .note(:d0_mass_constraint_omitted,
                      "The D0 mass constraint is not emitted (three pi+ on the " \
                      "signal side make the invariant_mass_of sub-system ambiguous " \
                      "under the derived-participant rules); the D0 mass window is " \
                      "applied in ROOT.")
                .note(:background_veto,
                      "M(pi+ D-) > 2.08 GeV/c^2 and " \
                      "|M_recoil(D pi) - M_PDG(D*)| < 30 MeV/c^2, applied in ROOT.")

alg_b_kpipipipi0.tag_side(:Dplus) do |t|
  t.modes(*dm_tag_modes)
  t.charm -1
end

# D0 -> K- pi+ pi+ pi- pi0 plus the bachelor pi+
alg_b_kpipipipi0.signal_side do |s|
  s.charged(km: 1, pip: 3, pim: 1, at_least: true)
  s.require_charge 1
  s.photons 2                      # pi0 -> gamma gamma of the D0
  s.min_photon_angle 20.0
  s.min_photon_energy 0.025
  s.missing :pi0
end

alg_b_kpipipipi0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_b_kpipipipi0.apply
alg_b_kpipipipi0.execute_on(data_points + incMCs + exMC_zc_b + exMC_phsp_b + exMC_bkg_d1)
