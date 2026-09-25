# BESIII arXiv:1507.02404v3
# Study of e+e- -> (D* Dbar*)^0 pi^0 and observation of the neutral charmoniumlike
# state Z_c(4025)^0 decaying to (D* Dbar*)^0, using data at sqrt(s) = 4.23 GeV
# (1092 pb^-1) and 4.26 GeV (826 pb^-1).
#
# Method: a PARTIAL RECONSTRUCTION is applied. A D and a Dbar originating from the
# D* -> D pi / D gamma decays are detected, together with the pi0 from the primary
# production (the "bachelor" pi0). The soft pi0/gamma from the D* decay is NOT used
# in the recoil-mass computation; the D* is inferred from the recoil four-momentum.
# The signal is extracted from the recoil-mass spectra RM(D pi0) and RM(Dbar pi0).
#
# The D+ is required to decay into K- pi+ pi+, while the D0 is required to decay into
# K- pi+, K- pi+ pi0 and K- pi+ pi+ pi-. Each D decay mode has a different charged-track
# multiplicity and reconstruction hypothesis, so each gets its own Algorithm (Rule T1).

### Dataset description ###
data_4230  = DatasetManager.real_data.find("703_4230")     # sqrt(s) = 4.226 GeV, 1092 pb^-1
data_4260  = DatasetManager.real_data.find("703_4260")     # sqrt(s) = 4.258 GeV, 826 pb^-1
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")  # inclusive MC (kkmc + evtgen + lundcharm)
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

data_points = [data_4230, data_4260]
incMCs      = [incMC_4230, incMC_4260]

### Decay cards (EvtGen format) ###
# Non-resonant three-body phase-space process e+e- -> D*Dbar* pi0 with uniform
# momentum-phase-space distributions; the D* decays inclusively according to its PDG
# branching ratios (D* -> D pi and D* -> D gamma). The soft pi0/gamma from the D* is
# not used in the recoil mass.
# The top mother is psi(4260) (BESIII KKMC convention for e+e- annihilation).

# ---- Mode I: D0 -> K- pi+  and  Dbar0 -> K+ pi- ------------------------------
# RecID order (DecayCardResolver.rec_id_list):
#   0 psi(4260)   1 D*0      2 D0       3 K-     4 pi+
#   5 gamma(D*0)  6 anti-D*0 7 anti-D0  8 K+     9 pi-
#  10 gamma(anti-D*0)        11 pi0(bachelor)    12 gamma  13 gamma
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.0000  D*0  anti-D*0  pi0              PHSP;
    Enddecay

    Decay D*0
    1.0000  D0  gamma                       PHSP;
    Enddecay

    Decay anti-D*0
    1.0000  anti-D0  gamma                  PHSP;
    Enddecay

    Decay D0
    1.0000  K-  pi+                         PHSP;
    Enddecay

    Decay anti-D0
    1.0000  K+  pi-                         PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                    PHSP;
    Enddecay

    End
DECAYCARD

# ---- Mode II: D0 -> K- pi+ pi0  and  Dbar0 -> K+ pi- pi0 ---------------------
# RecID order:
#   0 psi(4260)  1 D*0  2 D0  3 K-  4 pi+  5 pi0(D0)  6 gamma  7 gamma  8 gamma(D*0)
#   9 anti-D*0  10 anti-D0  11 K+  12 pi-  13 pi0(anti-D0)  14 gamma  15 gamma
#  16 gamma(anti-D*0)  17 pi0(bachelor)  18 gamma  19 gamma
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.0000  D*0  anti-D*0  pi0              PHSP;
    Enddecay

    Decay D*0
    1.0000  D0  gamma                       PHSP;
    Enddecay

    Decay anti-D*0
    1.0000  anti-D0  gamma                  PHSP;
    Enddecay

    Decay D0
    1.0000  K-  pi+  pi0                    PHSP;
    Enddecay

    Decay anti-D0
    1.0000  K+  pi-  pi0                    PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                    PHSP;
    Enddecay

    End
DECAYCARD

# ---- Mode III: D0 -> K- pi+ pi+ pi-  and  Dbar0 -> K+ pi- pi+ pi- ------------
# RecID order:
#   0 psi(4260)  1 D*0  2 D0  3 K-  4 pi+  5 pi+  6 pi-  7 gamma(D*0)
#   8 anti-D*0  9 anti-D0  10 K+  11 pi-  12 pi+  13 pi-  14 gamma(anti-D*0)
#  15 pi0(bachelor)  16 gamma  17 gamma
decay_card_modeIII = <<~DECAYCARD
    Decay psi(4260)
    1.0000  D*0  anti-D*0  pi0              PHSP;
    Enddecay

    Decay D*0
    1.0000  D0  gamma                       PHSP;
    Enddecay

    Decay anti-D*0
    1.0000  anti-D0  gamma                  PHSP;
    Enddecay

    Decay D0
    1.0000  K-  pi+  pi+  pi-               PHSP;
    Enddecay

    Decay anti-D0
    1.0000  K+  pi-  pi+  pi-               PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                    PHSP;
    Enddecay

    End
DECAYCARD

# ---- Mode IV: D+ -> K- pi+ pi+  and  D- -> K+ pi- pi- -----------------------
# RecID order:
#   0 psi(4260)  1 D*+  2 D+  3 K-  4 pi+  5 pi+  6 gamma(D*+)
#   7 D*-  8 D-  9 K+  10 pi-  11 pi-  12 gamma(D*-)
#  13 pi0(bachelor)  14 gamma  15 gamma
decay_card_modeIV = <<~DECAYCARD
    Decay psi(4260)
    1.0000  D*+  D*-  pi0                   PHSP;
    Enddecay

    Decay D*+
    1.0000  D+  gamma                       PHSP;
    Enddecay

    Decay D*-
    1.0000  D-  gamma                       PHSP;
    Enddecay

    Decay D+
    1.0000  K-  pi+  pi+                    PHSP;
    Enddecay

    Decay D-
    1.0000  K+  pi-  pi-                    PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                    PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive signal MC — one sample per D decay mode and per energy point ###
exMCs_modeI = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Zc4025_D0KPi"
  config.events        = 200000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default   # ISR follows the measured D*Dbar*pi0 line shape
end

exMCs_modeII = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Zc4025_D0KPiPi0"
  config.events        = 200000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

exMCs_modeIII = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Zc4025_D0KPiPiPi"
  config.events        = 200000
  config.decay_card    = decay_card_modeIII
  config.cross_section = :default
end

exMCs_modeIV = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Zc4025_DpKPiPi"
  config.events        = 200000
  config.decay_card    = decay_card_modeIV
  config.cross_section = :default
end

### ---------------------------------------------------------------------------
### Event selection (BOSS) — Mode I: D0 -> K- pi+ (both sides)
### ---------------------------------------------------------------------------
alg_name_I = "Zc4025D0KPi"
alg_I = Algorithm.new(alg_name_I)
alg_I.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
     .set_alias({"std::vector<double>" => "Vdouble"})
     # Multi-energy analysis: ECMS is injected per dataset at run time.

sel_I = Selection.new
sel_I.select_track {
        cos_theta 0.93    # |cos(theta)| < 0.93, theta defined w.r.t. the e+ beam
        Vz        10.0    # closest approach to the IP within +-10 cm along the beam
        Vr        1.0     # closest approach within 1 cm in the plane perpendicular to the beam
        nChrp     "==2"   # K+ and pi+ (from the two D mesons)
        nChrn     "==2"   # K- and pi-
        nNet      "==0"
      }
      .select_photon {
        tdc_emc_start     0     # EMC cluster time within a 700 ns window around the
        tdc_emc_end       14    # event start time (electronic-noise suppression)
        angle_to_track    10.0  # photon must not be associated with a charged track
        energyThreshold_b 0.025 # E > 25 MeV in the EMC barrel region
        energyThreshold_e 0.050 # E > 50 MeV in the EMC end-cap region
        nGam              ">=2" # two photons from the bachelor pi0 -> gamma gamma
      }
      # A track is a K (pi) when the PID probabilities satisfy L(K) > L(pi)
      # (L(K) < L(pi)), using the dE/dx and TOF information.
      .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion]
        identify :pion, against: [:kaon]
        nkp  "==1"
        nkm  "==1"
        npip "==1"
        npim "==1"
      }
      # pi0 candidates from photon pairs constrained to the nominal pi0 mass
      # (the two-photon invariant mass is preselected in (0.120, 0.145) GeV/c^2).
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      }
      # Partial reconstruction: tag the D0 (from D*0 -> D0 gamma) and the bachelor pi0;
      # the D*0bar (together with the soft photon from the tagged D*0) is inferred from
      # the recoil four-momentum. The window is the D* recoil-mass signal region.
      .partial_rec([2, 11]) {
        best_combination_by_mass :D0, 1.86484   # M(K- pi+) closest to m(D0)
        require_recoil_mass 2.10, 2.22          # RM(D pi0) signal region (Fig. 1)
      }

alg_I.note(:pi0_preselection,
      "The invariant mass of any photon pair is required to be within (0.120, 0.145) GeV/c^2 " \
      "and is then constrained to the nominal pi0 mass; the kinematics of the two photons " \
      "are updated according to the constraint fit.")
     .note(:d_vertex_fit,
      "The charged tracks from a D decay candidate are required to originate from a common " \
      "vertex with chi2_VF < 100.")
     .note(:d_mass_constraint,
      "The reconstructed masses of the final-state particles are constrained to the " \
      "corresponding D nominal masses; chi2_KF(D) < 15 for D final states containing charged " \
      "tracks only and chi2_KF(D) < 20 for final states containing a pi0. A loose chi2_cut is " \
      "used in BOSS and the published values are applied in ROOT.")
     .note(:ddbar_pair,
      "Signal candidates must contain at least one pair of D Dbar candidates that do not share " \
      "final-state particles; if more than one pair exists, only the pair with the minimum " \
      "chi2_KF(D) + chi2_KF(Dbar) is retained.")
     .note(:bachelor_pi0,
      "The bachelor pi0 is reconstructed from the remaining photon showers that are not " \
      "assigned to the D Dbar pair. Each photon from the bachelor pi0 must not form a pi0 with " \
      "any other photon in the event; the two-photon mass constraint requires chi2_KF(pi0) < 20.")
     .note(:background_veto,
      "To reject the bachelor pi0 from D* -> D pi0 decays, the D pi0 invariant mass is required " \
      "to be greater than 2.02 GeV/c^2.")
     .note(:recoil_mass_oval,
      "Events are kept within the two-dimensional oval signal regions in the " \
      "RM(D pi0) versus RM(Dbar pi0) distributions; the dimensions differ between the two " \
      "energy points and are determined from MC simulation. Applied at ROOT level.")
     .note(:fit_observable,
      "The Z_c(4025)^0 signal is extracted from a simultaneous unbinned maximum-likelihood fit " \
      "to the RM(pi0) spectra at the two energies, using an S-wave Breit-Wigner with " \
      "mass-dependent width convolved with the detector resolution (4 MeV at 4.23 GeV and " \
      "4.5 MeV at 4.26 GeV) and weighted by the efficiency, with kernel-estimated " \
      "non-parametric background shapes. Performed at ROOT level.")
     .note(:signal_model,
      "The signal MC assumes Z_c(4025)^0 with J^P = 1^+ decaying to (D* Dbar*)^0 in an S wave, " \
      "with equal decay rates to the neutral and charged D* Dbar* channels (f_k = 0.5) from " \
      "isospin symmetry; the D* decays inclusively according to its PDG branching ratios.")

alg_I.with_decay_card(decay_card_modeI).apply(sel_I)
alg_I.execute_on(data_points + incMCs + exMCs_modeI)

### ---------------------------------------------------------------------------
### Event selection (BOSS) — Mode II: D0 -> K- pi+ pi0 (both sides)
### ---------------------------------------------------------------------------
alg_name_II = "Zc4025D0KPiPi0"
alg_II = Algorithm.new(alg_name_II)
alg_II.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
      .set_alias({"std::vector<double>" => "Vdouble"})

sel_II = Selection.new
sel_II.select_track {
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==2"   # K+ and pi+ (from the two D mesons)
       nChrn     "==2"   # K- and pi-
       nNet      "==0"
     }
     .select_photon {
       tdc_emc_start     0
       tdc_emc_end       14
       angle_to_track    10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam              ">=4"  # two photons per pi0: one per D and one for the bachelor pi0
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :kaon, against: [:pion]
       identify :pion, against: [:kaon]
       nkp  "==1"
       nkm  "==1"
       npip "==1"
       npim "==1"
     }
     # pi0 candidates (from the D0 -> K- pi+ pi0 decays and from the bachelor pi0)
     .kalman_kinematic_fit([:gamma, :gamma]) {
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
       chi2_cut 25
       npi0 ">=3"   # one pi0 in each D decay plus the bachelor pi0
     }
     # Partial reconstruction: tag the D0 (recID 2) and the bachelor pi0 (recID 17);
     # the D*0bar is inferred from the recoil four-momentum.
     .partial_rec([2, 17]) {
       best_combination_by_mass :D0, 1.86484   # M(K- pi+ pi0) closest to m(D0)
       require_recoil_mass 2.10, 2.22          # RM(D pi0) signal region
     }

alg_II.note(:pi0_preselection,
      "The invariant mass of any photon pair is required to be within (0.120, 0.145) GeV/c^2 " \
      "and is then constrained to the nominal pi0 mass.")
      .note(:d_vertex_fit,
      "The charged tracks from a D decay candidate are required to originate from a common " \
      "vertex with chi2_VF < 100.")
      .note(:d_mass_constraint,
      "The final-state particles are constrained to the corresponding D nominal masses; " \
      "chi2_KF(D) < 20 for D final states containing a pi0 (loose chi2_cut in BOSS, published " \
      "value applied in ROOT).")
      .note(:ddbar_pair,
      "At least one pair of D Dbar candidates sharing no final-state particles is required; " \
      "the pair with the minimum chi2_KF(D) + chi2_KF(Dbar) is retained.")
      .note(:bachelor_pi0,
      "The bachelor pi0 is built from the photon showers left over after the D Dbar assignment; " \
      "its photons must not form a pi0 with any other photon, and chi2_KF(pi0) < 20.")
      .note(:background_veto,
      "D pi0 invariant mass > 2.02 GeV/c^2 to reject the bachelor pi0 from D* -> D pi0 decays.")
      .note(:recoil_mass_oval,
      "Events are kept within the two-dimensional oval signal regions in RM(D pi0) versus " \
      "RM(Dbar pi0); applied at ROOT level.")
      .note(:signal_model,
      "The Z_c(4025)^0 is generated with J^P = 1^+ and equal neutral/charged D* Dbar* decay " \
      "rates (f_k = 0.5); the D* decays inclusively per PDG.")

alg_II.with_decay_card(decay_card_modeII).apply(sel_II)
alg_II.execute_on(data_points + incMCs + exMCs_modeII)

### ---------------------------------------------------------------------------
### Event selection (BOSS) — Mode III: D0 -> K- pi+ pi+ pi- (both sides)
### ---------------------------------------------------------------------------
alg_name_III = "Zc4025D0KPiPiPi"
alg_III = Algorithm.new(alg_name_III)
alg_III.set_header(["#{alg_name_III}Alg/#{alg_name_III}.h"])
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_III = Selection.new
sel_III.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==4"   # three pi+ and one K+ from the two D mesons
        nChrn     "==4"   # three pi- and one K-
        nNet      "==0"
      }
      .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=2"  # photons from the bachelor pi0
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion]
        identify :pion, against: [:kaon]
        nkp  "==1"
        nkm  "==1"
        npip "==3"
        npim "==3"
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      }
      # Partial reconstruction: tag the D0 (recID 2) and the bachelor pi0 (recID 15)
      .partial_rec([2, 15]) {
        best_combination_by_mass :D0, 1.86484   # M(K- pi+ pi+ pi-) closest to m(D0)
        require_recoil_mass 2.10, 2.22          # RM(D pi0) signal region
      }

alg_III.note(:pi0_preselection,
      "The two-photon invariant mass is required to be within (0.120, 0.145) GeV/c^2 and " \
      "constrained to the nominal pi0 mass.")
       .note(:d_vertex_fit,
      "The charged tracks from a D decay candidate must originate from a common vertex with " \
      "chi2_VF < 100.")
       .note(:d_mass_constraint,
      "The final-state particles are constrained to the D nominal masses; chi2_KF(D) < 15 for " \
      "the all-charged D decay modes (loose chi2_cut in BOSS, published value applied in ROOT).")
       .note(:ddbar_pair,
      "At least one D Dbar pair sharing no final-state particles is required; the minimum " \
      "chi2_KF(D) + chi2_KF(Dbar) combination is kept.")
       .note(:bachelor_pi0,
      "The bachelor pi0 is built from the leftover photon showers; its photons must not form a " \
      "pi0 with any other photon and chi2_KF(pi0) < 20.")
       .note(:background_veto,
      "D pi0 invariant mass > 2.02 GeV/c^2.")
       .note(:mixed_d_modes,
      "The D0 and the Dbar0 are allowed to be reconstructed in different decay modes; mixed-mode " \
      "D Dbar combinations are handled together with the same-mode combinations by considering " \
      "all possible combinations of selected charged tracks and pi0.")
       .note(:recoil_mass_oval,
      "Two-dimensional oval signal regions in RM(D pi0) versus RM(Dbar pi0); applied at ROOT level.")

alg_III.with_decay_card(decay_card_modeIII).apply(sel_III)
alg_III.execute_on(data_points + incMCs + exMCs_modeIII)

### ---------------------------------------------------------------------------
### Event selection (BOSS) — Mode IV: D+ -> K- pi+ pi+ (charged channel)
### ---------------------------------------------------------------------------
alg_name_IV = "Zc4025DpKPiPi"
alg_IV = Algorithm.new(alg_name_IV)
alg_IV.set_header(["#{alg_name_IV}Alg/#{alg_name_IV}.h"])
      .set_alias({"std::vector<double>" => "Vdouble"})

sel_IV = Selection.new
sel_IV.select_track {
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==3"   # two pi+ and one K+ from the two D mesons
       nChrn     "==3"   # two pi- and one K-
       nNet      "==0"
     }
     .select_photon {
       tdc_emc_start     0
       tdc_emc_end       14
       angle_to_track    10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam              ">=2"  # photons from the bachelor pi0
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :kaon, against: [:pion]
       identify :pion, against: [:kaon]
       nkp  "==1"
       nkm  "==1"
       npip "==2"
       npim "==2"
     }
     .kalman_kinematic_fit([:gamma, :gamma]) {
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
       chi2_cut 25
       npi0 ">=1"
     }
     # Partial reconstruction: tag the D+ (recID 2) and the bachelor pi0 (recID 13)
     .partial_rec([2, 13]) {
       best_combination_by_mass :D_plus, 1.86962   # M(K- pi+ pi+) closest to m(D+)
       require_recoil_mass 2.10, 2.22              # RM(D pi0) signal region
     }

alg_IV.note(:pi0_preselection,
      "The two-photon invariant mass is required to be within (0.120, 0.145) GeV/c^2 and " \
      "constrained to the nominal pi0 mass.")
      .note(:d_vertex_fit,
      "The charged tracks from a D decay candidate must originate from a common vertex with " \
      "chi2_VF < 100.")
      .note(:d_mass_constraint,
      "The final-state particles are constrained to the D nominal masses; chi2_KF(D) < 15 for " \
      "the all-charged D+ -> K- pi+ pi+ mode (loose chi2_cut in BOSS, published value applied " \
      "in ROOT).")
      .note(:ddbar_pair,
      "At least one D+ D- pair sharing no final-state particles; the minimum " \
      "chi2_KF(D+) + chi2_KF(D-) combination is retained.")
      .note(:bachelor_pi0,
      "The bachelor pi0 is built from the leftover photon showers, with the requirement that " \
      "its photons do not form a pi0 with any other photon and chi2_KF(pi0) < 20.")
      .note(:background_veto,
      "D pi0 invariant mass > 2.02 GeV/c^2 to reject D* -> D pi0 feed-down.")
      .note(:recoil_mass_oval,
      "Two-dimensional oval signal regions in RM(D pi0) versus RM(Dbar pi0); applied at ROOT level.")

alg_IV.with_decay_card(decay_card_modeIV).apply(sel_IV)
alg_IV.execute_on(data_points + incMCs + exMCs_modeIV)
