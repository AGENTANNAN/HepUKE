# Observation of three charmoniumlike states with JPC = 1-- in e+e- -> D*0 D*- pi+
#   [arXiv:2301.07321]
#
# Born cross section measurement at 86 CMS energies from 4.189 to 4.951 GeV
# Partial reconstruction technique with two tagging methods: D0-tag and D--tag
# 17.9 fb-1 total integrated luminosity
# Multi-energy scan: NO single ECMS (skip in set_constant)
# Ordinary analysis: Algorithm + Selection with partial_rec

### Dataset preparation ###
# Data spans BOSS 703, 706, 707 across many energy points.
# The analysis uses both "XYZ data" (37 points, large luminosity) and
# "scan data" (49 points, small luminosity).
#
# XYZ data energy points (from Table I of Supplemental Material):
# 4.189, 4.199, 4.209, 4.219, 4.226, 4.236, 4.242*, 4.244, 4.258,
# 4.267, 4.278*, 4.287, 4.308*, 4.311, 4.337, 4.358, 4.377, 4.387*,
# 4.395, 4.416, 4.436, 4.467*, 4.527*, 4.575*, 4.600, 4.613*,
# 4.628, 4.641, 4.661, 4.682, 4.699, 4.740*, 4.750*, 4.781, 4.843,
# 4.918*, 4.951*
# (* = lower luminosity point)
#
# BOSS 703 covers: 4.189-4.226 (and ), 4.236-4.278, 4.287-4.387, 4.395-4.467
# BOSS 706 covers: 4.600-4.699
# BOSS 707 covers: 4.740-4.951
#
# For brevity, we enumerate representative energy points and note that the
# full 86-point list should be loaded from a config file.
#
# Key representative datasets:
data_4189 = DatasetManager.real_data.find("703_4180")   # 4.178 + 4.189 group
data_4199 = DatasetManager.real_data.find("703_4190")   # 4.189 + scan
data_4209 = DatasetManager.real_data.find("703_4200")   # 4.199 + scan
data_4219 = DatasetManager.real_data.find("703_4210")   # 4.209 + scan
data_4226 = DatasetManager.real_data.find("703_4220")   # 4.219 + scan
data_4230 = DatasetManager.real_data.find("703_4230")   # 4.226 + scan
data_4600 = DatasetManager.real_data.find("703_4600")   # 4.600
data_4620 = DatasetManager.real_data.find("706_4620")   # 4.628
data_4640 = DatasetManager.real_data.find("706_4640")   # 4.641
data_4660 = DatasetManager.real_data.find("706_4660")   # 4.661
data_4680 = DatasetManager.real_data.find("706_4680")   # 4.682
data_4700 = DatasetManager.real_data.find("706_4700")   # 4.699
data_4740 = DatasetManager.real_data.find("707_4740")   # 4.740
data_4750 = DatasetManager.real_data.find("707_4750")   # 4.750
data_4780 = DatasetManager.real_data.find("707_4780")   # 4.781
data_4840 = DatasetManager.real_data.find("707_4840")   # 4.843
# ...remaining scan points enumerated similarly from the config

incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
# ... (all energy points have corresponding inclusive MC)

all_data = [data_4189, data_4199, data_4209, data_4219, data_4226, data_4230,
            data_4600, data_4620, data_4640, data_4660, data_4680, data_4700,
            data_4740, data_4750, data_4780, data_4840]
all_incMC = [incMC_4180, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4230,
             incMC_4600]
# Note: the full 86-point dataset list should be loaded from a YAML config file;
# this is a representative subset.

# ---- Decay card -----------------------------------------------------------
# Signal: e+ e- -> D*0 D*- pi+, with partial reconstruction
# Sub-decays for the D0/D- tag modes
decay_card = <<~DECAYCARD
  Decay vpho
  1.0000 D*0 D*- pi+   PHSP;
  Enddecay

  Decay D*0
  1.0000 D0 pi0   VSS;
  Enddecay

  Decay D*-
  1.0000 D- pi0   VSS;
  Enddecay

  Decay D0
  1.0000 K- pi+   PHSP;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi-   PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay

  End
DECAYCARD

# Signal MC with PWA results at each energy point
exMC = DatasetManager.create_exclusive_mc_for(all_data) do |c|
  c.sample_name     = "ee_to_Dstar0DstarPi"
  c.events          = 200_000
  c.decay_card      = decay_card
  c.cross_section   = :default
end

### Event selection (BOSS) — Ordinary analysis with partial reconstruction ###
alg = Algorithm.new("Dstar0DstarPi_CrossSection")

# Multi-energy: skip ECMS (energy injected at runtime via MeasuredEcmsSvc)
alg.set_header(["Dstar0DstarPiAlg/Dstar0DstarPi.h"])

sel = Selection.new

# Charged tracks: bachelor pi+ + D0/D- daughters
# D0-tag: pi+ + K- + pi+ (Kpi), or pi+ + K- + pi+ + pi0 (Kpipi0), or pi+ + K- + 2pi+ + pi- (K3pi)
# D--tag: pi+ + K+ + 2pi- (Kpipi)
# Minimum: pi+ + K- + K+ + pi- + pi+ = 5 charged tracks (one tag at a time)
sel.select_track {
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     ">=2"
       nChrn     ">=2"
       nTot      ">=5"
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :kp, :km, against: [:pion]
       identify :pip, :pim, against: [:kaon]
     }
     # At least 2 photons for the pi0 from D* decay
     .select_photon {
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       angle_to_track    10.0
       tdc_emc_start     0
       tdc_emc_end       14
       nGam              ">=2"
     }

# Reconstruct pi0 from gamma gamma (1C Kalman fit)
sel.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
end

# Partial reconstruction: reconstruct D*0 (or D*-) + pi+, infer the missing
# D*- (or D*0) as recoil. The decay card has:
#   recID 0 = vpho (skip)
#   recID 1 = D*0, recID 2 = D*-, recID 3 = pi+_bachelor
#   recID 4 = D0, recID 5 = pi0_1 (from D*0)
#   recID 6 = D-, recID 7 = pi0_2 (from D*-)
#   recID 8 = K-, recID 9 = pi+_D0, recID 10 = K+, recID 11 = pi-_D-, recID 12 = pi-_D-2
#   recID 13,14 = gamma, gamma (from pi0)
#
# D0-tag: reconstruct D*0 (= D0 + pi0) + bachelor pi+, miss D*- (= D- + pi0)
# D--tag: reconstruct D*- (= D- + pi0) + bachelor pi+, miss D*0 (= D0 + pi0)
# Both tagging methods are handled simultaneously.
#
# For the D0-tag case:
#   tagged IDs: D*0(1), pi+(3)
#   missed IDs: D*-(2) and its daughters
sel.partial_rec([1, 3]) do
  # D0-tag: D*0 mass constrained via kalman_kinematic_fit + D0 mass window
  # D--tag: D*- mass constrained via kalman_kinematic_fit + D- mass window
  # The best_combination_by_mass selects the D* candidate closest to nominal mass
  require_recoil_mass 1.80, 2.20   # wide window on missing D* mass
end

# 3C kinematic fit: constrain pi0, D, D* masses to their known values
# This is performed after partial reconstruction to improve resolution.
# Participants: pi0 + D0/D- + D*0/D*- with the reconstructed/fitted momenta
# But since partial_rec replaces kinematic_fit, we note this in .note()
sel.note(:kinematic_fit_3c,
         "A 3C kinematic fit is performed constraining the reconstructed pi0, " \
         "D (D0 or D-), and D* (D*0 or D*-) mesons to their known masses. " \
         "Events with chi2_3C < 50 are retained. The fitted four-momenta are " \
         "used for all subsequent analysis. The combination with minimum chi2_3C " \
         "is kept. Note: this 3C fit is performed in the C++ algorithm code, not " \
         "at the DSL level, since partial_rec replaces the standard kinematic_fit block.")

alg.with_decay_card(decay_card).apply(sel)

# ---- Notes for inexpressible BOSS procedures ----
alg.note(:energy_points,
        "The analysis covers 86 CMS energy points from 4.189 to 4.951 GeV " \
        "(37 XYZ points with large luminosity and 49 scan points). " \
        "The full dataset list should be loaded from a config file. " \
        "The representative subset shown in this spec covers the main XYZ points. " \
        "Energy-dependent beam conditions (ECMS, beam energy spread) are read " \
        "from the conditions database at runtime via MeasuredEcmsSvc.")
  .note(:two_tagging_methods,
        "Two tagging methods are used: D0-tag (D0 -> K-pi+, K-pi+pi0, K-pi+pi+pi-) " \
        "and D--tag (D- -> K+pi-pi-). " \
        "In the D0-tag method, the D0, a bachelor pi+, and at least one soft pi0 " \
        "from D*0 -> D0 pi0 are reconstructed. " \
        "In the D--tag method, the D-, a bachelor pi+, and at least one soft pi0 " \
        "from D*- -> D- pi0 are reconstructed. " \
        "If an event survives both tag methods, only the D0-tag combination is kept.")
  .note(:d_mass_windows,
        "D0 mass windows: K-pi+ [1.835, 1.887], K-pi+pi0 [1.827, 1.882], " \
        "K-pi+pi+pi- [1.855, 1.874] GeV/c^2. " \
        "D- mass window: K+pi-pi- [1.856, 1.883] GeV/c^2. " \
        "These are applied at the ROOT level after reconstruction.")
  .note(:dstar_mass_and_pi0_veto,
        "D0-tag: M(D0 pi0) in [2.004, 2.009] GeV/c^2 with P*(pi0) not in [0.025, 0.050] GeV/c. " \
        "D--tag: M(D- pi0) in [2.008, 2.013] GeV/c^2 with P*(pi0) not in [0.030, 0.055] GeV/c. " \
        "The P*(pi0) veto removes pi0 from the missing D* candidate. " \
        "The pi+ D0 invariant mass must be > 2.02 GeV/c^2 in D0-tag to reject " \
        "bachelor pi+ from D*+ -> pi+ D0. " \
        "These cuts are applied at the ROOT level.")
  .note(:signal_extraction,
        "Signal yields extracted from simultaneous unbinned maximum likelihood fits " \
        "to RM(pi+ pi0 D0) and RM(pi+ pi0 D-) distributions. " \
        "Signal shape: MC convolved with Gaussian. Background: 2nd-order Chebyshev " \
        "plus peaking background from mis-combination (fixed from inclusive MC). " \
        "The Born cross section is calculated iteratively with ISR and vacuum " \
        "polarization corrections (ROOT level).")
  .note(:chic2_3c_cut,
        "chi2_3C < 50 applied. If multiple pi0 D candidates exist, " \
        "the one with minimum chi2_3C is retained.")
  .note(:systematics,
        "Systematic uncertainties (6.7-9.6%): tracking (1.0%/track), " \
        "PID (1.0%/track), pi0 reconstruction (2.0%), signal region requirements, " \
        "signal decay model (PWA results), ISR correction factor, luminosity (1.0%), " \
        "quoted BFs. All evaluated per tag method and combined by signal yield.")
  .note(:lineshape_fit,
        "The dressed cross section line shape is fitted with a coherent sum of " \
        "a continuum amplitude and three relativistic Breit-Wigner functions. " \
        "Eight degenerate solutions found. Three resonances observed: " \
        "psi(4210), psi(4470), psi(4660) with > 10 sigma significance.")

alg.execute_on(all_data + all_incMC + exMC)