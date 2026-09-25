# Observation of D0 -> K1(1270)- e+ nu_e semileptonic decay
# sqrt(s) = 3.773 GeV, 2.93 fb^-1; double-tag (DT) analysis with pre-stored DTag candidates.
# Tag side: anti-D0 (D0bar) -> K+pi-, K+pi-pi0, K+pi-pi-pi+.
# Signal side: D0 -> K1(1270)- e+ nu_e, K1(1270)- -> K-pi+pi-.

### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ---------------------------------------------------------------------
# Decay card: psi(3770) -> D0 D0bar.
# Tag side: anti-D0 -> hadronic ST modes (equal fractions, per-mode
# efficiency weighting applied in ROOT).
# Signal side: D0 -> K1(1270)- e+ nu_e, K1- -> K- pi+ pi-.
# ---------------------------------------------------------------------
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0  PHSP;
    Enddecay

    Decay D0
    1.0000 K_1- e+ nu_e  PHOTOS ISGW2;
    Enddecay

    Decay anti-D0
    0.3333 K+ pi-        PHSP;
    0.3333 K+ pi- pi0    PHSP;
    0.3334 K+ pi- pi- pi+ PHSP;
    Enddecay

    Decay K_1-
    1.0000 K- pi+ pi-    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma   PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "D0_K1_enu"
  config.related_dataset = data_3773
  config.events          = 400_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# ---------------------------------------------------------------------
# TagAnalysis: single-tag anti-D0 + signal-side D0 -> K1- e+ nu_e
# ---------------------------------------------------------------------
alg_name = "D0SemileptonicK1"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(decay_card_signal)

# Tag side: anti-D0 reconstructed via three hadronic D0bar flavour modes.
# charm -1 selects the D0bar charge conjugate of the declared D0 modes:
#   D0toKPi      -> K+ pi-
#   D0toKPiPi0   -> K+ pi- pi0
#   D0toKPiPiPi  -> K+ pi- pi- pi+
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: D0 -> K1(1270)- e+ nu_e with K1- -> K- pi+ pi-.
# The four signal-side tracks are: K- (km), pi+ (pip), pi- (pim), e+ (ep).
# Total charge = -1 + 1 - 1 + 1 = 0.
alg.signal_side do |s|
  s.charged(km: 1, pip: 1, pim: 1, ep: 1)
  s.require_charge 0
  s.missing :nu_e       # massless neutrino (semileptonic)
end

# 4C kinematic fit: tag D + K- + pi+ + pi- + e+ + nu_e = measured CMS 4-vector.
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# ---------------------------------------------------------------------
# Notes: selection criteria applied downstream in ROOT (no DSL constructs)
# ---------------------------------------------------------------------
alg.note(:tag_deltaE_windows,
         "ST deltaE windows (mode-dependent, applied in ROOT): " \
         "K+pi-       [-0.029, 0.027] GeV; " \
         "K+pi-pi0    [-0.069, 0.038] GeV; " \
         "K+pi-pi-pi+ [-0.031, 0.028] GeV. " \
         "These are ~3-sigma windows; the asymmetry in the K+pi-pi0 window " \
         "reflects the tail from the ISR and FSR photon emission.")
   .note(:tag_mBC_window,
         "ST beam-constrained mass M_BC in (1.858, 1.874) GeV/c^2, " \
         "applied in ROOT. The M_BC resolution is mode-dependent but the " \
         "same window covers all three tag modes with > 99% efficiency.")
   .note(:track_selection,
         "All charged tracks satisfy |cos(theta)| < 0.93, distance of closest " \
         "approach to the IP < 1 cm in the transverse plane and < 10 cm along " \
         "the beam. Tracks from K_S0 decays are not used.")
   .note(:pid_requirements,
         "Kaon and pion PID combine dE/dx and TOF: K candidates require " \
         "L(K) > L(pi); pi candidates require L(pi) > L(K). Electron PID " \
         "combines dE/dx, TOF and EMC: E/p > 0.8 and L'(e)/(L'(e)+L'(pi)+L'(K)) > 0.8.")
   .note(:photon_selection,
         "EMC showers with E > 25 MeV (barrel, |cos(theta)| < 0.80) or " \
         "E > 50 MeV (end cap, 0.86 < |cos(theta)| < 0.92), separated from " \
         "extrapolated charged tracks by > 10 degrees. Shower time within 700 ns " \
         "of the event start time.")
   .note(:pi0_reconstruction,
         "pi0 candidates from photon pairs with 0.115 < M(gamma gamma) < 0.150 GeV/c^2, " \
         "followed by a 1C mass-constrained fit to the nominal pi0 mass.")
   .note(:candidate_ranking,
         "When multiple ST candidates exist in an event, only the one with the " \
         "smallest |deltaE| is kept per tag mode.")
   .note(:k1_mass_window,
         "The K1(1270)- resonance is identified via the K-pi+pi- invariant mass: " \
         "the signal region and sidebands are defined in the ROOT analysis. The " \
         "K1 mass and width are extracted from a simultaneous fit to the M(Kpipi) " \
         "and U_miss distributions.")
   .note(:extra_energy_veto,
         "The total energy of unused EMC showers E_extra is required to be below a " \
         "threshold (typically 0.3 GeV) to suppress backgrounds with additional pi0 " \
         "or radiative photons. Applied in the ROOT analysis.")
   .note(:umiss_extraction,
         "The DT signal yield is extracted from an unbinned maximum-likelihood fit " \
         "to U_miss = E_miss - c|p_miss|, where E_miss and p_miss are computed from " \
         "the tagged D, the four signal-side tracks and the beam constraint. The " \
         "signal peaks at zero. Backgrounds are modelled with MC-derived shapes.")
   .note(:fit_extraction,
         "The M(Kpipi) spectrum is simultaneously fitted with the U_miss distribution " \
         "to extract the K1(1270)- signal yield, mass and width. The mass resolution " \
         "function is taken from MC and convoluted with a relativistic Breit-Wigner.")

alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_signal])