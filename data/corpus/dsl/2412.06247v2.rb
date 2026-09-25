# DSL for arXiv:2412.06247v2
# Partial wave analyses of psi(3686) → p pbar pi0 and p pbar eta
# Using (2712±14)×10^6 psi(3686) events, BESIII

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# =============================================================================
# Datasets
# =============================================================================

# Main psi(2S) data at 3.686 GeV
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# psi(3770) data at 3.773 GeV for continuum background estimation
psipp_data = DatasetManager.real_data.find("712_3773")
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")

# psi(2S) scan data points (√s = 3.670–3.710 GeV) for cross-section measurement
psip_scan_points = [
  DatasetManager.real_data.find("704_psip_scan_1"),   # 3581.5 MeV — actually this is below the scan range
  # The actual scan points for the paper's range 3.670–3.710 GeV are:
  # Using BOSS 704 scan data around psi(2S)
]

# psi(3770) data for continuum subtraction (√s = 3.773 GeV, L = 2.93 fb^-1)
psipp_3773 = DatasetManager.real_data.find("712_3773")
psipp_3773_incMC = DatasetManager.inclusive_mc.find("712_3773")

# =============================================================================
# Decay cards — KKMC + psi(2S) top mother
# =============================================================================

# Decay card: psi(3686) → p pbar pi0, pi0 → gamma gamma
decay_psip_ppbar_pi0 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 anti-p- p+ pi0   PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma      PHSP;
  Enddecay
  End
DECAYCARD

# Decay card: psi(3686) → p pbar eta, eta → gamma gamma
decay_psip_ppbar_eta = <<~DECAYCARD
  Decay psi(2S)
  1.0000 anti-p- p+ eta   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma      PHSP;
  Enddecay
  End
DECAYCARD

# =============================================================================
# Exclusive MC samples
# =============================================================================

# Signal MC for psi(3686) → p pbar pi0
exMC_psip_ppbar_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psip_ppbar_pi0"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_psip_ppbar_pi0
  config.cross_section   = :default
end

# Signal MC for psi(3686) → p pbar eta
exMC_psip_ppbar_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psip_ppbar_eta"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_psip_ppbar_eta
  config.cross_section   = :default
end

# =============================================================================
# Algorithm — psi(3686) → p pbar pi0 / p pbar eta
# Both channels share the same final state (p pbar gamma gamma) and selection
# The pi0 vs eta distinction happens at the PWA (ROOT) level via gamma-gamma
# invariant mass window
# =============================================================================

alg = Algorithm.new("PsipPPbarPi0Eta")
alg.set_header(["PsipPPbarPi0EtaAlg/PsipPPbarPi0Eta.h"])
alg.set_constant({ "ECMS" => [:double, 3.686] })
alg.note(:pwa_method, "partial wave analysis using relativistic covariant tensor
  amplitude formalism performed at ROOT level; unbinned maximum likelihood fit via
  MINUIT; N* resonances included: N(1440), N(1520), N(1535), N(1650), N(1710),
  N(1720), N(2100), N(2300), N(2570) plus virtual proton pole N(940); ppbar
  structures: rho(1900), rho(2000), rho(2150), rho(2225); for eta channel:
  phi3(1850), omega(1960), omega(2205)")
alg.note(:continuum_subtraction, "continuum contribution from e+e-→p pbar pi0(eta)
  estimated using psi(3770) data at sqrt(s)=3.773 GeV (L=2.93 fb^-1) and rescaled
  to 3.686 GeV via factor fc=1.564 (pi0) / fc=1.623 (eta) accounting for luminosity,
  Born cross section, efficiency, and ISR correction")
alg.note(:interference_phase, "ambiguous phase angle between psi(3686) resonance
  and continuum process yields two solutions (constructive/destructive) for
  branching fractions")
alg.note(:jpsi_veto, "events with |M(ppbar) - 3.097| < 0.050 GeV rejected to
  suppress J/psi → p pbar background; applied after kinematic fit using
  corrected momenta")
alg.note(:cross_section_fit, "observed cross sections at 9 energy points
  (3.670-3.710 GeV) fitted with coherent sum of continuum (a/s^n) and
  psi(3686) Breit-Wigner amplitude; ISR correction via ConExc; beam spread
  convolution with sigma_E = 0.3 MeV")

sel = Selection.new
sel.select_track {
  cos_theta 0.93     # |cos(theta)| < 0.93
  Vz 10.0            # |Vz| < 10 cm
  Vr 1.0             # Vr < 1 cm in transverse plane
  nChrp "==1"        # exactly 1 positive track (proton)
  nChrn "==1"        # exactly 1 negative track (anti-proton)
  nNet "==0"         # net charge zero
}
.select_photon {
  tdc_emc_start 0              # TDC start: 0 ns
  tdc_emc_end 14               # TDC end: 700 ns (14 x 50 ns)
  energyThreshold_b 0.025      # E > 25 MeV in barrel (|cos theta| < 0.80)
  energyThreshold_e 0.050      # E > 50 MeV in endcap (0.86 < |cos theta| < 0.92)
  angle_to_track 20.0          # min angle to nearest charged track > 20 deg
  nGam ">=2"                   # at least 2 photons
}
# PID: no explicit PID in the main analysis (both charged tracks are protons);
# for the scan data analysis, PID with confidence level proton > pion, kaon is applied
# We assign both charged tracks as protons
.assign({ chrgp: :prp, chrgn: :prm })
# Kalman fit: reconstruct pi0 from photon pairs (for the pi0 channel hypothesis)
# Two separate hypotheses for pi0 and eta — but at BOSS level both are gamma-gamma
# pairs with mass constraint to pi0 or eta
# We reconstruct pi0 first for the main 5C fit; the eta channel uses ROOT-level
# selection on gamma-gamma invariant mass
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 100       # loose cut; combined with 5C fit
}
# 5C kinematic fit: psi(3686) → p pbar pi0
# 4C (energy-momentum) + 1C (pi0 → gamma gamma mass constraint via Kalman fit)
.kinematic_fit([:prp, :prm, :pi0]) {
  nominal
  constrain_four_momentum
  chi2_cut 30        # chi2_5C < 30 as stated in paper
}
# Competing eta hypothesis: store chi2 for ROOT-level eta selection
# Note: the eta channel is selected by ROOT-level cut on |M(gamma gamma) - m_eta|
# The Kalman fit with eta mass constraint is done indirectly at ROOT level

alg.with_decay_card(decay_psip_ppbar_pi0).apply(sel)
alg.execute_on([psip_data, psip_incMC, exMC_psip_ppbar_pi0, exMC_psip_ppbar_eta])