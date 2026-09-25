# DSL for arXiv:2412.10803v1
# Study of the semileptonic decay D0 → Kbar0 pi- e+ nu_e
# Using 7.9 fb^-1 at sqrt(s)=3.773 GeV, BESIII
# TagAnalysis: ST + missing (semileptonic)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# =============================================================================
# Datasets — psi(3770) at 3.773 GeV
# =============================================================================

psipp_data = DatasetManager.real_data.find("712_3773")
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")

# =============================================================================
# Decay card — KKMC + psi(3770) top mother
# D0 → Kbar0 pi- e+ nu_e
# Kbar0 → K_S0; K_S0 → pi+ pi-
# =============================================================================

decay_D0_Kbar0_pi_e_nu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0   PHSP;
  Enddecay
  Decay D0
  1.0000 anti-K0 pi- e+ nu_e   ISGW2;
  Enddecay
  Decay anti-D0
  1.0000 K+ pi-   PHSP;
  Enddecay
  Decay anti-K0
  1.0000 K_S0   PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

# =============================================================================
# Exclusive MC sample — semileptonic signal
# =============================================================================

exMC_D0_Kbar0_pi_e_nu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_D0_Kbar0_pi_e_nu"
  config.related_dataset = psipp_data
  config.events          = 500_000
  config.decay_card      = decay_D0_Kbar0_pi_e_nu
  config.cross_section   = :default
end

# =============================================================================
# TagAnalysis — ST anti-D0 tags + semileptonic signal D0 → Kbar0 pi- e+ nu_e
# Tag modes (anti-D0): D0bar → K+ pi-, K+ pi- pi- pi+, K+ pi- pi0
# Signal: K_S0 (pi+ pi-), pi-, e+, nu_e (missing)
# =============================================================================

alg = TagAnalysis.new("D0TagKbar0PiENu")
alg.set_header(["D0TagKbar0PiENuAlg/D0TagKbar0PiENu.h"])
alg.set_constant({ "ECMS" => [:double, 3.773] })
alg.note(:st_yield, "ST yields extracted from M_BC fits with MC shape convolved
  with double-Gaussian; ARGUS background; summed over 3 tag modes:
  N_ST = 6306.7±2.9 × 10^3; WS peaking backgrounds subtracted from simulation")
alg.note(:dt_yield, "DT semileptonic yield N_DT = 8752±132 from U_miss fit;
  U_miss = E_miss - |p_miss|c peaks at zero for signal; signal shape from
  simulated events convolved with Gaussian, background from inclusive MC")
alg.note(:K_S0_reco, "K_S0 reconstructed via pi+ pi- pair with invariant mass
  (0.485, 0.510) GeV/c^2; secondary vertex fit constrains to common vertex with
  decay length significance > 2 sigma — applied at ROOT level")
alg.note(:electron_pid, "Electron PID: L'_e > 0.001 and L'_e/(L'_e+L'_pi+L'_K) > 0.8;
  E/p > 0.7; bremsstrahlung recovery adding EMC showers within 5 deg of electron
  direction; E_gamma_max < 0.25 GeV to suppress pi0 backgrounds")
alg.note(:hadron_veto, "M(Kbar0 pi- e+) < 1.80 GeV/c^2 rejects D0 → Kbar0 pi+ pi-
  background; D0bar → K+ pi- pi- pi+ reconstructed as wrong-sign background estimated
  from simulation")
alg.note(:form_factors, "Five-dimensional unbinned ML fit to m(Kbar0 pi-), q^2,
  cos(theta_e), cos(theta_Kbar0), chi; S-wave + P-wave (K*(892)-) model;
  form factor ratios r_V = V(0)/A1(0) = 1.48±0.05±0.02,
  r_2 = A2(0)/A1(0) = 0.70±0.04±0.02 from ROOT-level fit")
alg.note(:bf, "B(D0 → Kbar0 pi- e+ nu_e) = (1.444±0.022±0.024)%;
  B fractions: K*(892)- = (94.15±0.32±0.16)%, S-wave = (5.87±0.32±0.16)%")

# Tag side: anti-D0 reconstructed in three hadronic modes
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPiPi, :D0toKPiPi0
  t.charm -1            # anti-D0 tag
end

# Signal side: K_S0 (pi+ pi-), pi- (from D0), e+ (from D0), nu_e (missing)
alg.signal_side do |s|
  s.photons 0           # no photon required on signal side
  s.charged(pip: 1, pim: 2, ep: 1)   # K_S0 daughters + signal pi-, e+
  s.require_charge 0    # D0 → Kbar0 pi- e+ nu_e: 0 + (-1) + 1 + 0 = 0
  s.missing :nu_e       # massless neutrino (semileptonic)
end

# Kinematic fit: 4C (tag + signal + nu = ecms_lab)
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200        # loose cut; tight cut applied in ROOT
end

alg.with_decay_card(decay_D0_Kbar0_pi_e_nu).apply
alg.execute_on([psipp_data, psipp_incMC, exMC_D0_Kbar0_pi_e_nu])