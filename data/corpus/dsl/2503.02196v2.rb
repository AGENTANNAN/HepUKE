# First Measurement of the Decay Dynamics in Semileptonic D+(0) → K1(1270) e+ νe
# Paper: 2503.02196v2
# DT method at √s = 3.773 GeV using 712_3773 data.
# Two independent channels (D0 and D+) — two TagAnalysis objects (Rule T1).
# Signal side: semileptonic K-π+π0(-)e+νe with missing νe.

DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

psipp_data  = DatasetManager.real_data.find("712_3773")
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================================
# Channel 1: D0 → K- π+ π- e+ νe  (tagging D0bar)
# ============================================================================

decay_card_d0 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0000 K- pi+ pi- e+ nu_e PHSP;
  Enddecay
  Decay anti-D0
  1.0000 K+ pi-    PHSP;
  1.0000 K+ pi- pi0  PHSP;
  1.0000 K+ pi- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

exMC_d0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_d0_k1enu_d0tag"
  config.related_dataset = psipp_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_d0
  config.cross_section   = :default
end

algo_d0 = TagAnalysis.new("D0toKPiPieNu")
algo_d0
  .set_header(["D0toKPiPieNuAlg/D0toKPiPieNu.h"])
  .set_constant(ECMS: 3.773)

algo_d0.tag_side(:D0) do
  modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  charm -1
end

algo_d0.signal_side do
  charged km: 1, pip: 1, pim: 1, ep: 1
  missing :nu_e
  require_charge 0
end

algo_d0.fit do
  constrain_four_momentum
  chi2_cut 200
end

algo_d0
  .note(:electron_pid,
    "e+ PID: CL_e > 0.001 and CL_e/(CL_e+CL_π+CL_K) > 0.8. " \
    "Hadron suppression for D0: E/p − 0.05 × χ²_(e−dE/dx) > 0.60. " \
    "FSR recovery: neutral showers within 5° of e+ merged into e+ four-momentum.")
  .note(:d0_vetoes,
    "M(K-π+π-π+_(e→π)) < 1.80 GeV/c² (hadronic D0→K-π+π-π+ veto). " \
    "cosθ(e+,π-) < 0.93 (D0→K-π+π0, π0→e+e-γ veto). " \
    "cos(νe,γ_extra) < 0.78 (D0→K-π+π-π+π0 veto). " \
    "M(π+π-) < 0.31 GeV/c² (D0→K-π0e+νe, π0→e+e-γ veto). " \
    "|ΔE[(K+π-)tag π-_sig]| > 7 MeV (tag misreconstruction veto).")
  .note(:d0_no_extra_tracks, "No additional charged tracks allowed beyond the signal candidates (suppresses hadronic background).")
  .note(:umiss_cut, "|U_miss| < 0.03 GeV for signal region, where U_miss ≡ E_miss − |p_miss|.")
  .note(:d0_fit, "Unbinned ML fit to U_miss distribution. Signal shape: MC convolved with Gaussian. " \
    "Peaking background (D0→K-π+π-π+) fixed. Other backgrounds floating.")
  .note(:amplitude_analysis,
    "Amplitude analysis of K-π+π- system fitted with K1(1270)→ρK and K1(1270)→K*π components. " \
    "Amplitude fit, angular analysis for up-down asymmetry A'_ud, and FFs extraction are ROOT-side only.")
  .with_decay_card(decay_card_d0)
  .apply
  .execute_on([psipp_data, psipp_incMC, exMC_d0])

# ============================================================================
# Channel 2: D+ → K- π+ π0 e+ νe  (tagging D-)
# ============================================================================

decay_card_dp = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D- PHSP;
  Enddecay
  Decay D+
  1.0000 K- pi+ pi0 e+ nu_e PHSP;
  Enddecay
  Decay D-
  1.0000 K+ pi- pi-    PHSP;
  1.0000 K_S0 pi-    PHSP;
  1.0000 K+ pi- pi- pi0  PHSP;
  1.0000 K_S0 pi- pi0  PHSP;
  1.0000 K_S0 pi- pi+ pi-  PHSP;
  1.0000 K+ K- pi-    PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

exMC_dp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_dp_k1enu_dptag"
  config.related_dataset = psipp_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_dp
  config.cross_section   = :default
end

algo_dp = TagAnalysis.new("DptoKPiPi0eNu")
algo_dp
  .set_header(["DptoKPiPi0eNuAlg/DptoKPiPi0eNu.h"])
  .set_constant(ECMS: 3.773)

algo_dp.tag_side(:Dplus) do
  modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
        :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  charm -1
end

algo_dp.signal_side do
  photons :pi0
  charged km: 1, pip: 1, ep: 1
  missing :nu_e
  require_charge 1
end

algo_dp.fit do
  constrain_four_momentum
  chi2_cut 200
end

algo_dp
  .note(:electron_pid,
    "e+ PID: CL_e > 0.001 and CL_e/(CL_e+CL_π+CL_K) > 0.8. " \
    "Hadron suppression for D+: E/p − 0.05 × χ²_(e−dE/dx) > 0.53. " \
    "FSR recovery: neutral showers within 5° of e+ merged into e+ four-momentum.")
  .note(:dp_vetoes,
    "p(π0) > 0.2 GeV/c to suppress fake π0. " \
    "M(K-π+π0e+) < 1.78 GeV/c² (hadronic D+→K-π+π0π+ veto). " \
    "U'_miss < 0.03 GeV rejected (D+→K-π+e+νe background with fake π0).")
  .note(:dp_no_extra_tracks, "No additional charged tracks allowed beyond the signal candidates (suppresses hadronic background).")
  .note(:umiss_cut, "|U_miss| < 0.03 GeV for signal region, where U_miss ≡ E_miss − |p_miss|.")
  .note(:dp_fit, "Unbinned ML fit to U_miss distribution. Signal shape: MC convolved with Gaussian. " \
    "All background components floating (no significant peaking background).")
  .note(:dp_d_momentum_constraint,
    "D momentum constrained as p_D = −p̂_Dbar · sqrt(E_beam² − m_Dbar²) to improve resolution.")
  .note(:amplitude_analysis,
    "Amplitude analysis of K-π+π0 system fitted with K1(1270)→ρK and K1(1270)→K*π components. " \
    "Simultaneous fit across D+ and D0 channels imposing isospin relations. " \
    "Angular analysis on cosθ_L vs cosθ_K distributions. All ROOT-side.")
  .with_decay_card(decay_card_dp)
  .apply
  .execute_on([psipp_data, psipp_incMC, exMC_dp])