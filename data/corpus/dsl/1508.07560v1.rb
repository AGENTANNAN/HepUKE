# BESIII: absolute branching fractions and hadronic form factors of the semileptonic
# decays D0 -> K- e+ nu_e and D0 -> pi- e+ nu_e.
# Data: 2.92 fb^-1 taken at sqrt(s) = 3.773 GeV (psi(3770) resonance), 2010-2011.
# Technique: double tag -- a single anti-D0 (D0bar) tag is reconstructed in five
# hadronic modes and the semileptonic D0 decay is reconstructed from the recoiling
# tracks; the undetected neutrino is inferred from U_miss = E_miss - |p_miss|.
# The two semileptonic channels are independent final states (different hadron
# hypothesis) -> two TagAnalysis algorithms sharing the same tag side.

### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # psi(3770), 2.92 fb^-1 @ 3.773 GeV
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # inclusive MC (D0D0, D+D-, ditau, nonDD, qq, RR2S, RR1S, QED)

# --- Decay card: psi(3770) -> D0 D0bar, signal D0 -> K- e+ nu_e ---
# The tagged anti-D0 is simulated as a cocktail of the five hadronic tag modes
# weighted by their PDG branching fractions.
decay_card_Kenu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0                          PHSP;
  Enddecay

  Decay D0
  1.0000 K- e+ nu_e                          PHOTOS PHSP;
  Enddecay

  Decay anti-D0
  0.0393 K+ pi-                              PHSP;
  0.1400 K+ pi- pi0                          PHSP;
  0.0803 K+ pi- pi- pi+                      PHSP;
  0.0500 K+ pi- pi+ pi- pi0                  PHSP;
  0.0230 K+ pi- pi0 pi0                      PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                          PHSP;
  Enddecay

  End
DECAYCARD

# --- Decay card: psi(3770) -> D0 D0bar, signal D0 -> pi- e+ nu_e ---
decay_card_pienu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0                          PHSP;
  Enddecay

  Decay D0
  1.0000 pi- e+ nu_e                         PHOTOS PHSP;
  Enddecay

  Decay anti-D0
  0.0393 K+ pi-                              PHSP;
  0.1400 K+ pi- pi0                          PHSP;
  0.0803 K+ pi- pi- pi+                      PHSP;
  0.0500 K+ pi- pi+ pi- pi0                  PHSP;
  0.0230 K+ pi- pi0 pi0                      PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                          PHSP;
  Enddecay

  End
DECAYCARD

exMC_Kenu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "D0toKenu_tag"
  config.related_dataset = psi3770_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_Kenu
  config.cross_section   = :default
end
exMC_Kenu.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_pienu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "D0topienu_tag"
  config.related_dataset = psi3770_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_pienu
  config.cross_section   = :default
end
exMC_pienu.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) — TagAnalysis: single anti-D0 tag + missing neutrino ###

# Five hadronic tag modes of the anti-D0 (charge conjugate of the paper's D0 modes).
tag_modes = [:D0toKPi,          # K+ pi-
             :D0toKPiPi0,       # K+ pi- pi0
             :D0toKPiPiPi,      # K+ pi- pi- pi+
             :D0toKPiPiPiPi0,   # K+ pi- pi+ pi- pi0
             :D0toKPiPi0Pi0]    # K+ pi- pi0 pi0

# ================================================================
# Channel 1: D0 -> K- e+ nu_e
# ================================================================
alg_Kenu = TagAnalysis.new("D0ToKenuDoubleTag")
alg_Kenu.set_header(["D0ToKenuDoubleTagAlg/D0ToKenuDoubleTag.h"])
        .set_constant({ "ECMS" => [:double, 3.773] })      # sqrt(s) = 3.773 GeV
        .with_decay_card(decay_card_Kenu)

alg_Kenu.tag_side(:D0) do |t|
  t.modes(*tag_modes)
  t.charm -1                                               # tagged anti-D0 (D0bar) side
  t.window :mBC, min: 1.858, max: 1.875                    # single D0 tag M_BC signal region
end

# Signal side: the two charged tracks recoiling against the tag, K- and e+, plus nu_e
alg_Kenu.signal_side do |s|
  s.charged(km: 1, ep: 1)                                  # exactly two signal-side tracks
  s.require_charge 0                                       # -1 (K-) + 1 (e+) = 0
  s.missing :nu_e                                          # massless missing neutrino
end

# Kinematic fit: total four-momentum of tag + K- + e+ + nu = CMS energy.
# (The neutrino four-momentum and the U_miss / q^2 observables are auto-stored.)
alg_Kenu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200                                           # loose in BOSS; tight cut applied in ROOT
end

alg_Kenu
  .note(:tag_side_selection,
        "Single D0bar tags: at least two good helix tracks with |cos(theta)| < 0.93 " \
        "originating from the interaction region (|Vr| < 1.0 cm, |Vz| < 15.0 cm). " \
        "Kaon (pion) PID from combined dE/dx + TOF confidence levels: CL_K > CL_pi " \
        "(CL_pi > CL_K) for p < 0.75 GeV/c and CL_K > 0.1% (CL_pi > 0.1%) for " \
        "p >= 0.75 GeV/c. pi0 -> gamma gamma: E > 0.025 (0.050) GeV in the EMC barrel " \
        "(endcap), in time with the beam crossing, photon-to-nearest-track angle > 10 deg, " \
        "1C mass-constrained fit with chi2 < 50. For K+ pi-: |t1 - t2| < 5 ns, opening " \
        "angle < 176 deg and sum of E/p < 1.4 to reject cosmic-ray, Bhabha and dimuon events.")
  .note(:tag_deltaE_windows,
        "Per-mode |Delta E| = |E_Knpi - E_beam| windows, applied within (2-3) sigma of the " \
        "measured E_Knpi resolution: K pi (-0.049,+0.044), K pi pi0 (-0.071,+0.052), " \
        "K pi pi pi (-0.043,+0.043), K pi pi pi pi0 (-0.067,+0.066), K pi pi0 pi0 " \
        "(-0.082,+0.050) GeV. If several Knpi combinations satisfy the energy requirement, " \
        "the one with the smallest |Delta E| is retained (one candidate per mode).")
  .note(:tag_mbc_fit,
        "Single tag yields are extracted from a fit to the M_BC distribution: signal shape " \
        "from simulation convolved with a double Gaussian plus an ARGUS function times a " \
        "third-order polynomial for the combinatorial background. Small wrong-sign peaking " \
        "backgrounds (doubly Cabibbo suppressed decays and K_S0 modes) are estimated from " \
        "MC and subtracted. Total: 2,793,317 +- 3,684 single D0bar tags.")
  .note(:signal_side_selection,
        "Exactly two oppositely charged tracks recoil against the tag, identified as one " \
        "positron and one kaon. Positron: combined CL_e from dE/dx, TOF and EMC (deposited " \
        "energy and shower shape) with CL_e > 0.1% and CL_e/(CL_e + CL_pi + CL_K) > 0.8. " \
        "Kaon: CL_K > CL_pi. Signals are extracted from a fit to the U_miss distribution " \
        "with the CLEO empirical function (Gaussian core with asymmetric power-law tails, " \
        "tail parameters fixed from signal MC and convolved with a Gaussian of free mean " \
        "and width); observed yields 70,727.0 +- 278.3 (K- e+ nu_e).")
  .note(:fsr_recovery,
        "Final-state-radiation recovery: the four-momenta of nearby photons within the " \
        "positron direction are added to the positron four-momentum; the difference between " \
        "the branching fraction with and without FSR recovery is taken as the systematic " \
        "uncertainty (0.30%).")
  .note(:egamma_max_veto,
        "Fake photon background suppressed by requiring the maximum energy of any unused " \
        "photon in the recoil system, E_gamma_max, to be less than 300 MeV (0.10% systematic " \
        "uncertainty from this requirement).")
  .note(:umiss_reconstruction,
        "U_miss = E_miss - |p_miss| with E_miss = E_beam - E_h- - E_e+ and " \
        "p_miss = p_D0 - p_h- - p_e+, where p_D0 = -p_hat_tag * sqrt(E_beam^2 - m_D0^2) " \
        "uses the direction of the single D0bar tag. U_miss peaks at zero for a correctly " \
        "identified semileptonic decay with one missing neutrino.")
  .note(:efficiency_correction,
        "Overall reconstruction efficiencies are corrected for data/MC differences in " \
        "tracking and PID by f_corr = 1.0118 (K- e+ nu_e); corrected efficiency " \
        "eps = 0.7224 +- 0.0012. The absolute branching fraction is obtained as " \
        "B = N_observed / (N_tag * eps) (double-tag technique).")
  .note(:q2_binning_and_efficiency_matrix,
        "Differential decay rates: q^2 = (E_e + E_nu)^2 - (p_e + p_nu)^2 with E_nu = E_miss " \
        "and p_nu = E_miss * p_hat_miss. 18 q^2 bins of 0.1 GeV^2 from 0.0 to 1.7 GeV^2/c^4 " \
        "plus a last bin from 1.7 to q^2_max; an efficiency matrix eps_ij (weighted over the " \
        "five tag modes by their yields) unfolds N_observed^i = sum_j eps_ij N_produced^j, " \
        "and Delta Gamma_i = N_produced^i / (tau_D0 * N_tag).")
  .note(:form_factor_fit,
        "Partial decay rates are fitted with single pole, modified pole (BK), two-parameter " \
        "and three-parameter z-series expansion form-factor models; the main result uses the " \
        "two-parameter series expansion: f_+^K(0)|V_cs| = 0.7172 +- 0.0025 +- 0.0035.")
  .apply

root_files_K = alg_Kenu.execute_on([psi3770_data, psi3770_incMC, exMC_Kenu])

# ================================================================
# Channel 2: D0 -> pi- e+ nu_e
# ================================================================
alg_pienu = TagAnalysis.new("D0ToPienuDoubleTag")
alg_pienu.set_header(["D0ToPienuDoubleTagAlg/D0ToPienuDoubleTag.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })     # sqrt(s) = 3.773 GeV
         .with_decay_card(decay_card_pienu)

alg_pienu.tag_side(:D0) do |t|
  t.modes(*tag_modes)
  t.charm -1                                               # tagged anti-D0 (D0bar) side
  t.window :mBC, min: 1.858, max: 1.875                    # single D0 tag M_BC signal region
end

# Signal side: the two charged tracks recoiling against the tag, pi- and e+, plus nu_e
alg_pienu.signal_side do |s|
  s.charged(pim: 1, ep: 1)                                 # exactly two signal-side tracks
  s.require_charge 0                                       # -1 (pi-) + 1 (e+) = 0
  s.missing :nu_e                                          # massless missing neutrino
end

alg_pienu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_pienu
  .note(:tag_side_selection,
        "Same single D0bar tag selection as the D0 -> K- e+ nu_e channel: five hadronic " \
        "modes with |cos(theta)| < 0.93 tracks, |Vr| < 1.0 cm, |Vz| < 15.0 cm, " \
        "kaon/pion confidence-level PID, pi0 from gamma gamma with a 1C mass-constrained " \
        "fit (chi2 < 50), and the M_BC fit used to extract the tag yields.")
  .note(:tag_deltaE_windows,
        "Per-mode Delta E windows identical to the K- e+ nu_e channel " \
        "(K pi +/-0.04-0.05, K pi pi0 (-0.071,+0.052), K pi pi pi (-0.043,+0.043), " \
        "K pi pi pi pi0 (-0.067,+0.066), K pi pi0 pi0 (-0.082,+0.050) GeV); the candidate " \
        "with the smallest |Delta E| is retained.")
  .note(:signal_side_selection,
        "Exactly two oppositely charged tracks recoil against the tag, identified as one " \
        "positron and one pion. Positron: CL_e > 0.1% and CL_e/(CL_e + CL_pi + CL_K) > 0.8 " \
        "from combined dE/dx, TOF and EMC information. Pion: CL_pi > CL_K. The signal yield " \
        "is extracted from a fit to the U_miss distribution with the CLEO empirical " \
        "function; observed yield 6,297.1 +- 86.8 (pi- e+ nu_e).")
  .note(:fsr_recovery,
        "Final-state-radiation recovery adds nearby photon four-momenta to the positron " \
        "candidate; the resulting 0.30% difference in the branching fraction is taken as " \
        "the systematic uncertainty.")
  .note(:egamma_max_veto,
        "Maximum energy of any unused photon in the recoil system required to be less than " \
        "300 MeV to suppress fake photon background (0.10% systematic uncertainty).")
  .note(:umiss_reconstruction,
        "U_miss = E_miss - |p_miss| computed with the measured tag direction, " \
        "p_D0 = -p_hat_tag * sqrt(E_beam^2 - m_D0^2); it is zero for a correctly " \
        "reconstructed D0 -> pi- e+ nu_e decay with one missing neutrino.")
  .note(:efficiency_correction,
        "Data/MC tracking and PID efficiencies are corrected by f_corr = 0.9814 " \
        "(pi- e+ nu_e); corrected efficiency eps = 0.7643 +- 0.0013, and " \
        "B = N_observed / (N_tag * eps).")
  .note(:q2_binning_and_efficiency_matrix,
        "14 q^2 bins of 0.2 GeV^2 from 0.0 to 2.6 GeV^2/c^4 plus a last bin from 2.6 to " \
        "q^2_max; the measured rates are unfolded with the q^2 efficiency matrix and " \
        "converted to partial widths via Delta Gamma_i = N_produced^i / (tau_D0 * N_tag).")
  .note(:form_factor_fit,
        "The partial decay rates are fitted with the single pole, modified pole (BK) and " \
        "z-series expansion models; the main result uses the two-parameter series " \
        "expansion: f_+^pi(0)|V_cd| = 0.1435 +- 0.0018 +- 0.0009.")
  .apply

root_files_pi = alg_pienu.execute_on([psi3770_data, psi3770_incMC, exMC_pienu])
