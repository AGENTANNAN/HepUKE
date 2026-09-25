# First measurement of the absolute branching fractions of Sigma+ nonleptonic decays
# and test of the Delta I = 1/2 rule
#   [arXiv:2512.09628]
#
# Double-tag (DT) analysis at sqrt(s) = 3.097 GeV with (10087 +/- 44) x 10^6 J/psi events.
# ST: Sigma- -> pbar pi0 (reconstructed)
# DT signal modes:
#   Mode 1: Sigma+ -> p pi0  (reconstruct p on signal side, pi0 not reconstructed)
#   Mode 2: Sigma+ -> n pi+  (reconstruct pi+ on signal side, n not reconstructed)
# DT yields extracted from fits to RM_bar_Sigma distributions; BFs from N_DT / (N_ST * epsilon_sig).

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # Inclusive J/psi MC (10B events)

# ---------------------------------------------------------------- decay cards
sub_pi0_gg = <<~DECAYCARD
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
DECAYCARD

# Signal decay card for Mode 1: J/psi -> Sigma+ (-> p pi0) Sigma- (-> anti-p- pi0)
decay_card_sig_mode1 = <<~DECAYCARD
  Decay J/psi
  1.0000 Sigma+ anti-Sigma-  PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0  PHSP;
  Enddecay

  Decay anti-Sigma-
  1.0000 anti-p- pi0  PHSP;
  Enddecay

  #{sub_pi0_gg}

  End
DECAYCARD

# Signal decay card for Mode 2: J/psi -> Sigma+ (-> n0 pi+) Sigma- (-> anti-p- pi0)
decay_card_sig_mode2 = <<~DECAYCARD
  Decay J/psi
  1.0000 Sigma+ anti-Sigma-  PHSP;
  Enddecay

  Decay Sigma+
  1.0000 n0 pi+  PHSP;
  Enddecay

  Decay anti-Sigma-
  1.0000 anti-p- pi0  PHSP;
  Enddecay

  #{sub_pi0_gg}

  End
DECAYCARD

# Exclusive MC for Mode 1: Sigma+ -> p pi0
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Jpsi_SigmaP_pPi0_SigmaM_pbarPi0"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_sig_mode1
  config.cross_section   = :default
end

# Exclusive MC for Mode 2: Sigma+ -> n pi+
exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Jpsi_SigmaP_nPi_SigmaM_pbarPi0"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_sig_mode2
  config.cross_section   = :default
end

### Event selection (BOSS) — Mode 1: Sigma+ -> p pi0 + ST Sigma- -> pbar pi0 ###
alg_mode1 = Algorithm.new("SigmaPPi0_SigmaMPbarPi0")
alg_mode1.set_header(["SigmaPPi0_SigmaMPbarPi0Alg/SigmaPPi0_SigmaMPbarPi0.h"])
          .set_constant({"ECMS" => [:double, 3.097]})

sel_mode1 = Selection.new
sel_mode1.select_track {
            cos_theta 0.93          # |cos(theta)| < 0.93
            Vz        100.0         # |Vz| < 10 cm (paper: < 10 cm along beam)
            Vr        20.0          # |Vr| < 2 cm in transverse plane (paper: < 2 cm)
            nChrp     ">=2"         # at least pbar (anti-p-) and p
            nChrn     ">=1"         # at least pbar negative charge
            nNet      ">=0"         # net charge from proton + anti-proton
          }
         .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    10.0     # > 10 deg from charged tracks
            angle_to_prm_track 20.0   # > 20 deg from anti-proton (to suppress pbar annihilation)
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam              ">=2"    # at least 2 photons for pi0 -> gamma gamma
          }
         .pid(method: :probability) {
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]
            identify :pion, against: [:kaon]
            nprp     ">=1"            # at least one proton (signal Sigma+ -> p pi0)
            nprm     ">=1"            # at least one anti-proton (ST Sigma- -> pbar pi0)
         }
         .remove([:prp <= :chrgp])
         .remove([:prm <= :chrgn])
         .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"
         }

alg_mode1
  .note(:st_selection,
        "ST Sigma- selection: combine anti-proton with pi0 candidates. " \
        "If multiple pbar-pi0 combinations, retain the one with minimum " \
        "|M(pbar pi0) - m_Sigma-|. Require M(pbar pi0) in [1.169, 1.205] GeV/c^2. " \
        "ST yield N_ST = 4504566 +/- 3320 determined by unbinned ML fit to " \
        "RM_bar_Sigma = sqrt((E_cms - E_pbar_pi0)^2/c^4 - P_pbar_pi0^2/c^2) + " \
        "M(pbar pi0) - m_Sigma-. Signal shape: MC convoluted with Gaussian; " \
        "background: J/psi -> Delta+ Delta- MC shape + 2nd-order Chebyshev.")
  .note(:dt_mode1_signal_side,
        "DT Mode 1 (Sigma+ -> p pi0): additionally select a proton on the signal side. " \
        "If both Sigma+ and Sigma- can be reconstructed (double-counting), keep only " \
        "the combination with tagged Sigma closest in mass to m_Sigma+. " \
        "Require RM_bar_Sigma_p = sqrt((E_cms - E_pbar_pi0 - E_p)^2/c^4 - " \
        "(P_pbar_pi0 + P_p)^2/c^2) + M(pbar pi0) - m_Sigma- in " \
        "[0.034, 0.231] GeV/c^2 (3 sigma around pi0 mass). " \
        "DT yield N_DT = 1437863 +/- 1672 from fit to RM_bar_Sigma for the " \
        "Sigma- p sample. DT efficiency epsilon_DT = 26.71%. " \
        "BF = N_DT/(N_ST * epsilon_sig) with epsilon_sig = epsilon_DT/epsilon_ST.")
  .note(:st_efficiency,
        "ST efficiency epsilon_ST = 41.66%. Signal efficiency " \
        "epsilon_sig = epsilon_DT/epsilon_ST. Most ST systematic uncertainties " \
        "cancel in the DT method.")
  .note(:delta_peaks,
        "Peaking backgrounds in DT Mode 1 from other Sigma+ decays with a proton: " \
        "Sigma+ -> p gamma, Sigma+ -> Lambda e+ nu_e, Sigma+ -> n pi+ gamma. " \
        "Systematic uncertainty from these: 0.12%.")
  .note(:signal_shape_sys,
        "Signal shape systematic estimated by varying signal MC shape convoluted " \
        "with Gaussian to signal MC shape only. Changes in BF: 0.18% (mode 1), " \
        "0.26% (mode 2).")
  .note(:background_shape_sys,
        "Background shape systematic: refit with 4th-order Chebyshev polynomial " \
        "for all background components. Changes: 0.30% (mode 1), 0.26% (mode 2).")
  .note(:mc_generator_sys,
        "MC generator systematic: vary decay parameters by +/- 1 sigma of measured " \
        "values [Refs. 15, 24]. Efficiency changes: 0.05% (mode 1), 0.04% (mode 2).")
  .note(:rm_window_sys,
        "RM window systematic: change from +/- 3 sigma to +/- 2.5, 3.5, or 4 sigma. " \
        "Largest deviation: 0.14% (mode 1), 0.16% (mode 2).")
  .note(:deltaI_half_rule,
        "Delta I = 1/2 rule tested via amplitudes: " \
        "A(Sigma- -> n pi-) - A(Sigma+ -> n pi+) + sqrt(2) A(Sigma+ -> p pi0) " \
        "= -0.127 +/- 0.014 (this work) vs -0.180 +/- 0.063 (PDG). " \
        "B amplitudes: -2.78 +/- 0.16 (this work) vs -2.50 +/- 0.76 (PDG). " \
        "Both deviate from zero by > 5 sigma, indicating Delta I = 3/2 amplitude.")
  .note(:bf_results,
        "BF(Sigma+ -> p pi0) = (49.79 +/- 0.06 +/- 0.22)%, " \
        "BF(Sigma+ -> n pi+) = (49.87 +/- 0.05 +/- 0.29)%. " \
        "Ratio = 0.9984 +/- 0.0016 +/- 0.0073. " \
        "Deviations from PDG: 4.4 sigma for p pi0, 3.4 sigma for n pi+, " \
        "5.5 sigma for ratio. Total systematic: 0.44% (mode 1), 0.59% (mode 2).")
  .with_decay_card(decay_card_sig_mode1)
  .apply(sel_mode1)

### Event selection (BOSS) — Mode 2: Sigma+ -> n pi+ + ST Sigma- -> pbar pi0 ###
alg_mode2 = Algorithm.new("SigmaPNPi_SigmaMPbarPi0")
alg_mode2.set_header(["SigmaPNPi_SigmaMPbarPi0Alg/SigmaPNPi_SigmaMPbarPi0.h"])
          .set_constant({"ECMS" => [:double, 3.097]})

sel_mode2 = Selection.new
sel_mode2.select_track {
            cos_theta 0.93
            Vz        100.0
            Vr        20.0
            nChrp     ">=1"         # at least one positive track (pi+ on signal side)
            nChrn     ">=1"         # at least one negative track (pbar for ST)
            nNet      "==0"         # pi+ + pbar sum to 0
          }
         .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    10.0
            angle_to_prm_track 20.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam              ">=2"
          }
         .pid(method: :probability) {
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]
            identify :pion, against: [:kaon]
            nprm     ">=1"            # anti-proton for ST
            npip     ">=1"            # pi+ for signal side
         }
         .remove([:prm <= :chrgn])
         .remove([:pip <= :chrgp])
         .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"
         }

alg_mode2
  .note(:st_selection,
        "ST Sigma- selection: identical to Mode 1. " \
        "Combine anti-proton with pi0 candidates, retain minimum |M(pbar pi0) - m_Sigma-|. " \
        "Require M(pbar pi0) in [1.169, 1.205] GeV/c^2.")
  .note(:dt_mode2_signal_side,
        "DT Mode 2 (Sigma+ -> n pi+): additionally select a pi+ on the signal side. " \
        "Require RM_bar_Sigma_pi+ = sqrt((E_cms - E_pbar_pi0 - E_pi)^2/c^4 - " \
        "(P_pbar_pi0 + P_pi)^2/c^2) + M(pbar pi0) - m_Sigma- in " \
        "[0.881, 0.997] GeV/c^2 (3 sigma around neutron mass). " \
        "DT yield N_DT = 1638422 +/- 1780 from fit to RM_bar_Sigma for the " \
        "Sigma- pi+ sample. DT efficiency epsilon_DT = 30.38%.")
  .note(:bf_results,
        "BF(Sigma+ -> p pi0) = (49.79 +/- 0.06 +/- 0.22)%, " \
        "BF(Sigma+ -> n pi+) = (49.87 +/- 0.05 +/- 0.29)%. " \
        "Total systematic: 0.44% (mode 1), 0.59% (mode 2).")
  .note(:deltaI_half_rule,
        "Delta I = 1/2 rule test: both A and B amplitude combinations deviate from zero " \
        "by > 5 sigma, providing first compelling evidence for non-negligible " \
        "Delta I = 3/2 transition amplitude in Sigma hyperon decays.")
  .with_decay_card(decay_card_sig_mode2)
  .apply(sel_mode2)

alg_mode1.execute_on([jpsi_data, jpsi_incMC, exMC_mode1])
alg_mode2.execute_on([jpsi_data, jpsi_incMC, exMC_mode2])