# =============================================================================
# BESIII: Search for eta_c' -> rho0 rho0, K*0 K*0bar and phi phi
# in psi' -> gamma eta_c'  (psi' -> gamma V V)
# arXiv:1110.0949v2,  1.06 x 10^8 psi' events
# =============================================================================

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# ---- Decay cards (three exclusive final states) ------------------------------
# psi' -> gamma eta_c', eta_c' -> rho0 rho0 -> 2(pi+ pi-)
decay_card_rho = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta_c(2S)   PHSP;
  Enddecay

  Decay eta_c(2S)
  1.000 rho0 rho0         SVV_0;
  Enddecay

  Decay rho0
  1.000 pi+ pi-           VSS;
  Enddecay

  End
DECAYCARD

# psi' -> gamma eta_c', eta_c' -> K*0 K*0bar -> pi+ pi- K+ K-
decay_card_kstar = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta_c(2S)      PHSP;
  Enddecay

  Decay eta_c(2S)
  1.000 K*0 anti-K*0         SVV_0;
  Enddecay

  Decay K*0
  1.000 K+ pi-               VSS;
  Enddecay

  Decay anti-K*0
  1.000 K- pi+               VSS;
  Enddecay

  End
DECAYCARD

# psi' -> gamma eta_c', eta_c' -> phi phi -> 2(K+ K-)
decay_card_phi = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta_c(2S)   PHSP;
  Enddecay

  Decay eta_c(2S)
  1.000 phi phi           SVV_0;
  Enddecay

  Decay phi
  1.000 K+ K-             VSS;
  Enddecay

  End
DECAYCARD

# ---- Exclusive MC samples ----------------------------------------------------
exMC_rho = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psip_gamma_etacp_rho0rho0"
  c.related_dataset = psip_data
  c.events          = 200000
  c.decay_card      = decay_card_rho
  c.cross_section   = :default
end

exMC_kstar = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psip_gamma_etacp_Kstar0Kstar0bar"
  c.related_dataset = psip_data
  c.events          = 200000
  c.decay_card      = decay_card_kstar
  c.cross_section   = :default
end

exMC_phi = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psip_gamma_etacp_phiphi"
  c.related_dataset = psip_data
  c.events          = 200000
  c.decay_card      = decay_card_phi
  c.cross_section   = :default
end

# =============================================================================
# ALGORITHM 1: psi' -> gamma eta_c', eta_c' -> rho0 rho0 -> gamma 2(pi+ pi-)
# =============================================================================
alg_rho = Algorithm.new("PsipGammaEtacpRhoRho")
alg_rho.set_header(["PsipGammaEtacpRhoRhoAlg/PsipGammaEtacpRhoRho.h"])
       .set_constant({"ECMS" => [:double, 3.686]})   # psi(2S) center-of-mass energy (GeV)

sel_rho = Selection.new
sel_rho.select_track {
          cos_theta  0.93    # |cos(theta)| < 0.93 with respect to the e+ beam direction
          Vz         10.0    # within 10 cm of the IP along the beam axis
          Vr         1.0     # within 1 cm of the IP transverse to the beam line
          nChrp      "==2"   # four charged tracks, zero net charge
          nChrn      "==2"
          nNet       "==0"
        }
        .select_photon {
          # EMC clusters with E > 25 MeV in the barrel (|cos|<0.8) / endcap
          # (0.86<|cos|<0.92) active area; timing suppresses noise
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    10.0
          energyThreshold_b 0.025
          energyThreshold_e 0.025
          nGam              ">=1"   # one radiative photon
        }
        # dE/dx + TOF information; tracks matched to the decay-channel particles by
        # minimizing the summed chi2_PID (combinatorial assignment, pi/K/p hypotheses)
        .pid(method: :chi2_sum) {
          chi_min_cut 4
          identify :pion, :kaon, :proton
        }
        # 4C kinematic fit: the four-momenta of the four charged tracks and the photon
        # candidate are constrained to the initial psi' four-momentum. The photon giving
        # the smallest chi2_4C is taken as the radiative photon; chi2_4C < 40.
        .kinematic_fit([:pip, :pim, :pip, :pim, :gamma]) {
          constrain_four_momentum
          chi2_cut 40
        }
        # 3C kinematic fit (fitted photon energy not used) provides the final mass
        # spectrum M_X^3C, which separates the eta_c' signal from psi' -> X background
        .kinematic_fit([:pip, :pim, :pip, :pim, :gamma]) {
          nominal
          constrain_three_momentum
        }

alg_rho
  .note(:track_quality,
        "each charged track is required to have good quality in the track fit, in addition to " \
        "the |cos(theta)| < 0.93, Vz and Vr requirements.")
  .note(:pid_channel_assignment,
        "chi2_PID(i) is computed for each charged track under the pion, kaon and proton hypotheses; " \
        "the total chi2_PID is the sum over the tracks and the track-to-final-state matching with the " \
        "minimum total chi2_PID is adopted. The decay channel of a reconstructed event is selected as " \
        "the one with the minimum chi2_PID among the possible decay channels.")
  .note(:brem_recovery,
        "final-state radiation and bremsstrahlung recovery are applied to the charged tracks.")
  .note(:v_mass_windows,
        "background from non-VV production is reduced by requiring 0.67 < M(pi+pi-) < 0.87 GeV/c^2 " \
        "for the rho0 candidates (and 0.85 < M(pi+-K-+) < 0.95 GeV/c^2 for K*0, " \
        "1.00 < M(K+K-) < 1.03 GeV/c^2 for phi in the other channels); the windows are determined by " \
        "fitting the mass distributions in the chi_cJ mass region and are applied at the ROOT level.")
  .note(:jpsi_veto,
        "background from psi' -> pi+pi- J/psi (J/psi -> lepton pair) and psi' -> eta J/psi is removed " \
        "by requiring the recoil mass of any pi+pi- pair to satisfy M_recoil(pi+pi-) < 3.05 GeV/c^2.")
  .note(:pi0x_background,
        "the psi' -> pi0 X background is measured from data by reconstructing pi0 -> gamma gamma; when " \
        "more than two photons are present the pi0 candidate with the minimum chi2 from a 5C fit " \
        "(4C plus a pi0 mass constraint) is chosen, and chi2_5C < 30 is required as a veto. The shape is " \
        "described by a Novosibirsk function and fixed in the final fit.")
  .note(:fsr_background,
        "the psi' -> (gamma_FSR) X background shape is taken from MC simulation with the FSR photon " \
        "generated by PHOTOS and scaled by f_FSR = 1.70 +- 0.10 (X = 2(pi+pi-)) measured from the " \
        "psi' -> gamma chi_c0, chi_c0 -> (gamma_FSR) X control sample.")
  .note(:continuum_background,
        "the e+e- -> gamma* -> (gamma_FSR) X continuum and ISR backgrounds are estimated from a " \
        "923 pb^-1 data sample at sqrt(s) = 3.773 GeV (sample 712_3773) with luminosity and " \
        "cross-section energy-dependence normalisation; 46 +- 3 events are expected for V = rho0.")
  .note(:signal_extraction,
        "the signal yield is extracted from an unbinned maximum-likelihood fit to the M_VV^3C " \
        "distribution, using a signal shape BW(m0,Gamma) x E_gamma^3 x damping convoluted with a " \
        "Gaussian resolution function; this fit is performed at the ROOT level.")

alg_rho.with_decay_card(decay_card_rho).apply(sel_rho)
root_files_rho = alg_rho.execute_on([psip_data, psip_incMC, exMC_rho])

# =============================================================================
# ALGORITHM 2: psi' -> gamma eta_c', eta_c' -> K*0 K*0bar -> gamma pi+ pi- K+ K-
# =============================================================================
alg_kstar = Algorithm.new("PsipGammaEtacpKstarKstar")
alg_kstar.set_header(["PsipGammaEtacpKstarKstarAlg/PsipGammaEtacpKstarKstar.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_kstar = Selection.new
sel_kstar.select_track {
            cos_theta  0.93
            Vz         10.0
            Vr         1.0
            nChrp      "==2"
            nChrn      "==2"
            nNet       "==0"
          }
          .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.025
            nGam              ">=1"
          }
          .pid(method: :chi2_sum) {
            chi_min_cut 4
            identify :pion, :kaon, :proton
          }
          .kinematic_fit([:pip, :pim, :kp, :km, :gamma]) {
            constrain_four_momentum
            chi2_cut 40
          }
          .kinematic_fit([:pip, :pim, :kp, :km, :gamma]) {
            nominal
            constrain_three_momentum
          }

alg_kstar
  .note(:track_quality,
        "each charged track is required to have good quality in the track fit in addition to the " \
        "|cos(theta)| < 0.93, Vz and Vr requirements.")
  .note(:pid_channel_assignment,
        "chi2_PID(i) is computed for each charged track under pion, kaon and proton hypotheses; the " \
        "track-to-final-state matching with the minimum total chi2_PID is adopted and the decay channel " \
        "with the minimum chi2_PID among the possible channels is selected.")
  .note(:brem_recovery,
        "final-state radiation and bremsstrahlung recovery are applied to the charged tracks.")
  .note(:v_mass_windows,
        "K*0 candidates are required to satisfy 0.85 < M(pi+-K-+) < 0.95 GeV/c^2 (in addition to " \
        "0.67 < M(pi+pi-) < 0.87 GeV/c^2 and 1.00 < M(K+K-) < 1.03 GeV/c^2 used in the other " \
        "channels); the windows are determined by fitting the mass distributions in the chi_cJ region " \
        "and are applied at the ROOT level.")
  .note(:jpsi_veto,
        "background from psi' -> pi+pi- J/psi and psi' -> eta J/psi is removed by requiring the recoil " \
        "mass of any pi+pi- pair to satisfy M_recoil(pi+pi-) < 3.05 GeV/c^2.")
  .note(:pi0x_background,
        "the psi' -> pi0 X background is measured from data with the pi0 reconstructed from two photons " \
        "and chi2_5C < 30 from the 5C fit (4C plus pi0 mass constraint); the shape is fitted with a " \
        "Novosibirsk function and fixed in the final fit.")
  .note(:fsr_background,
        "the psi' -> (gamma_FSR) X background shape is taken from MC simulation with PHOTOS and scaled " \
        "by the measured factor f_FSR = 1.39 +- 0.08 for X = pi+ pi- K+ K-.")
  .note(:continuum_background,
        "continuum and ISR backgrounds are estimated from the 923 pb^-1 data sample at sqrt(s) = " \
        "3.773 GeV; 8 +- 2 events are expected for V = K*0.")
  .note(:signal_extraction,
        "the yield is extracted from an unbinned maximum-likelihood fit to M_VV^3C using the same " \
        "signal parametrisation as the rho0 rho0 channel; performed at the ROOT level.")

alg_kstar.with_decay_card(decay_card_kstar).apply(sel_kstar)
root_files_kstar = alg_kstar.execute_on([psip_data, psip_incMC, exMC_kstar])

# =============================================================================
# ALGORITHM 3: psi' -> gamma eta_c', eta_c' -> phi phi -> gamma 2(K+ K-)
# =============================================================================
alg_phi = Algorithm.new("PsipGammaEtacpPhiPhi")
alg_phi.set_header(["PsipGammaEtacpPhiPhiAlg/PsipGammaEtacpPhiPhi.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

sel_phi = Selection.new
sel_phi.select_track {
          cos_theta  0.93
          Vz         10.0
          Vr         1.0
          nChrp      "==2"
          nChrn      "==2"
          nNet       "==0"
        }
        .select_photon {
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    10.0
          energyThreshold_b 0.025
          energyThreshold_e 0.025
          nGam              ">=1"
        }
        .pid(method: :chi2_sum) {
          chi_min_cut 4
          identify :pion, :kaon, :proton
        }
        .kinematic_fit([:kp, :km, :kp, :km, :gamma]) {
          constrain_four_momentum
          chi2_cut 40
        }
        .kinematic_fit([:kp, :km, :kp, :km, :gamma]) {
          nominal
          constrain_three_momentum
        }

alg_phi
  .note(:track_quality,
        "each charged track is required to have good quality in the track fit in addition to the " \
        "|cos(theta)| < 0.93, Vz and Vr requirements.")
  .note(:pid_channel_assignment,
        "chi2_PID(i) is computed for each charged track under pion, kaon and proton hypotheses; the " \
        "track-to-final-state matching with the minimum total chi2_PID is adopted and the decay channel " \
        "with the minimum chi2_PID among the possible channels is selected.")
  .note(:brem_recovery,
        "final-state radiation and bremsstrahlung recovery are applied to the charged tracks.")
  .note(:v_mass_windows,
        "phi candidates are required to satisfy 1.00 < M(K+K-) < 1.03 GeV/c^2; the window is determined " \
        "by fitting the mass distribution in the chi_cJ mass region and is applied at the ROOT level.")
  .note(:jpsi_veto,
        "background from psi' -> pi+pi- J/psi and psi' -> eta J/psi is removed by requiring the recoil " \
        "mass of any pi+pi- pair to satisfy M_recoil(pi+pi-) < 3.05 GeV/c^2.")
  .note(:continuum_background,
        "no events from the 923 pb^-1 data sample at sqrt(s) = 3.773 GeV survive the phi phi selection.")
  .note(:signal_extraction,
        "only one eta_c' -> phi phi candidate event is found in the signal region; no fit is performed " \
        "and a Poisson PDF is used for the upper-limit calculation (ROOT level).")

alg_phi.with_decay_card(decay_card_phi).apply(sel_phi)
root_files_phi = alg_phi.execute_on([psip_data, psip_incMC, exMC_phi])
