# BESIII paper 1106.5118v1 — Search for CP and P violating pseudoscalar decays into pi pi
# J/psi -> gamma eta / eta' / eta_c, with eta/eta'/eta_c -> pi+ pi- or pi0 pi0.
# Sample: (225.2 +- 2.8) x 10^6 J/psi events collected with BESIII.
# Two final-state topologies are searched:
#   Charged mode : gamma pi+ pi-   (eta/eta'/eta_c -> pi+ pi-)
#   Neutral mode : gamma pi0 pi0   (eta/eta'/eta_c -> pi0 pi0)
# The three pseudoscalars share the same topology within each mode, so one
# Algorithm per topology is used; the eta_c-specific EMC/dE/dx requirements are
# recorded as notes.

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi data, 225.2e6 events
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # inclusive J/psi MC

# Charged mode: J/psi -> gamma R, R -> pi+ pi-
decay_card_gam_pipi = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta                 PHSP;
    Enddecay

    Decay eta
    1.0000 pi+ pi-                   PHSP;
    Enddecay

    End
DECAYCARD

# Neutral mode: J/psi -> gamma R, R -> pi0 pi0
decay_card_gam_pi0pi0 = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta                 PHSP;
    Enddecay

    Decay eta
    1.0000 pi0 pi0                   PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma              PHSP;
    Enddecay

    End
DECAYCARD

exMC_pipi = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "jpsi_gamma_eta_pipim"
    config.related_dataset = jpsi_data
    config.events          = 200000
    config.decay_card      = decay_card_gam_pipi
    config.cross_section   = :default
end

exMC_pi0pi0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "jpsi_gamma_eta_pi0pi0"
    config.related_dataset = jpsi_data
    config.events          = 200000
    config.decay_card      = decay_card_gam_pi0pi0
    config.cross_section   = :default
end

### Event selection (BOSS) — Charged mode: J/psi -> gamma pi+ pi- ###
alg_name_pipi = "JpsiGammaPiPi"
alg_pipi = Algorithm.new(alg_name_pipi)
alg_pipi.set_header(["#{alg_name_pipi}Alg/#{alg_name_pipi}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})

sel_pipi = Selection.new
sel_pipi.select_track {
           cos_theta 0.93   # |cos(theta)| < 0.93
           Vz        10.0   # within +-10 cm of the IP along the beam direction
           Vr         1.0   # within +-1 cm of the beamline in the transverse plane
           nChrp    "==1"   # two charged tracks with zero net charge
           nChrn    "==1"
           nNet     "==0"
         }
         .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14    # EMC cluster timing suppresses electronic noise
           angle_to_track    20.0  # photon separated by >= 20 deg from any charged track
           energyThreshold_b 0.025 # barrel (|cos theta| < 0.8) : E > 25 MeV
           energyThreshold_e 0.050 # endcap (0.86 < |cos theta| < 0.92) : E > 50 MeV
           nGam              ">=1" # at least one photon; the most energetic is the radiative photon
         }
         .pid(method: :probability) {
           prob_cut 0.001
           # For pion candidates Prob_PID(pi) > Prob_PID(K) and Prob_PID(pi) > 0.001;
           # at least one charged track must be identified as a pion.
           identify :pion, against: [:kaon, :proton]
           npip ">=1"
           npim ">=1"
         }
         # 4C kinematic fit under energy-momentum conservation for the signal hypothesis
         # J/psi -> gamma pi+ pi-. The photon with the maximum energy in the e+e- C.M. frame
         # is taken as the radiative photon; the combination with the smallest chi2 is kept.
         .kinematic_fit([:gamma, :pip, :pim]) {
           nominal
           constrain_four_momentum
           chi2_cut 200   # loose BOSS cut; the paper requires chi2_4C(gamma pi+pi-) < 30 (applied in ROOT)
         }
         # Competing hypothesis J/psi -> gamma K+ K-: the chi2 is stored for the ROOT-level
         # requirement chi2_4C(gamma pi+pi-) < chi2_4C(gamma K+K-).
         .kinematic_fit([:gamma, :pip, :pim]) {
           constrain_four_momentum
         }

alg_pipi.note(:pion_identification,
              "Pion candidates require Prob_PID(pi) > Prob_PID(K) and Prob_PID(pi) > 0.001, where " \
              "Prob_PID(i) is the confidence level from combined TOF and dE/dx information. " \
              "Only one of the two pions is required to be identified.")
         .note(:competing_hypothesis_veto,
              "Candidate events are also fitted to J/psi -> gamma K+ K-; the signal hypothesis " \
              "J/psi -> gamma pi+ pi- must satisfy chi2_4C(gamma pi+pi-) < chi2_4C(gamma K+K-). " \
              "The competing chi2 is stored by the second kinematic fit; the comparison is applied " \
              "in ROOT.")
         .note(:pi0_veto,
              "To suppress the dominant J/psi -> rho pi -> pi+ pi- pi0 background, for events with " \
              "two or more photons all pairings of the radiative photon with the remaining photons " \
              "are used to form pi0 candidates and the pairing closest to the nominal pi0 mass is " \
              "selected; the event is rejected if 0.12 < m(gamma gamma_rad) < 0.15 GeV/c^2. " \
              "Efficiency 99% for eta(eta') -> pi+pi- and 96% for eta_c -> pi+pi-. Applied in ROOT.")
         .note(:electron_background,
              "J/psi -> e+ e- background suppressed by E_pi^EMC < 1.2 GeV and " \
              "|chi_dE/dx(pi)| < 3 for both pi+ and pi- (E_pi^EMC is the EMC deposited energy). " \
              "Applied in ROOT.")
         .note(:muon_background,
              "For the eta_c -> pi+ pi- search, the J/psi -> mu+ mu- background is suppressed by " \
              "requiring E_pi+^EMC > 0.4 GeV or E_pi-^EMC > 0.4 GeV (removes 99.96% of the " \
              "J/psi -> mu+mu- background, efficiency 62.4% for eta_c -> pi+pi-). This requirement " \
              "is not applied for eta/eta' -> pi+pi- where the contamination is very small. " \
              "Applied in ROOT.")
         .note(:signal_regions,
              "Blinded signal regions in m(pi+pi-): 0.53-0.56 GeV/c^2 (eta -> pi+pi-), " \
              "0.95-0.97 GeV/c^2 (eta' -> pi+pi-) and 2.95-3.02 GeV/c^2 (eta_c -> pi+pi-).")
         .note(:background_study,
              "Backgrounds studied with an inclusive J/psi MC sample and dedicated exclusive MC. " \
              "Main non-peaking backgrounds: J/psi -> rho pi, J/psi -> mu+mu-, J/psi -> e+e-, " \
              "J/psi -> a_2(1320) pi -> gamma pi+pi-, J/psi -> b_1(1235) pi -> gamma pi+pi-, " \
              "J/psi -> pi+pi-, J/psi -> gamma sigma/f_2(1270)/f_0(1500)/f_0(1710) -> gamma pi+pi- " \
              "and ISR e+e- -> gamma_ISR pi+pi-. Possible peaking backgrounds: " \
              "J/psi -> gamma eta with eta -> gamma pi+pi- (less than 11 events in +-2 sigma around " \
              "the eta peak, neglected), J/psi -> gamma eta' with eta' -> gamma rho0 -> gamma pi+pi- " \
              "(fixed in the fit), and the OZI-suppressed J/psi -> gamma eta_c with " \
              "eta_c -> gamma pi+pi- (negligible).")
         .note(:mass_spectrum_fit,
              "The pi+ pi- invariant mass is fitted with an MC signal shape plus a 2nd-order " \
              "Chebyshev background (for eta) or with the MC eta' signal shape, the normalised " \
              "J/psi -> gamma eta' -> gamma gamma pi+pi- peaking background and a 2nd-order " \
              "Chebyshev (for eta'), and with an acceptance-corrected MC eta_c signal shape plus a " \
              "3rd-order Chebyshev (for eta_c). Fitted yields: 17 +- 23 (eta, 0.8 sigma), " \
              "0.1 +- 15 (eta', 0.1 sigma), 52 +- 35 (eta_c, 1.5 sigma). ROOT-level.")
         .note(:upper_limits,
              "Bayesian 90% C.L. upper limits on the signal yields: N_sig^UP = 48 (eta -> pi+pi-), " \
              "32 (eta' -> pi+pi-) and 92 (eta_c -> pi+pi-). Quoted upper limits on the branching " \
              "fractions: 3.9e-4, 5.5e-5 and 1.3e-4, respectively (efficiencies 54.28%, 53.81% and " \
              "25.27%, with the efficiencies lowered by 1 - sigma_sys for conservatism).")
         .note(:systematics,
              "Tracking 2% per track (4% total), kinematic fit 2%, photon efficiency 1% per photon, " \
              "MC statistics 1.4%, EMC energy and |chi_dE/dx| cuts 0.2%, background shape 4.2% / " \
              "6.3% / 6.6% for eta / eta' / eta_c, PDG branching fractions 3.1% / 2.9% / 24%, " \
              "resonance parameter 8.0% (eta_c only), N_J/psi 1.3%; total in quadrature " \
              "7.3% / 8.6% / 27%. Pion identification uncertainty neglected since only one pion is " \
              "required.")

alg_pipi.with_decay_card(decay_card_gam_pipi).apply(sel_pipi)
root_files_pipi = alg_pipi.execute_on([jpsi_data, jpsi_incMC, exMC_pipi])

### Event selection (BOSS) — Neutral mode: J/psi -> gamma pi0 pi0 ###
alg_name_pi0pi0 = "JpsiGammaPi0Pi0"
alg_pi0pi0 = Algorithm.new(alg_name_pi0pi0)
alg_pi0pi0.set_header(["#{alg_name_pi0pi0}Alg/#{alg_name_pi0pi0}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_pi0pi0 = Selection.new
sel_pi0pi0.select_track {
            cos_theta 0.93
            Vz        10.0
            Vr         1.0
            nTot     "==0"   # no charged tracks (neutral final state)
            nNet     "==0"
          }
          .select_photon {
            # no EMC cluster timing requirement in this mode
            energyThreshold_b 0.025  # barrel : E > 25 MeV
            energyThreshold_e 0.050  # endcap : E > 50 MeV
            nGam              ">=5"  # five or six photons
          }
          # Two pi0 candidates are formed from the photons; the combination with the minimum
          # chi = sqrt((m_gg1 - m_pi0)^2 + (m_gg2 - m_pi0)^2) is selected, realised by the
          # chi2-minimising combination choice of two Kalman pi0 reconstructions.
          .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=2"
          }
          # 4C kinematic fit under the hypothesis J/psi -> gamma pi0 pi0, which improves the
          # pi0 pi0 mass resolution.
          .kinematic_fit([:gamma, :pi0, :pi0]) {
            nominal
            constrain_four_momentum
            chi2_cut 200   # default loose requirement; the analysis uses the very loose default 200
          }

alg_pi0pi0.note(:neutral_topology,
               "An event must have 5 or 6 photons and no charged tracks. The EMC cluster timing " \
               "requirement is not used in photon selection for this mode.")
         .note(:radiative_photon_and_pairing,
               "For eta/eta' -> pi0 pi0 the photon with the maximum energy is taken as the " \
               "radiative photon and all remaining photons form the two pi0 candidates. For " \
               "eta_c -> pi0 pi0, photons with E < 0.3 GeV (potential radiative photons) are " \
               "tried; if more than one radiative-photon candidate exists the one giving the " \
               "smallest |m(5 gamma) - m(J/psi)| is used. The chosen pairing is that minimising " \
               "chi = sqrt((m_gg1 - m_pi0)^2 + (m_gg2 - m_pi0)^2). Endcap-endcap photon " \
               "combinations (about 0.7% of all pi0 candidates) are removed. Applied in ROOT.")
         .note(:pi0_mass_windows,
               "The gamma gamma invariant masses must satisfy |m_gg - m_pi0| < 0.01625 GeV/c^2 " \
               "(~2.5 sigma) for eta/eta' -> pi0 pi0 and |m_gg - m_pi0| < 0.0175 GeV/c^2 " \
               "(~2.5 sigma) for eta_c -> pi0 pi0. Applied in ROOT.")
         .note(:omega_veto,
               "Backgrounds with omega -> gamma pi0 are suppressed by rejecting events with " \
               "0.72 < m(gamma pi0) < 0.82 GeV/c^2. Applied in ROOT.")
         .note(:signal_regions,
               "Blinded signal regions in m(pi0 pi0): 0.52-0.57 GeV/c^2 (eta -> pi0 pi0), " \
               "0.93-0.98 GeV/c^2 (eta' -> pi0 pi0) and 2.95-3.02 GeV/c^2 (eta_c -> pi0 pi0).")
         .note(:background_study,
               "Main backgrounds with pi0 signals: J/psi -> a_2(1320) pi -> gamma pi0 pi0, " \
               "J/psi -> b_1(1235) pi -> gamma pi0 pi0, J/psi -> omega pi0 (pi0) with " \
               "omega -> gamma pi0, and J/psi -> gamma sigma / f_2(1270) / f_0(1500) / " \
               "f_0(1710) / f_4(2050) -> gamma pi0 pi0. Non-pi0 background is represented by the " \
               "normalised number of events in the two pi0 mass sidebands " \
               "0.084 < m_gg < 0.098 GeV/c^2 or 0.168 < m_gg < 0.182 GeV/c^2.")
         .note(:mass_spectrum_fit,
               "Since there are no peaking backgrounds, the pi0 pi0 invariant mass distributions " \
               "are fitted with eta, eta' and acceptance-corrected eta_c MC signal shapes plus " \
               "2nd-order Chebyshev backgrounds. Fitted yields: 11 +- 18 (eta, 0.6 sigma), " \
               "75 +- 30 (eta', 2.6 sigma), 0.1 +- 14 (eta_c, 0.1 sigma). The small accumulation " \
               "in the eta' signal region may come from J/psi -> gamma f_0(980) -> gamma pi0 pi0 " \
               "but is not included in the fit; the peak may be a statistical fluctuation. " \
               "ROOT-level.")
         .note(:upper_limits,
               "Bayesian 90% C.L. upper limits on the signal yields: N_sig^UP = 36 (eta -> pi0pi0), " \
               "110 (eta' -> pi0pi0) and 40 (eta_c -> pi0pi0). Quoted upper limits on the " \
               "branching fractions: 6.9e-4, 4.5e-4 and 4.2e-5, respectively (efficiencies 23.75%, " \
               "23.18% and 35.70%, with the efficiencies lowered by 1 - sigma_sys for " \
               "conservatism); for R -> pi0 pi0 the branching fraction is divided by " \
               "B(pi0 -> gamma gamma)^2.")
         .note(:systematics,
               "Photon efficiency 1% per photon (5% for the neutral mode), pi0 selection 2.0%, " \
               "MC statistics 1.4%, trigger efficiency <0.1%, background shape 5.6% / 5.5% / 8.8% " \
               "for eta / eta' / eta_c, PDG branching fractions 3.1% / 2.9% / 24%, resonance " \
               "parameter 9.5% (eta_c only), N_J/psi 1.3%; total in quadrature 8.6% / 8.5% / 28%. " \
               "The kinematic-fit uncertainty is neglected because only the very loose default " \
               "chi2 requirement (200) is used.")

alg_pi0pi0.with_decay_card(decay_card_gam_pi0pi0).apply(sel_pi0pi0)
root_files_pi0pi0 = alg_pi0pi0.execute_on([jpsi_data, jpsi_incMC, exMC_pi0pi0])
