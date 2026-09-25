# BESIII paper 1103.2661v1 — chi_cJ -> p pbar K+ K- (J = 0,1,2)
# psi' -> gamma chi_cJ -> gamma p pbar K+ K-, 106 million psi' events.
# Intermediate states studied: Lambda(1520) -> p K-, phi -> K+ K-, rho etc.
# The three chi_cJ states share the identical final state
# (gamma p pbar K+ K-) and the identical event selection, so a single
# Algorithm object covers all three.

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) data, 106e6 events
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # 100e6 inclusive psi' MC

# Signal decay cards: psi' -> gamma chi_cJ -> gamma p pbar K+ K-
decay_card_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0             PHSP;
    Enddecay

    Decay chi_c0
    1.0000 p+ anti-p- K+ K-         PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1             PHSP;
    Enddecay

    Decay chi_c1
    1.0000 p+ anti-p- K+ K-         PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2             PHSP;
    Enddecay

    Decay chi_c2
    1.0000 p+ anti-p- K+ K-         PHSP;
    Enddecay

    End
DECAYCARD

# Detection efficiencies are determined separately for the three chi_cJ states
exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic0_ppbarKpKm"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_chic0
    config.cross_section   = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic1_ppbarKpKm"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_chic1
    config.cross_section   = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic2_ppbarKpKm"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_chic2
    config.cross_section   = :default
end

# Exclusive MC for the intermediate / resonant sub-channels
decay_card_lambda1520 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1             PHSP;
    Enddecay

    Decay chi_c1
    1.0000 anti-p- K+ Lambda0       PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ K-                    PHSP;
    Enddecay

    End
DECAYCARD

decay_card_pphip = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1             PHSP;
    Enddecay

    Decay chi_c1
    1.0000 p+ anti-p- phi           PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-                    VSS;
    Enddecay

    End
DECAYCARD

exMC_lambda1520 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic1_pbarKpLambda1520"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_lambda1520
    config.cross_section   = :default
end

exMC_pphi = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic1_ppbarphi"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_pphip
    config.cross_section   = :default
end

### Event selection (BOSS) — psi' -> gamma chi_cJ -> gamma p pbar K+ K- ###
alg_name = "PsipGamChicJPPbarKK"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta  0.93   # |cos(theta)| < 0.93
                  Vz        10.0    # closest approach within +-10 cm of the IP along the beam
                  Vr         1.0    # within 1 cm of the beamline in the transverse plane
                  nChrp    "==2"
                  nChrn    "==2"
                  nNet     "==0"    # four tracks identified as p, pbar, K+, K-
                }
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0    # shower kept > 10 deg from the nearest charged track
                  energyThreshold_b 0.080   # EMC energy deposition > 80 MeV
                  energyThreshold_e 0.080
                  nGam              ">=1"   # at least one photon candidate (radiative photon)
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]   # p+ and anti-p-
                  identify :kaon,   against: [:pion, :proton] # K+ and K-
                  nprp "==1"
                  nprm "==1"
                  nkp  "==1"
                  nkm  "==1"
                }
               # 4C kinematic fit to psi' -> gamma p pbar K+ K-
               # (energy-momentum conservation). The photon combination with the
               # smallest chi2_4C is retained when more than one photon candidate exists.
               .kinematic_fit([:gamma, :prp, :prm, :kp, :km]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200   # loose BOSS cut; the analysis-level chi2_4C cut is applied in ROOT
                }

algorithm.note(:background_veto,
               "The p pbar K+ K- final state is generated with exactly four charged tracks " \
               "(p, pbar, K+, K-) plus the radiative photon; the Lambda(1520)/phi intermediate " \
               "states are not reconstructed in the BOSS selection, they are extracted from " \
               "the pK-, pbarK+ and K+K- invariant mass spectra at the ROOT level.")
         .note(:intermediate_vetoes,
               "For the direct chi_cJ -> p pbar K+ K- branching fraction, the Lambda(1520), " \
               "Lambda(1520)bar and phi intermediate states are removed with the mass windows " \
               "|M(pK-) - 1.52| > 0.07 GeV/c^2, |M(pbarK+) - 1.52| > 0.07 GeV/c^2 and " \
               "|M(K+K-) - 1.02| > 0.03 GeV/c^2. Applied in ROOT after the 4C fit.")
         .note(:chi_cJ_mass_windows,
               "chi_cJ invariant-mass windows used for the intermediate-state analyses: " \
               "chi_c0 : 3.365 - 3.455 GeV/c^2; chi_c1 : 3.490 - 3.530 GeV/c^2; " \
               "chi_c2 : 3.530 - 3.580 GeV/c^2. Applied in ROOT.")
         .note(:lambda1520_selection,
               "Lambda(1520)Lambda(1520) candidates require |M(pbarK+) - 1.520| < 0.05 GeV/c^2 " \
               "and |M(pK-) - 1.520| < 0.05 GeV/c^2 after rejecting phi p pbar with " \
               "|M(K+K-) - 1.02| > 0.03 GeV/c^2. Two-dimensional sideband regions A, B, C " \
               "estimate the background in the signal region S; all applied in ROOT.")
         .note(:background_studies,
               "Backgrounds checked with an inclusive psi' MC sample (100e6 events) and with " \
               "exclusive MC for psi' -> pi0 p pbar K+ K- (2e5 events), and " \
               "psi' -> gamma chi_cJ -> gamma K+K-K+K-, gamma K+K-pi+pi-, gamma p pbar pi+pi-, " \
               "psi' -> p pbar K+ K- (1e5 events each). After event selection 265 events survive, " \
               "all from psi' -> pi0 p pbar K+ K-; the estimated background contribution is " \
               "about 1.4 events. A 42.9 pb^-1 continuum data sample at 3.65 GeV is used to " \
               "investigate continuum backgrounds: no events survive the selection.")
         .note(:mass_spectrum_fit,
               "The p pbar K+ K- mass distribution is fitted with Breit-Wigner functions " \
               "convolved with Gaussian resolution functions for the chi_cJ signals plus a flat " \
               "background; the chi_cJ widths are fixed to the PDG values and the instrumental " \
               "resolution is about 4 MeV/c^2. The intermediate-state spectra are fitted with " \
               "Breit-Wigner * Gaussian signals on Chebyshev-polynomial backgrounds. " \
               "Fitted yields (statistical only): chi_cJ -> p pbar K+ K-: " \
               "48.2+-7.7 (chi_c0), 81.5+-9.2 (chi_c1), 131+-12 (chi_c2); " \
               "chi_cJ -> pbar K+ Lambda(1520)+c.c.: 62+-12, 48+-10, 79+-13; " \
               "chi_cJ -> Lambda(1520)Lambda(1520): 28.1+-9.8 (chi_c0), <6.9 (chi_c1, 90% C.L.), " \
               "28.9+-7.4 (chi_c2); chi_cJ -> p pbar phi: 42.4+-8.2 (chi_c0), <13.3 (chi_c1), " \
               "24.4+-6.8 (chi_c2). Performed at the ROOT level.")
         .note(:branching_fraction_normalization,
               "B = N_obs / [N_psi' * B(psi' -> gamma chi_cJ) * (intermediate Br) * epsilon], " \
               "with N_psi' = 106e6 (4% uncertainty) and the PDG values " \
               "B(psi' -> gamma chi_c0) = (9.62+-0.31)%, B(psi' -> gamma chi_c1) = (9.2+-0.4)%, " \
               "B(psi' -> gamma chi_c2) = (8.74+-0.35)%; " \
               "B(Lambda(1520) -> pK-) = 22.5% and B(phi -> K+K-) = 48.9% are used for the " \
               "resonant channels. No BOSS-side action.")
         .note(:efficiency_angular_model,
               "Detection efficiencies are determined from MC that assumes an angular " \
               "distribution 1 + alpha cos^2(theta) for the two-body decays; alpha is obtained " \
               "from fits to the data cos(theta) distribution separately for chi_c0, chi_c1 and " \
               "chi_c2.")
         .note(:systematics,
               "Tracking 2% per charged track (8% for the four-track final state), PID 2% per " \
               "particle (8% total), photon reconstruction 1% per photon, kinematic fit 1.4% / " \
               "1.6% / 2.3% for chi_c0 / chi_c1 / chi_c2, fitting procedure, mass-window " \
               "variations, alpha value (decay-model) variation from -1 to 1, intermediate " \
               "branching fractions from the PDG, and 4% on N_psi'.")

algorithm.with_decay_card(decay_card_chic1).apply(event_selection)
root_files = algorithm.execute_on([
  psip_data, psip_incMC,
  exMC_chic0, exMC_chic1, exMC_chic2,
  exMC_lambda1520, exMC_pphi
])
