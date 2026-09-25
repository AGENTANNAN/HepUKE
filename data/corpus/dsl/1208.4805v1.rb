# =============================================================================
# BESIII paper arXiv:1208.4805v1
# Search for the hadronic transition chi_cJ -> eta_c pi+ pi- and observation of
# chi_cJ -> K Kbar pi pi pi (J = 0, 1, 2), produced via psi(3686) -> gamma chi_cJ.
#
# Data: 1.06 x 10^8 psi(3686) events + 42 pb^-1 continuum data at sqrt(s) = 3.65 GeV.
#
# Two distinct final states are reconstructed, each covering both the direct
# K Kbar pi pi pi decay and the eta_c pi+ pi- transition with the corresponding
# eta_c decay mode:
#   Mode A: psi(3686) -> gamma chi_cJ, chi_cJ -> K_S0 K+- pi-+ pi+ pi-
#                                     chi_cJ -> eta_c pi+ pi-, eta_c -> K_S0 K+- pi-+
#           (>= 6 charged tracks, >= 1 photon)
#   Mode B: psi(3686) -> gamma chi_cJ, chi_cJ -> K+ K- pi+ pi- pi0
#                                     chi_cJ -> eta_c pi+ pi-, eta_c -> K+ K- pi0
#           (4 charged tracks, >= 3 photons)
# Within one mode the chi_c0, chi_c1 and chi_c2 decays share the same final state
# and the same selection; only the decay-card mother differs.
# =============================================================================

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) data, 1.06e8 events
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # psi(3686) inclusive MC
cont_data  = DatasetManager.real_data.find("709_3650")     # continuum data at 3.65 GeV
cont_incMC = DatasetManager.inclusive_mc.find("709_3650")  # continuum inclusive MC

# -----------------------------------------------------------------------------
# Decay cards (EvtGen format)
# -----------------------------------------------------------------------------
# --- Mode A: chi_cJ -> K_S0 K+- pi-+ pi+ pi- ---
decay_card_kskpipmpipm_c0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0                 PHSP;
    Enddecay

    Decay chi_c0
    1.0000 K_S0 K+ pi- pi+ pi-          PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                      PHSP;
    Enddecay

    End
DECAYCARD

decay_card_kskpipmpipm_c1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1                 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 K_S0 K+ pi- pi+ pi-          PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                      PHSP;
    Enddecay

    End
DECAYCARD

decay_card_kskpipmpipm_c2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2                 PHSP;
    Enddecay

    Decay chi_c2
    1.0000 K_S0 K+ pi- pi+ pi-          PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                      PHSP;
    Enddecay

    End
DECAYCARD

# --- Mode B: chi_cJ -> K+ K- pi+ pi- pi0 ---
decay_card_kkpipipi0_c0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0                 PHSP;
    Enddecay

    Decay chi_c0
    1.0000 K+ K- pi+ pi- pi0            PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                  PHSP;
    Enddecay

    End
DECAYCARD

decay_card_kkpipipi0_c1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1                 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 K+ K- pi+ pi- pi0            PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                  PHSP;
    Enddecay

    End
DECAYCARD

decay_card_kkpipipi0_c2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2                 PHSP;
    Enddecay

    Decay chi_c2
    1.0000 K+ K- pi+ pi- pi0            PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                  PHSP;
    Enddecay

    End
DECAYCARD

# --- Mode A (eta_c channel): chi_cJ -> eta_c pi+ pi-, eta_c -> K_S0 K+- pi-+ ---
decay_card_etac_kskpim_c0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0                 PHSP;
    Enddecay

    Decay chi_c0
    1.0000 eta_c pi+ pi-                PHSP;
    Enddecay

    Decay eta_c
    1.0000 K_S0 K+ pi-                  PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                      PHSP;
    Enddecay

    End
DECAYCARD

decay_card_etac_kskpim_c1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1                 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 eta_c pi+ pi-                PHSP;
    Enddecay

    Decay eta_c
    1.0000 K_S0 K+ pi-                  PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                      PHSP;
    Enddecay

    End
DECAYCARD

decay_card_etac_kskpim_c2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2                 PHSP;
    Enddecay

    Decay chi_c2
    1.0000 eta_c pi+ pi-                PHSP;
    Enddecay

    Decay eta_c
    1.0000 K_S0 K+ pi-                  PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                      PHSP;
    Enddecay

    End
DECAYCARD

# --- Mode B (eta_c channel): chi_cJ -> eta_c pi+ pi-, eta_c -> K+ K- pi0 ---
decay_card_etac_kkpi0_c0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0                 PHSP;
    Enddecay

    Decay chi_c0
    1.0000 eta_c pi+ pi-                PHSP;
    Enddecay

    Decay eta_c
    1.0000 K+ K- pi0                    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                  PHSP;
    Enddecay

    End
DECAYCARD

decay_card_etac_kkpi0_c1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1                 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 eta_c pi+ pi-                PHSP;
    Enddecay

    Decay eta_c
    1.0000 K+ K- pi0                    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                  PHSP;
    Enddecay

    End
DECAYCARD

decay_card_etac_kkpi0_c2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2                 PHSP;
    Enddecay

    Decay chi_c2
    1.0000 eta_c pi+ pi-                PHSP;
    Enddecay

    Decay eta_c
    1.0000 K+ K- pi0                    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                  PHSP;
    Enddecay

    End
DECAYCARD

# -----------------------------------------------------------------------------
# Exclusive MC samples (one per chi_cJ state and per final state)
# -----------------------------------------------------------------------------
exMC_kskpipmpipm_c0 = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_gamma_chic0_KSKpipimpipm"; c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_kskpipmpipm_c0; c.cross_section = :default }
exMC_kskpipmpipm_c1 = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_gamma_chic1_KSKpipimpipm"; c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_kskpipmpipm_c1; c.cross_section = :default }
exMC_kskpipmpipm_c2 = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_gamma_chic2_KSKpipimpipm"; c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_kskpipmpipm_c2; c.cross_section = :default }

exMC_kkpipipi0_c0 = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_gamma_chic0_KKpipipi0"; c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_kkpipipi0_c0; c.cross_section = :default }
exMC_kkpipipi0_c1 = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_gamma_chic1_KKpipipi0"; c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_kkpipipi0_c1; c.cross_section = :default }
exMC_kkpipipi0_c2 = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_gamma_chic2_KKpipipi0"; c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_kkpipipi0_c2; c.cross_section = :default }

exMC_etac_kskpim_c0 = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_gamma_chic0_etac_KSKpi"; c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_etac_kskpim_c0; c.cross_section = :default }
exMC_etac_kskpim_c1 = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_gamma_chic1_etac_KSKpi"; c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_etac_kskpim_c1; c.cross_section = :default }
exMC_etac_kskpim_c2 = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_gamma_chic2_etac_KSKpi"; c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_etac_kskpim_c2; c.cross_section = :default }

exMC_etac_kkpi0_c0 = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_gamma_chic0_etac_KKpi0"; c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_etac_kkpi0_c0; c.cross_section = :default }
exMC_etac_kkpi0_c1 = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_gamma_chic1_etac_KKpi0"; c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_etac_kkpi0_c1; c.cross_section = :default }
exMC_etac_kkpi0_c2 = DatasetManager.create_exclusive_mc { |c| c.sample_name = "psip_gamma_chic2_etac_KKpi0"; c.related_dataset = psip_data; c.events = 500_000; c.decay_card = decay_card_etac_kkpi0_c2; c.cross_section = :default }

# =============================================================================
# ALGORITHM A: psi(3686) -> gamma chi_cJ,
#              chi_cJ -> K_S0 K+- pi-+ pi+ pi-  (also eta_c -> K_S0 K+- pi-+)
# =============================================================================
alg_name_A = "PsipGammaChicJKSKpipimpipm"
alg_A = Algorithm.new(alg_name_A)
alg_A.set_header(["#{alg_name_A}Alg/#{alg_name_A}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})

sel_A = Selection.new
sel_A.select_track {
        cos_theta 0.93    # within the MDC angular coverage |cos(theta)| < 0.93
        Vz        10.0    # within 10 cm of the IP along the beam axis
        Vr         1.0    # within 1 cm of the IP transverse to the beam line
        nChrp     "==3"   # at least six charged tracks: K_S0 daughters + 4 more
        nChrn     "==3"
        nNet      "==0"   # zero net charge for the four additional tracks
      }
      .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14    # EMC timing requirement suppresses noise
        angle_to_track    20.0  # isolated showers >= 20 deg away from any charged track
        energyThreshold_b 0.025 # barrel (|cos(theta)| < 0.8)     : E > 25 MeV
        energyThreshold_e 0.050 # endcaps (0.86 < |cos(theta)| < 0.92) : E > 50 MeV
        nGam              ">=1" # at least one good photon
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :pion,   against: [:kaon, :proton]
        identify :kaon,   against: [:pion, :proton]
        identify :proton, against: [:kaon, :pion]
        nkp  ">=1"   # the single charged kaon (K+ or K- of the charge-conjugate pair)
        npip ">=2"
        npim ">=2"
      }
      # K_S0 candidates: secondary vertex fit to the charged-track pairs treated as
      # pions; the combination with the best fit quality is kept.
      .secondary_vertex_fit([:pip, :pim]) {
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
      # 4C kinematic fit: psi(3686) -> gamma K_S0 K+- pi-+ pi+ pi-, constrained by
      # the initial e+e- four-momentum. The reconstructed K_S0 information from the
      # secondary vertex fit is used as input.
      .kinematic_fit([:gamma, :K_S0, :kp, :pip, :pim, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200   # loose BOSS cut; the paper's chi2_4C < 50 is applied in ROOT
      }
      # Competing-hypothesis fit (no chi2_cut, no nominal): stores chi2_4C for the
      # gamma gamma K_S0 K+- pi-+ pi+ pi- hypothesis, used for the ROOT-level veto
      # against background with one extra photon.
      .kinematic_fit([:gamma, :gamma, :K_S0, :kp, :pip, :pim, :pim]) {
        constrain_four_momentum
      }

alg_A
  .note(:ks0_selection,
        "K_S0 candidates are reconstructed from secondary vertex fits to all charged-track " \
        "pairs in the event (assuming the tracks to be pions). The combination with the best " \
        "fit quality is kept, and the K_S0 candidate must have an invariant mass within " \
        "7 MeV/c^2 of the nominal K_S0 mass and the secondary vertex must be at least 0.5 cm " \
        "away from the IP. The reconstructed K_S0 information is used as input for the " \
        "subsequent kinematic fit.")
  .note(:charge_conjugate_channel,
        "the K_S0 K+ pi- pi+ pi- and K_S0 K- pi+ pi+ pi- final states are summed and quoted " \
        "together as K_S0 K+- pi-+ pi+ pi-; the selection and the fit are identical for the " \
        "two charge combinations.")
  .note(:particle_species_selection,
        "the species of the final-state particles, the best photon and the best pi0 candidate " \
        "are determined by selecting the combination with the minimum value of " \
        "chi^2 = chi^2_4C + sum_j chi^2_PID(j), where chi^2_PID(j) is built from the dE/dx and " \
        "TOF measurements; events with chi^2_4C < 50 are kept. Implementing the combined " \
        "chi^2_4C + sum chi^2_PID minimisation requires the PID chi^2 values to be stored " \
        "alongside the kinematic-fit chi^2 (ROOT level).")
  .note(:jpsi_veto,
        "background from psi(3686) -> X + J/psi is removed by requiring the recoil mass of any " \
        "pi+ pi- pair to lie outside a +-3 sigma window around the nominal J/psi mass, where " \
        "sigma is the resolution of the pi+ pi- recoil mass. Uses the 4C-fit-corrected " \
        "four-momenta; applied at the ROOT level.")
  .note(:eta_veto,
        "psi(3686) -> eta J/psi with eta -> gamma pi+ pi- is rejected by discarding events with " \
        "0.535 GeV/c^2 < M(gamma pi+ pi-) < 0.555 GeV/c^2. Applied at the ROOT level.")
  .note(:extra_photon_veto,
        "to reject background with one more photon than the signal, chi^2_4C > 10 is required " \
        "when a 4C kinematic fit to gamma gamma K_S0 K+- pi-+ pi+ pi- is applied to the event; " \
        "the competing-hypothesis chi^2 is stored by the additional non-nominal kinematic fit " \
        "and the comparison is applied at the ROOT level.")
  .note(:helicity_correction,
        "the helix parameters of the charged tracks in the MC simulation are corrected (smeared " \
        "with a Gaussian whose mean and width are taken from the data/MC pull distributions of " \
        "the control sample J/psi -> phi f_0(980), phi -> K+ K-, f_0(980) -> pi+ pi-) to reduce " \
        "the data-MC difference of the chi^2_4C distribution. The efficiency obtained from the " \
        "track-parameter-corrected MC is taken as the nominal value and half of the difference " \
        "between the corrected and uncorrected MC samples is assigned as the systematic " \
        "uncertainty from the kinematic fitting. The correction is applied to the MC samples " \
        "outside the DSL grammar.")
  .note(:ks0_efficiency_correction,
        "the K_S0 reconstruction efficiency is corrected as a function of the K_S0 momentum " \
        "using the data/MC difference measured with J/psi -> K*+- K-+, K*+- -> K0 pi+-; the " \
        "uncertainty of this correction is taken as a systematic error.")
  .note(:fitting_range_and_background,
        "the M(K_S0 K+- pi-+ pi+ pi-) spectrum (3.3 - 3.6 GeV/c^2) is fitted with the MC signal " \
        "shapes of the three chi_cJ states convolved with a Gaussian to account for the " \
        "data-MC difference in mass scale and resolution, and a second-order Chebyshev function " \
        "for the background. In the chi_c0 MC generation the E1 radiative transition factor " \
        "E_gamma^3 is included together with the KEDR damping function. Peaking backgrounds " \
        "(chi_c0 -> K_S0 K_S0 pi+ pi-, chi_c0 -> K_S0 K_S0 K+ K-, etc.) are estimated with the " \
        "inclusive MC and subtracted. ROOT level.")
  .note(:narrow_resonance_rejection,
        "for the measurement of B(chi_cJ -> K+ K- pi+ pi- pi0) the contributions from the " \
        "narrow intermediate resonances eta, omega and phi are excluded by requiring " \
        "|M(pi+ pi- pi0) - m_eta| > 15 MeV/c^2, |M(pi+ pi- pi0) - m_omega| > 40 MeV/c^2, " \
        "|M(pi+ pi- pi0) - m_phi| > 15 MeV/c^2 and |M(K+ K-) - m_phi| > 15 MeV/c^2. These " \
        "requirements also remove chi_cJ -> phi phi, phi -> K+ K- and phi -> pi+ pi- pi0. " \
        "Applied at the ROOT level.")
  .note(:continuum_background,
        "background from the continuum process e+e- -> q qbar is estimated with the 42 pb^-1 " \
        "data sample taken at sqrt(s) = 3.65 GeV; it is small and uniformly distributed over " \
        "the mass region of interest and is absorbed in the smooth background term of the fit.")
  .note(:systematic_uncertainties,
        "MDC tracking 8.0% (four tracks), photon reconstruction 1.0%, MC statistics 1.1-1.5%, " \
        "K_S0 reconstruction 1.4-1.7%, kinematic fit 1.5-1.9%, damping function 0.1-0.5%, " \
        "intermediate states 1.0%, fitting range 0.2-1.0%, background shape 0.6-1.4%, " \
        "B(psi(3686) -> gamma chi_cJ) 3.2-4.4%, number of psi(3686) events 4.0%; total " \
        "10.1-10.5%.")

alg_A.with_decay_card(decay_card_kskpipmpipm_c0).apply(sel_A)
alg_A.execute_on([psip_data, psip_incMC, cont_data, cont_incMC,
                  exMC_kskpipmpipm_c0, exMC_kskpipmpipm_c1, exMC_kskpipmpipm_c2,
                  exMC_etac_kskpim_c0, exMC_etac_kskpim_c1, exMC_etac_kskpim_c2])

# =============================================================================
# ALGORITHM B: psi(3686) -> gamma chi_cJ,
#              chi_cJ -> K+ K- pi+ pi- pi0  (also eta_c -> K+ K- pi0)
# =============================================================================
alg_name_B = "PsipGammaChicJKKpipipi0"
alg_B = Algorithm.new(alg_name_B)
alg_B.set_header(["#{alg_name_B}Alg/#{alg_name_B}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})

sel_B = Selection.new
sel_B.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr         1.0
        nChrp     "==2"   # four good charged tracks with zero net charge
        nChrn     "==2"
        nNet      "==0"
      }
      .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    20.0   # isolated showers >= 20 deg away from any charged track
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=3"  # at least three good photons
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :pion,   against: [:kaon, :proton]
        identify :kaon,   against: [:pion, :proton]
        identify :proton, against: [:kaon, :pion]
        nkp  "==1"
        nkm  "==1"
        npip "==1"
        npim "==1"
      }
      # pi0 candidates: photon pairs with 0.120 GeV/c^2 < M(gamma gamma) < 0.145 GeV/c^2.
      # The pair with invariant mass closest to the nominal pi0 mass is taken as the
      # pi0 candidate and its mass is constrained by the Kalman fit.
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      }
      # 4C kinematic fit: psi(3686) -> gamma K+ K- pi+ pi- pi0, constrained by the
      # initial e+e- four-momentum.
      .kinematic_fit([:gamma, :kp, :km, :pip, :pim, :pi0]) {
        nominal
        constrain_four_momentum
        chi2_cut 200   # loose BOSS cut; the paper's chi2_4C < 50 is applied in ROOT
      }
      # Competing-hypothesis fit (no chi2_cut, no nominal): stores chi2_4C for the
      # gamma gamma K+ K- pi+ pi- pi0 hypothesis, used for the ROOT-level veto
      # against background with one extra photon.
      .kinematic_fit([:gamma, :gamma, :kp, :km, :pip, :pim, :pi0]) {
        constrain_four_momentum
      }

alg_B
  .note(:pi0_selection,
        "pi0 candidates are reconstructed from pairs of photons with an invariant mass in the " \
        "range 0.120 GeV/c^2 < M(gamma gamma) < 0.145 GeV/c^2.")
  .note(:particle_species_selection,
        "the species of the final-state particles, the best photon and the best pi0 candidate " \
        "among the extra photons are determined by selecting the combination with the minimum " \
        "value of chi^2 = chi^2_4C + sum_j chi^2_PID(j); events with chi^2_4C < 50 are kept. " \
        "The combined minimisation is completed at the ROOT level using the stored chi^2 values.")
  .note(:jpsi_veto,
        "background from psi(3686) -> X + J/psi is removed by requiring the recoil mass of any " \
        "pi+ pi- pair to lie outside a +-3 sigma window around the nominal J/psi mass. Applied " \
        "at the ROOT level.")
  .note(:eta_veto,
        "psi(3686) -> eta J/psi with eta -> gamma pi+ pi- is rejected by discarding events with " \
        "0.535 GeV/c^2 < M(gamma pi+ pi-) < 0.555 GeV/c^2. Applied at the ROOT level.")
  .note(:extra_photon_veto,
        "to reject background with one more photon than the signal, chi^2_4C > 10 is required " \
        "when a 4C kinematic fit to gamma gamma K+ K- pi+ pi- pi0 is applied to the event; the " \
        "competing-hypothesis chi^2 is stored by the additional non-nominal kinematic fit and " \
        "the comparison is applied at the ROOT level.")
  .note(:psi_to_pi0pi0jpsi_veto,
        "to suppress psi(3686) -> pi0 pi0 J/psi with J/psi -> K+ K- pi+ pi-, " \
        "M(K+ K- pi+ pi-) > 3.125 GeV/c^2 or M(K+ K- pi+ pi-) < 3.095 GeV/c^2 is required. " \
        "Applied at the ROOT level.")
  .note(:helicity_correction,
        "the helix parameters of the charged tracks in the MC simulation are corrected using the " \
        "control sample J/psi -> phi f_0(980), phi -> K+ K-, f_0(980) -> pi+ pi-; the efficiency " \
        "from the corrected MC is the nominal value and half of the corrected/uncorrected " \
        "difference is the systematic uncertainty of the kinematic fit. Applied outside the " \
        "DSL grammar.")
  .note(:narrow_resonance_rejection,
        "for the measurement of B(chi_cJ -> K+ K- pi+ pi- pi0) the contributions from the " \
        "narrow intermediate resonances eta, omega and phi are excluded by requiring " \
        "|M(pi+ pi- pi0) - m_eta| > 15 MeV/c^2, |M(pi+ pi- pi0) - m_omega| > 40 MeV/c^2, " \
        "|M(pi+ pi- pi0) - m_phi| > 15 MeV/c^2 and |M(K+ K-) - m_phi| > 15 MeV/c^2; this also " \
        "removes chi_cJ -> phi phi, phi -> K+ K- and phi -> pi+ pi- pi0. The branching fractions " \
        "of chi_cJ -> eta K+ K-, omega K+ K-, phi K+ K- and phi pi+ pi- pi0 are extracted from " \
        "fits to the M(pi+ pi- pi0) and M(K+ K-) spectra in the three chi_cJ signal regions. " \
        "Applied at the ROOT level.")
  .note(:fitting_range_and_background,
        "the M(K+ K- pi+ pi- pi0) spectrum is fitted with the MC signal shapes of the three " \
        "chi_cJ states convolved with a Gaussian plus a second-order Chebyshev background. " \
        "Peaking backgrounds from chi_cJ -> K_S0 K+- pi-+ pi0 are estimated with the inclusive " \
        "MC and subtracted. ROOT level.")
  .note(:continuum_background,
        "background from the continuum process e+e- -> q qbar is estimated with the 42 pb^-1 " \
        "data sample at sqrt(s) = 3.65 GeV; it is small and uniformly distributed and is " \
        "absorbed in the smooth background term of the fit.")
  .note(:systematic_uncertainties,
        "MDC tracking 8.0% (four tracks), photon reconstruction 3.0% (three photons), MC " \
        "statistics 1.0-1.4%, pi0 reconstruction 1.0%, kinematic fit 0.2-0.4%, damping function " \
        "0.1-0.4%, intermediate states 4.0%, fitting range 0.4-0.7%, background shape 0.4-1.3%, " \
        "B(psi(3686) -> gamma chi_cJ) 3.2-4.4%, number of psi(3686) events 4.0%; total " \
        "10.9-11.3%.")

alg_B.with_decay_card(decay_card_kkpipipi0_c0).apply(sel_B)
alg_B.execute_on([psip_data, psip_incMC, cont_data, cont_incMC,
                  exMC_kkpipipi0_c0, exMC_kkpipipi0_c1, exMC_kkpipipi0_c2,
                  exMC_etac_kkpi0_c0, exMC_etac_kkpi0_c1, exMC_etac_kkpi0_c2])
