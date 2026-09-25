# BESIII paper arXiv:1205.5103v2
# First observation of the M1 transition psi(3686) -> gamma eta_c(2S) with
# eta_c(2S) -> K_S0 K+- pi-+ and eta_c(2S) -> K+ K- pi0.
# Data: 106 M psi(3686) events (156 pb^-1) at sqrt(s) = 3.686 GeV,
# plus 42 pb^-1 at sqrt(s) = 3.65 GeV for continuum background studies.
# The two final states have different track/photon multiplicities and different
# kinematic-fit hypotheses, so each charge-conjugate mode gets its own Algorithm.
# Nominal fits: 4C (gamma K_S0 K pi) and 5C (gamma K+ K- pi0).

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) data, 106 M events (156 pb^-1)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # psi(3686) inclusive MC (lundcharm, for background studies)
off_data   = DatasetManager.real_data.find("709_3650")     # continuum data at sqrt(s) = 3.65 GeV (42 pb^-1)
off_incMC  = DatasetManager.inclusive_mc.find("709_3650")

### Decay cards (EvtGen format) ###
# Signal: psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K_S0 K+ pi-
decay_card_ks_kp_pim = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  eta_c(2S)         PHSP;
    Enddecay

    Decay eta_c(2S)
    1.0000  K_S0  K+  pi-           PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                PHSP;
    Enddecay

    End
DECAYCARD

# Signal: psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K_S0 K- pi+ (charge conjugate)
decay_card_ks_km_pip = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  eta_c(2S)         PHSP;
    Enddecay

    Decay eta_c(2S)
    1.0000  K_S0  K-  pi+           PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                PHSP;
    Enddecay

    End
DECAYCARD

# Signal: psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K+ K- pi0, pi0 -> gamma gamma
decay_card_kk_pi0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  eta_c(2S)         PHSP;
    Enddecay

    Decay eta_c(2S)
    1.0000  K+  K-  pi0             PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma            PHSP;
    Enddecay

    End
DECAYCARD

# Background (1): psi(3686) -> K_S0 K+- pi-+ with a fake photon candidate
decay_card_bkg_ks_kpi = <<~DECAYCARD
    Decay psi(2S)
    1.0000  K_S0  K+  pi-           PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                PHSP;
    Enddecay

    End
DECAYCARD

# Background (1): psi(3686) -> K+ K- pi0 with a fake photon candidate
decay_card_bkg_kk_pi0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000  K+  K-  pi0             PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma            PHSP;
    Enddecay

    End
DECAYCARD

# Background (2): same final states with an extra ISR/FSR photon (FSR via PHOTOS),
# and psi(3686) -> omega K+ K- with omega -> gamma pi0
decay_card_bkg_ks_kpi_fsr = <<~DECAYCARD
    Decay psi(2S)
    1.0000  K_S0  K+  pi-           PHOTOS PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                PHSP;
    Enddecay

    End
DECAYCARD

decay_card_bkg_omega_kk = <<~DECAYCARD
    Decay psi(2S)
    1.0000  omega  K+  K-           PHSP;
    Enddecay

    Decay omega
    1.0000  gamma  pi0              PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma            PHSP;
    Enddecay

    End
DECAYCARD

# Background (3): an extra photon, primarily psi(3686) -> pi0 K_S0 K+- pi-+ and
# psi(3686) -> pi0 K+ K- pi0 with pi0 -> gamma gamma
decay_card_bkg_pi0_ks_kpi = <<~DECAYCARD
    Decay psi(2S)
    1.0000  pi0  K_S0  K+  pi-      PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma            PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                PHSP;
    Enddecay

    End
DECAYCARD

decay_card_bkg_pi0_kk_pi0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000  pi0  K+  K-  pi0        PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma            PHSP;
    Enddecay

    End
DECAYCARD

# Background: psi(3686) -> pi+ pi- J/psi with J/psi -> K+ K- (and leptonic J/psi decays)
decay_card_bkg_pipi_jpsi_kk = <<~DECAYCARD
    Decay psi(2S)
    1.0000  pi+  pi-  J/psi         PHSP;
    Enddecay

    Decay J/psi
    1.0000  K+  K-                  PHSP;
    Enddecay

    End
DECAYCARD

decay_card_bkg_pipi_jpsi_ll = <<~DECAYCARD
    Decay psi(2S)
    1.0000  pi+  pi-  J/psi         PHSP;
    Enddecay

    Decay J/psi
    1.0000  e+  e-                  PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Background: psi(3686) -> eta J/psi with eta -> pi+ pi- pi0 and eta -> gamma gamma
decay_card_bkg_eta_jpsi_pipipi0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000  eta  J/psi              PHSP;
    Enddecay

    Decay eta
    1.0000  pi+  pi-  pi0           PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma            PHSP;
    Enddecay

    End
DECAYCARD

decay_card_bkg_eta_jpsi_gg = <<~DECAYCARD
    Decay psi(2S)
    1.0000  eta  J/psi              PHSP;
    Enddecay

    Decay eta
    1.0000  gamma  gamma            PHSP;
    Enddecay

    End
DECAYCARD

# Background: continuum e+e- -> gamma* -> K_S0 K+- pi-+ (estimated from the 3.65 GeV data)
decay_card_bkg_continuum_ks_kpi = <<~DECAYCARD
    Decay psi(4260)
    1.0000  K_S0  K+  pi-           PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples ###
exMC_ks_kp_pim = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_etac2S_ks_kp_pim"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_ks_kp_pim
  config.cross_section   = :default
end
exMC_ks_kp_pim.save_to_config(format: :yaml, file_path: 'exMC_ks_kp_pim')

exMC_ks_km_pip = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_etac2S_ks_km_pip"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_ks_km_pip
  config.cross_section   = :default
end
exMC_ks_km_pip.save_to_config(format: :yaml, file_path: 'exMC_ks_km_pip')

exMC_kk_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_etac2S_kk_pi0"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_kk_pi0
  config.cross_section   = :default
end
exMC_kk_pi0.save_to_config(format: :yaml, file_path: 'exMC_kk_pi0')

exMC_bkg_ks_kpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_ks_kpi_fake_gamma_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_ks_kpi
  config.cross_section   = :default
end

exMC_bkg_kk_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_kk_pi0_fake_gamma_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_kk_pi0
  config.cross_section   = :default
end

exMC_bkg_ks_kpi_fsr = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_ks_kpi_fsr_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_ks_kpi_fsr
  config.cross_section   = :default
end

exMC_bkg_omega_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_omega_kk_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_omega_kk
  config.cross_section   = :default
end

exMC_bkg_pi0_ks_kpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_pi0_ks_kpi_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_pi0_ks_kpi
  config.cross_section   = :default
end

exMC_bkg_pi0_kk_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_pi0_kk_pi0_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_pi0_kk_pi0
  config.cross_section   = :default
end

exMC_bkg_pipi_jpsi_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_pipi_jpsi_jpsi_to_kk_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_pipi_jpsi_kk
  config.cross_section   = :default
end

exMC_bkg_pipi_jpsi_ll = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_pipi_jpsi_jpsi_to_ll_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_pipi_jpsi_ll
  config.cross_section   = :default
end

exMC_bkg_eta_jpsi_pipipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_eta_jpsi_eta_to_3pi_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_eta_jpsi_pipipi0
  config.cross_section   = :default
end

exMC_bkg_eta_jpsi_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_eta_jpsi_eta_to_gg_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_eta_jpsi_gg
  config.cross_section   = :default
end

exMC_bkg_continuum_ks_kpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "continuum_gamma_to_ks_kpi_bkg"
  config.related_dataset = off_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_continuum_ks_kpi
  config.cross_section   = :default
end

### Event selection (BOSS) — channel 1: psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K_S0 K+ pi- ###
alg1_name = "PsipToGammaEtac2S_KsKPi"
alg1 = Algorithm.new(alg1_name)
alg1.set_header(["#{alg1_name}Alg/#{alg1_name}.h"])
    .set_constant({"ECMS" => [:double, 3.686]})
    .set_alias({"std::vector<double>" => "Vdouble"})

sel1 = Selection.new
sel1.select_track {
       cos_theta 0.93     # |cos theta| < 0.93 for each charged track
       Vz        10.0     # track origin within 10 cm of the IP along the beam axis
       Vr        1.0      # track origin within 1 cm of the IP transverse to the beam line
       nChrp     "==2"    # exactly four charged tracks with zero net charge
       nChrn     "==2"
       nNet      "==0"
     }
     .select_photon {
       tdc_emc_start     0      # EMC timing requirement (noise / out-of-time suppression)
       tdc_emc_end      14
       energyThreshold_b 0.025  # E_min = 25 MeV in the barrel
       energyThreshold_e 0.050  # E_min = 25 MeV in the endcap
       angle_to_track   10.0
       nGam             ">=1"   # at least one good photon (the M1 transition photon)
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :kaon, against: [:pion, :proton]  # dE/dx + TOF chi2_PID(K) for the K+ track
       nkp ">=1"
       nkm ">=1"
     }
     .secondary_vertex_fit([:pip, :pim]) {          # K_S0 -> pi+ pi- from all opposite-charge pairs
       build_virtual_particle(:K_S0).by_minimizing_verfit_chi2  # combination with best vertex-fit quality
       remove_used_particle_from_candidate_list
     }
     .kinematic_fit([:K_S0, :kp, :pim, :gamma]) {
       nominal
       constrain_four_momentum   # 4C fit (four-momentum conservation)
       chi2_cut 200              # loose BOSS cut; the paper requires chi2_4C < 50 (applied in ROOT)
     }

alg1.note(:ks_mass_window,
          "K_S0 candidates are required to have |M(pi+ pi-) - m(K_S0)| < 7 MeV/c^2 and a secondary " \
          "vertex at least 0.5 cm from the interaction point. Only the combination with the best " \
          "vertex-fit quality is retained. The fitted K_S0 information is used as input to the " \
          "subsequent 4C kinematic fit of the complete event.")
    .note(:background_veto,
          "To suppress gamma K_S0 K_S0 events, the tracks remaining after the K_S0 candidate is " \
          "formed are required NOT to form a good K_S0 candidate.")
    .note(:recoil_mass_veto,
          "Background from psi(3686) -> pi+ pi- J/psi (J/psi -> K+ K- and leptonic J/psi decays) and " \
          "psi(3686) -> eta J/psi is suppressed by requiring the recoil mass of all pi+ pi- pairs to " \
          "be less than 3.05 GeV/c^2.")
    .note(:best_candidate_selection,
          "The discrimination of the two charge-conjugate channels (K_S0 K+ pi- versus K_S0 K- pi+, " \
          "implemented as two separate Algorithms) and the selection of the best photon among " \
          "multiple candidates are achieved by minimizing chi^2 = chi2_4C + chi2_PID(K) + chi2_PID(pi).")
    .note(:signal_mc_generation,
          "The signal is generated with the expected angular distribution for the M1 transition " \
          "psi(3686) -> gamma eta_c(2S); the subsequent eta_c(2S) -> K_S0 K+- pi-+ and " \
          "eta_c(2S) -> K+ K- pi0 decays are generated according to phase space. Inclusive MC " \
          "events for background studies are generated with lundcharm.")
    .note(:fit_details,
          "The eta_c(2S) mass spectrum is fitted over 3.46-3.71 GeV/c^2 simultaneously for the two " \
          "channels with four components: the eta_c(2S) signal, the chi_c1 and chi_c2 resonances, and " \
          "the summed background. The eta_c(2S) line shape is taken as " \
          "(E_gamma^3 x BW(m) x f_d(E_gamma) x epsilon(m)) convolved with a Gaussian resolution, with " \
          "the KEDR damping function f_d. Yields: 81 +/- 14 (K_S0 K+- pi-+) and 46 +/- 11 " \
          "(K+ K- pi0); efficiencies 25.6% and 20.2%; combined significance 11.1 sigma.")
    .note(:systematic_uncertainties,
          "Systematic uncertainties on the product branching fraction (relative): background shape " \
          "9.9%, damping function 19.6%, fitting range 1.3%, mass shift 0.4%, tracking 4.0%, photon " \
          "reconstruction 1.3%, particle identification 1.3%, K_S0 reconstruction 2.3%, kinematic " \
          "fitting 3.9%, eta_c(2S) decay dynamics 1.5%, number of psi(3686) events 4.0%; total 23.3%. " \
          "Mass and width systematic uncertainties: 1.6 MeV/c^2 and 4.8 MeV.")

alg1.with_decay_card(decay_card_ks_kp_pim).apply(sel1)

### Event selection (BOSS) — channel 2 (charge conjugate): eta_c(2S) -> K_S0 K- pi+ ###
alg2_name = "PsipToGammaEtac2S_KsKmPi"
alg2 = Algorithm.new(alg2_name)
alg2.set_header(["#{alg2_name}Alg/#{alg2_name}.h"])
    .set_constant({"ECMS" => [:double, 3.686]})
    .set_alias({"std::vector<double>" => "Vdouble"})

sel2 = Selection.new
sel2.select_track {
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==2"
       nChrn     "==2"
       nNet      "==0"
     }
     .select_photon {
       tdc_emc_start     0
       tdc_emc_end      14
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       angle_to_track   10.0
       nGam             ">=1"
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :kaon, against: [:pion, :proton]
       nkp ">=1"
       nkm ">=1"
     }
     .secondary_vertex_fit([:pip, :pim]) {
       build_virtual_particle(:K_S0).by_minimizing_verfit_chi2
       remove_used_particle_from_candidate_list
     }
     .kinematic_fit([:K_S0, :km, :pip, :gamma]) {
       nominal
       constrain_four_momentum   # 4C fit
       chi2_cut 200              # loose BOSS cut; the paper requires chi2_4C < 50 (applied in ROOT)
     }

alg2.note(:ks_mass_window,
          "K_S0 candidates are required to have |M(pi+ pi-) - m(K_S0)| < 7 MeV/c^2 and a secondary " \
          "vertex at least 0.5 cm from the interaction point; the best vertex-fit-quality " \
          "combination is retained and its fitted information feeds the 4C kinematic fit.")
    .note(:background_veto,
          "The remaining tracks are required not to form a good K_S0 candidate, suppressing " \
          "gamma K_S0 K_S0 events.")
    .note(:recoil_mass_veto,
          "The recoil mass of all pi+ pi- pairs is required to be less than 3.05 GeV/c^2 to suppress " \
          "psi(3686) -> pi+ pi- J/psi and psi(3686) -> eta J/psi contamination.")
    .note(:best_candidate_selection,
          "Charge-conjugate discrimination (K_S0 K- pi+ versus K_S0 K+ pi-) and best-photon selection " \
          "use the minimum of chi^2 = chi2_4C + chi2_PID(K) + chi2_PID(pi).")

alg2.with_decay_card(decay_card_ks_km_pip).apply(sel2)

### Event selection (BOSS) — channel 3: eta_c(2S) -> K+ K- pi0 (5C fit) ###
alg3_name = "PsipToGammaEtac2S_KKPi0"
alg3 = Algorithm.new(alg3_name)
alg3.set_header(["#{alg3_name}Alg/#{alg3_name}.h"])
    .set_constant({"ECMS" => [:double, 3.686]})
    .set_alias({"std::vector<double>" => "Vdouble"})

sel3 = Selection.new
sel3.select_track {
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==1"    # exactly two charged tracks with zero net charge
       nChrn     "==1"
       nNet      "==0"
     }
     .select_photon {
       tdc_emc_start     0
       tdc_emc_end      14
       energyThreshold_b 0.025  # E_min = 25 MeV in the barrel
       energyThreshold_e 0.050  # E_min = 25 MeV in the endcap
       angle_to_track   10.0
       nGam             ">=3"   # at least three good photons: the transition photon plus the pi0 daughters
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :kaon, against: [:pion, :proton]  # Prob_PID(K) > 0.001 and larger than any other hypothesis
       nkp "==1"
       nkm "==1"
     }
     .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma]) {
       nominal
       constrain_four_momentum                                       # 5C fit
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # pi0 mass as 5th constraint
       chi2_cut 200              # loose BOSS cut; the paper requires chi2_5C < 30 (applied in ROOT)
     }

alg3.note(:kaon_pid_requirement,
          "Both charged tracks are required to satisfy the kaon-hypothesis probability " \
          "Prob_PID(K) > 0.001 and Prob_PID(K) to be larger than the probability of any other " \
          "hypothesis (pion or proton). PID combines dE/dx and TOF information into chi2_PID(i).")
    .note(:best_candidate_selection,
          "The five-constraint (5C) kinematic fit, with the pi0 mass as the additional constraint, " \
          "is used to select the best transition photon and the best pi0 -> gamma gamma combination; " \
          "the combination with the smallest chi2_5C is retained.")
    .note(:dimuon_mass_veto,
          "Contamination from psi(3686) -> pi+ pi- J/psi with J/psi -> K+ K- (and leptonic J/psi " \
          "decays) and from psi(3686) -> eta J/psi is removed by requiring the invariant mass of the " \
          "two charged tracks, calculated under the muon mass hypothesis, to be less than " \
          "2.9 GeV/c^2.")
    .note(:photon_energy_float,
          "Because a fake photon adds no information to the fit, the eta_c(2S) mass is determined " \
          "from a modified kinematic fit in which the magnitude of the photon momentum is allowed to " \
          "float freely (3C for gamma K_S0 K+- pi-+ and 4C for gamma K+ K- pi0); a fake photon then " \
          "gives a photon momentum tending to zero. This modified fit is a separate, post-selection " \
          "fit and is not expressible in the DSL.")
    .note(:fit_details,
          "The K+ K- pi0 mass spectrum is fitted simultaneously with the K_S0 K+- pi-+ spectrum over " \
          "3.46-3.71 GeV/c^2, including the eta_c(2S) signal, the chi_c1 and chi_c2 resonances and " \
          "the summed background. The chi_c1/chi_c2 line shapes come from MC convolved with a " \
          "Gaussian; their widths are fixed to the PDG values, and the mass shift and resolution are " \
          "extrapolated linearly to the eta_c(2S) mass.")

alg3.with_decay_card(decay_card_kk_pi0).apply(sel3)

### Job generation ###
datasets = [psip_data, psip_incMC, off_data, off_incMC,
            exMC_ks_kp_pim, exMC_ks_km_pip, exMC_kk_pi0,
            exMC_bkg_ks_kpi, exMC_bkg_kk_pi0, exMC_bkg_ks_kpi_fsr, exMC_bkg_omega_kk,
            exMC_bkg_pi0_ks_kpi, exMC_bkg_pi0_kk_pi0, exMC_bkg_pipi_jpsi_kk,
            exMC_bkg_pipi_jpsi_ll, exMC_bkg_eta_jpsi_pipipi0, exMC_bkg_eta_jpsi_gg,
            exMC_bkg_continuum_ks_kpi]

root_files_1 = alg1.execute_on(datasets)
root_files_2 = alg2.execute_on(datasets)
root_files_3 = alg3.execute_on(datasets)
