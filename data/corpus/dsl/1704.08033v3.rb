# =============================================================================
# BESIII : Observation of e+e- -> eta h_c at center-of-mass energies from
#          4.085 to 4.600 GeV
# arXiv:1704.08033v3
#
# The E1 transition h_c -> gamma eta_c is reconstructed with the eta_c decaying
# to 16 exclusive hadronic final states X_i.  The signal is extracted from the
# eta recoil mass spectrum M_recoil(eta).  The Born cross sections are measured
# at 15 c.m. energy points between 4.085 and 4.600 GeV.
# BOSS part only : dataset preparation + event selection up to the 4C fit.
# =============================================================================

# =============================================================================
# Datasets : the 15 c.m. energy points of Table I (BOSS 703, XYZ scan samples)
# =============================================================================
scan_points = [
  DatasetManager.real_data.find("703_4090"),   # 4085.4 MeV,  52.4  pb^-1
  DatasetManager.real_data.find("703_4190"),   # 4188.6 MeV,  43.1  pb^-1
  DatasetManager.real_data.find("703_4210"),   # 4207.7 MeV,  54.6  pb^-1
  DatasetManager.real_data.find("703_4220"),   # 4217.1 MeV,  54.1  pb^-1
  DatasetManager.real_data.find("703_4230"),   # 4226.3 MeV, 1091.7 pb^-1
  DatasetManager.real_data.find("703_4245"),   # 4241.7 MeV,  55.6  pb^-1
  DatasetManager.real_data.find("703_4260"),   # 4258.0 MeV, 825.7  pb^-1
  DatasetManager.real_data.find("703_4310"),   # 4307.9 MeV,  44.9  pb^-1
  DatasetManager.real_data.find("703_4360"),   # 4358.3 MeV, 539.8  pb^-1
  DatasetManager.real_data.find("703_4390"),   # 4387.4 MeV,  55.2  pb^-1
  DatasetManager.real_data.find("703_4420"),   # 4415.6 MeV, 1073.6 pb^-1
  DatasetManager.real_data.find("703_4470"),   # 4467.1 MeV, 109.9  pb^-1
  DatasetManager.real_data.find("703_4530"),   # 4527.1 MeV, 110.0  pb^-1
  DatasetManager.real_data.find("703_4575"),   # 4574.5 MeV,  47.7  pb^-1
  DatasetManager.real_data.find("703_4600")    # 4599.5 MeV, 566.9  pb^-1
]

# Inclusive MC is produced with luminosity equivalent to data at sqrt(s) = 4.23,
# 4.26 and 4.36 GeV.
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4230"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4360")
]

# =============================================================================
# Mode 01 : eta_c -> p pbar
# =============================================================================
decay_card_m01 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 p+ anti-p-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m01 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_ppbar"
  config.events        = 40_000
  config.decay_card    = decay_card_m01
  config.cross_section = :default
end

alg_m01 = Algorithm.new("EtaHcEtacPPbar")
alg_m01.set_header(["EtaHcEtacPPbarAlg/EtaHcEtacPPbar.h"])
sel_m01 = Selection.new
sel_m01.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=1"    # p
             nChrn     ">=1"    # pbar
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=3"   # E1 gamma + the two photons of the primary eta
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :proton, against: [:kaon, :pion]
             nprp ">=1"; nprm ">=1"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25            # |M(gamma gamma) - m(eta)| < 15 MeV/c^2
             neta ">=1"
           }
          .kinematic_fit([:gamma, :eta, :prp, :prm]) {
             nominal
             vertex_fit([2, 3])     # common production vertex for the two charged tracks
             constrain_four_momentum
             chi2_cut 200
           }
alg_m01.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25 to suppress background events with different particle assignments. The loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "If several candidates are found, the one minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex is retained, where chi2_1C is the mass constraint of the two photons of each pi0/eta, chi2_PID is the sum over charged tracks of the PID-hypothesis chi2 and chi2_vertex is the chi2 of the K_S0 secondary-vertex fit. BOSS resolves the combination by the minimum chi2_4C of the kinematic fit.")
      .note(:vertex_fit, "A vertex fit constrains the run-by-run production vertex and all charged tracks to a common vertex; its chi2 enters the combined selection variable above.")
      .note(:e1_photon_energy, "The E1 photon energy is required to be in (400, 600) MeV, enlarged to (350, 650) MeV for the data sets collected at sqrt(s) > 4.416 GeV. This energy window is applied in the ROOT analysis.")
      .note(:eta_recoil_mass, "The eta candidate must have a recoil mass in (3480, 3600) MeV/c^2 (h_c signal region); the yield is extracted from a fit to the M_recoil(eta) distribution in ROOT.")
      .note(:eta_c_mass_window, "The invariant mass of the hadronic system is required to be within (2940, 3020) MeV/c^2; if more than one eta candidate is found in the h_c signal region, the one giving an eta_c mass closest to the nominal eta_c mass is selected.")
      .note(:pid_selection, "TOF and dE/dx information are combined into PID confidence levels for the pion, kaon and proton hypotheses; both PID and kinematic-fit information are used to determine the particle type of each charged track.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i); the radiative correction (1+delta) and the vacuum-polarisation factor |1+Pi|^2 are obtained iteratively from dedicated MC samples.")
      .note(:fit_model, "The 16 eta recoil mass distributions are fitted simultaneously with an unbinned maximum-likelihood method (signal shape from MC, linear background); performed in ROOT.")
alg_m01.with_decay_card(decay_card_m01).apply(sel_m01)
alg_m01.execute_on(scan_points + incMC_points + exMC_m01)

# =============================================================================
# Mode 02 : eta_c -> 2 (pi+ pi-)
# =============================================================================
decay_card_m02 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m02 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_2pipi"
  config.events        = 40_000
  config.decay_card    = decay_card_m02
  config.cross_section = :default
end

alg_m02 = Algorithm.new("EtaHcEtac2PiPi")
alg_m02.set_header(["EtaHcEtac2PiPiAlg/EtaHcEtac2PiPi.h"])
sel_m02 = Selection.new
sel_m02.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=2"
             nChrn     ">=2"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=3"   # E1 gamma + the two photons of the primary eta
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :pion, against: [:kaon, :proton]
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=1"
           }
          .kinematic_fit([:gamma, :eta, :pip, :pim, :pip, :pim]) {
             nominal
             vertex_fit([2, 3, 4, 5])
             constrain_four_momentum
             chi2_cut 200
           }
alg_m02.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex; resolved in BOSS by the minimum chi2_4C.")
      .note(:vertex_fit, "Common production-vertex fit on all charged tracks.")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m02.with_decay_card(decay_card_m02).apply(sel_m02)
alg_m02.execute_on(scan_points + incMC_points + exMC_m02)

# =============================================================================
# Mode 03 : eta_c -> 2 (K+ K-)
# =============================================================================
decay_card_m03 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- K+ K-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m03 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_2KK"
  config.events        = 40_000
  config.decay_card    = decay_card_m03
  config.cross_section = :default
end

alg_m03 = Algorithm.new("EtaHcEtac2KK")
alg_m03.set_header(["EtaHcEtac2KKAlg/EtaHcEtac2KK.h"])
sel_m03 = Selection.new
sel_m03.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=2"
             nChrn     ">=2"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=3"
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion, :proton]
             nkp ">=1"; nkm ">=1"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=1"
           }
          .kinematic_fit([:gamma, :eta, :kp, :km, :kp, :km]) {
             nominal
             vertex_fit([2, 3, 4, 5])
             constrain_four_momentum
             chi2_cut 200
           }
alg_m03.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex; resolved in BOSS by the minimum chi2_4C.")
      .note(:vertex_fit, "Common production-vertex fit on all charged tracks.")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m03.with_decay_card(decay_card_m03).apply(sel_m03)
alg_m03.execute_on(scan_points + incMC_points + exMC_m03)

# =============================================================================
# Mode 04 : eta_c -> K+ K- pi+ pi-
# =============================================================================
decay_card_m04 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m04 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_KKPiPi"
  config.events        = 40_000
  config.decay_card    = decay_card_m04
  config.cross_section = :default
end

alg_m04 = Algorithm.new("EtaHcEtacKKPiPi")
alg_m04.set_header(["EtaHcEtacKKPiPiAlg/EtaHcEtacKKPiPi.h"])
sel_m04 = Selection.new
sel_m04.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=2"
             nChrn     ">=2"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=3"
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion, :proton]
             identify :pion, against: [:kaon, :proton]
             nkp ">=1"; nkm ">=1"
             npip ">=1"; npim ">=1"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=1"
           }
          .kinematic_fit([:gamma, :eta, :kp, :km, :pip, :pim]) {
             nominal
             vertex_fit([2, 3, 4, 5])
             constrain_four_momentum
             chi2_cut 200
           }
alg_m04.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex; resolved in BOSS by the minimum chi2_4C.")
      .note(:vertex_fit, "Common production-vertex fit on all charged tracks.")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information; the pi/K separation is resolved with the kinematic fit.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m04.with_decay_card(decay_card_m04).apply(sel_m04)
alg_m04.execute_on(scan_points + incMC_points + exMC_m04)

# =============================================================================
# Mode 05 : eta_c -> p pbar pi+ pi-
# =============================================================================
decay_card_m05 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 p+ anti-p- pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m05 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_ppbarpipi"
  config.events        = 40_000
  config.decay_card    = decay_card_m05
  config.cross_section = :default
end

alg_m05 = Algorithm.new("EtaHcEtacPPbarPiPi")
alg_m05.set_header(["EtaHcEtacPPbarPiPiAlg/EtaHcEtacPPbarPiPi.h"])
sel_m05 = Selection.new
sel_m05.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=2"
             nChrn     ">=2"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=3"
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :proton, against: [:kaon, :pion]
             identify :pion,   against: [:kaon, :proton]
             nprp ">=1"; nprm ">=1"
             npip ">=1"; npim ">=1"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=1"
           }
          .kinematic_fit([:gamma, :eta, :prp, :prm, :pip, :pim]) {
             nominal
             vertex_fit([2, 3, 4, 5])
             constrain_four_momentum
             chi2_cut 200
           }
alg_m05.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex; resolved in BOSS by the minimum chi2_4C.")
      .note(:vertex_fit, "Common production-vertex fit on all charged tracks.")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m05.with_decay_card(decay_card_m05).apply(sel_m05)
alg_m05.execute_on(scan_points + incMC_points + exMC_m05)

# =============================================================================
# Mode 06 : eta_c -> 3 (pi+ pi-)
# =============================================================================
decay_card_m06 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi+ pi- pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m06 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_3pipi"
  config.events        = 40_000
  config.decay_card    = decay_card_m06
  config.cross_section = :default
end

alg_m06 = Algorithm.new("EtaHcEtac3PiPi")
alg_m06.set_header(["EtaHcEtac3PiPiAlg/EtaHcEtac3PiPi.h"])
sel_m06 = Selection.new
sel_m06.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=3"
             nChrn     ">=3"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=3"
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :pion, against: [:kaon, :proton]
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=1"
           }
          .kinematic_fit([:gamma, :eta, :pip, :pim, :pip, :pim, :pip, :pim]) {
             nominal
             vertex_fit([2, 3, 4, 5, 6, 7])
             constrain_four_momentum
             chi2_cut 200
           }
alg_m06.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex; resolved in BOSS by the minimum chi2_4C.")
      .note(:vertex_fit, "Common production-vertex fit on all charged tracks.")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m06.with_decay_card(decay_card_m06).apply(sel_m06)
alg_m06.execute_on(scan_points + incMC_points + exMC_m06)

# =============================================================================
# Mode 07 : eta_c -> K+ K- 2 (pi+ pi-)
# =============================================================================
decay_card_m07 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- pi+ pi- pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m07 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_KK2pipi"
  config.events        = 40_000
  config.decay_card    = decay_card_m07
  config.cross_section = :default
end

alg_m07 = Algorithm.new("EtaHcEtacKK2PiPi")
alg_m07.set_header(["EtaHcEtacKK2PiPiAlg/EtaHcEtacKK2PiPi.h"])
sel_m07 = Selection.new
sel_m07.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=3"
             nChrn     ">=3"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=3"
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion, :proton]
             identify :pion, against: [:kaon, :proton]
             nkp ">=1"; nkm ">=1"
             npip ">=1"; npim ">=1"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=1"
           }
          .kinematic_fit([:gamma, :eta, :kp, :km, :pip, :pim, :pip, :pim]) {
             nominal
             vertex_fit([2, 3, 4, 5, 6, 7])
             constrain_four_momentum
             chi2_cut 200
           }
alg_m07.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex; resolved in BOSS by the minimum chi2_4C.")
      .note(:vertex_fit, "Common production-vertex fit on all charged tracks.")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information; the pi/K separation is resolved with the kinematic fit.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m07.with_decay_card(decay_card_m07).apply(sel_m07)
alg_m07.execute_on(scan_points + incMC_points + exMC_m07)

# =============================================================================
# Mode 08 : eta_c -> K+ K- pi0
# =============================================================================
decay_card_m08 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- pi0   PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m08 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_KKpi0"
  config.events        = 40_000
  config.decay_card    = decay_card_m08
  config.cross_section = :default
end

alg_m08 = Algorithm.new("EtaHcEtacKKPi0")
alg_m08.set_header(["EtaHcEtacKKPi0Alg/EtaHcEtacKKPi0.h"])
sel_m08 = Selection.new
sel_m08.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=1"    # K+
             nChrn     ">=1"    # K-
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=5"   # E1 gamma + primary eta + the two photons of pi0
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion, :proton]
             nkp ">=1"; nkm ">=1"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 25            # |M(gamma gamma) - m(pi0)| < 15 MeV/c^2
             npi0 ">=1"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=1"
           }
          .kinematic_fit([:gamma, :eta, :kp, :km, :pi0]) {
             nominal
             vertex_fit([2, 3])
             constrain_four_momentum
             chi2_cut 200
           }
alg_m08.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex, where chi2_1C is the mass constraint of the two photons of each pi0/eta; resolved in BOSS by the minimum chi2_4C.")
      .note(:vertex_fit, "Common production-vertex fit on the two charged tracks.")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:pi0_mass_window, "pi0 candidates are reconstructed from photon pairs with |M(gamma gamma) - m(pi0)| < 15 MeV/c^2, improved by a 1C mass-constrained fit (kalman_kinematic_fit).")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m08.with_decay_card(decay_card_m08).apply(sel_m08)
alg_m08.execute_on(scan_points + incMC_points + exMC_m08)

# =============================================================================
# Mode 09 : eta_c -> p pbar pi0
# =============================================================================
decay_card_m09 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 p+ anti-p- pi0   PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m09 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_ppbarpi0"
  config.events        = 40_000
  config.decay_card    = decay_card_m09
  config.cross_section = :default
end

alg_m09 = Algorithm.new("EtaHcEtacPPbarPi0")
alg_m09.set_header(["EtaHcEtacPPbarPi0Alg/EtaHcEtacPPbarPi0.h"])
sel_m09 = Selection.new
sel_m09.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=1"
             nChrn     ">=1"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=5"   # E1 gamma + primary eta + the two photons of pi0
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :proton, against: [:kaon, :pion]
             nprp ">=1"; nprm ">=1"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 25
             npi0 ">=1"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=1"
           }
          .kinematic_fit([:gamma, :eta, :prp, :prm, :pi0]) {
             nominal
             vertex_fit([2, 3])
             constrain_four_momentum
             chi2_cut 200
           }
alg_m09.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex; resolved in BOSS by the minimum chi2_4C.")
      .note(:vertex_fit, "Common production-vertex fit on the two charged tracks.")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:pi0_mass_window, "pi0 candidates: |M(gamma gamma) - m(pi0)| < 15 MeV/c^2 with a 1C mass-constrained fit.")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m09.with_decay_card(decay_card_m09).apply(sel_m09)
alg_m09.execute_on(scan_points + incMC_points + exMC_m09)

# =============================================================================
# Mode 10 : eta_c -> K_S0 K+ pi-
# =============================================================================
decay_card_m10 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K_S0 K+ pi-   PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m10 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_KSKpi"
  config.events        = 40_000
  config.decay_card    = decay_card_m10
  config.cross_section = :default
end

alg_m10 = Algorithm.new("EtaHcEtacKSKPi")
alg_m10.set_header(["EtaHcEtacKSKPiAlg/EtaHcEtacKSKPi.h"])
sel_m10 = Selection.new
sel_m10.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=2"    # K+ and the pi+ of K_S0
             nChrn     ">=2"    # pi- and the pi- of K_S0
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=3"
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion, :proton]
             nkp ">=1"
           }
          .secondary_vertex_fit([:pip, :pim]) {
             build_virtual_particle(:K_S0).by_minimizing_mass_difference
             remove_used_particle_from_candidate_list
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=1"
           }
          .kinematic_fit([:gamma, :eta, :K_S0, :kp, :pim]) {
             nominal
             vertex_fit([3, 4])     # only the prompt tracks; K_S0 daughters are excluded
             constrain_four_momentum
             chi2_cut 200
           }
alg_m10.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex; resolved in BOSS by the minimum chi2_4C.")
      .note(:ks0_selection, "K_S0 -> pi+ pi- candidates require |M(pi+ pi-) - m(K_S0)| < 20 MeV/c^2, a vertex fit constraining the two tracks to a common decay vertex, and a decay length greater than twice the vertex resolution; a secondary-vertex constraint between production and decay vertices rejects random pi+ pi- combinations. The POCA and production-vertex requirements are NOT applied to the K_S0 daughters.")
      .note(:vertex_fit, "Common production-vertex fit on the prompt charged tracks only (the K_S0 daughters are excluded).")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information; the momentum of the K+ is used to distinguish it from the pions.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m10.with_decay_card(decay_card_m10).apply(sel_m10)
alg_m10.execute_on(scan_points + incMC_points + exMC_m10)

# =============================================================================
# Mode 11 : eta_c -> K_S0 K+ pi- pi+ pi-
# =============================================================================
decay_card_m11 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K_S0 K+ pi- pi+ pi-   PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi-   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m11 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_KSKpipi"
  config.events        = 40_000
  config.decay_card    = decay_card_m11
  config.cross_section = :default
end

alg_m11 = Algorithm.new("EtaHcEtacKSKPiPi")
alg_m11.set_header(["EtaHcEtacKSKPiPiAlg/EtaHcEtacKSKPiPi.h"])
sel_m11 = Selection.new
sel_m11.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=3"    # K+, pi+ of K_S0 and pi+
             nChrn     ">=3"    # pi-, pi- of K_S0 and pi-
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=3"
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion, :proton]
             nkp ">=1"
           }
          .secondary_vertex_fit([:pip, :pim]) {
             build_virtual_particle(:K_S0).by_minimizing_mass_difference
             remove_used_particle_from_candidate_list
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=1"
           }
          .kinematic_fit([:gamma, :eta, :K_S0, :kp, :pim, :pip, :pim]) {
             nominal
             vertex_fit([3, 4, 5, 6])   # only the prompt tracks; K_S0 daughters are excluded
             constrain_four_momentum
             chi2_cut 200
           }
alg_m11.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex; resolved in BOSS by the minimum chi2_4C.")
      .note(:ks0_selection, "K_S0 -> pi+ pi- candidates require |M(pi+ pi-) - m(K_S0)| < 20 MeV/c^2, a common decay-vertex fit and a decay length greater than twice the vertex resolution, with the production/decay secondary-vertex constraint; the POCA and production-vertex requirements are not applied to the daughters.")
      .note(:vertex_fit, "Common production-vertex fit on the prompt charged tracks only (the K_S0 daughters are excluded).")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m11.with_decay_card(decay_card_m11).apply(sel_m11)
alg_m11.execute_on(scan_points + incMC_points + exMC_m11)

# =============================================================================
# Mode 12 : eta_c -> pi+ pi- eta
# =============================================================================
decay_card_m12 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- eta   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m12 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_pipieta"
  config.events        = 40_000
  config.decay_card    = decay_card_m12
  config.cross_section = :default
end

alg_m12 = Algorithm.new("EtaHcEtacPiPiEta")
alg_m12.set_header(["EtaHcEtacPiPiEtaAlg/EtaHcEtacPiPiEta.h"])
sel_m12 = Selection.new
sel_m12.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=1"
             nChrn     ">=1"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=5"   # E1 gamma + the two photons of each of the two eta
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :pion, against: [:kaon, :proton]
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=2"             # primary eta and the eta from eta_c
           }
          .kinematic_fit([:gamma, :eta, :pip, :pim, :eta]) {
             nominal
             vertex_fit([2, 3])
             constrain_four_momentum
             chi2_cut 200
           }
alg_m12.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex; resolved in BOSS by the minimum chi2_4C.")
      .note(:vertex_fit, "Common production-vertex fit on the two charged tracks.")
      .note(:two_eta_ambiguity, "The final state contains two eta mesons (the primary eta from e+e- -> eta h_c and the eta from eta_c decay); both are reconstructed from gamma gamma pairs with a 1C mass constraint to the nominal eta mass, and the pairing/combination is resolved by the combined chi2.")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:eta_mass_window, "eta candidates: |M(gamma gamma) - m(eta)| < 15 MeV/c^2 with a 1C mass-constrained fit.")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m12.with_decay_card(decay_card_m12).apply(sel_m12)
alg_m12.execute_on(scan_points + incMC_points + exMC_m12)

# =============================================================================
# Mode 13 : eta_c -> K+ K- eta
# =============================================================================
decay_card_m13 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- eta   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m13 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_KKeta"
  config.events        = 40_000
  config.decay_card    = decay_card_m13
  config.cross_section = :default
end

alg_m13 = Algorithm.new("EtaHcEtacKKEta")
alg_m13.set_header(["EtaHcEtacKKEtaAlg/EtaHcEtacKKEta.h"])
sel_m13 = Selection.new
sel_m13.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=1"
             nChrn     ">=1"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=5"   # E1 gamma + the two photons of each of the two eta
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion, :proton]
             nkp ">=1"; nkm ">=1"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=2"             # primary eta and the eta from eta_c
           }
          .kinematic_fit([:gamma, :eta, :kp, :km, :eta]) {
             nominal
             vertex_fit([2, 3])
             constrain_four_momentum
             chi2_cut 200
           }
alg_m13.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex; resolved in BOSS by the minimum chi2_4C.")
      .note(:vertex_fit, "Common production-vertex fit on the two charged tracks.")
      .note(:two_eta_ambiguity, "The final state contains two eta mesons (the primary eta and the eta from eta_c); both are reconstructed from gamma gamma pairs with a 1C mass constraint and the pairing is resolved by the combined chi2.")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:eta_mass_window, "eta candidates: |M(gamma gamma) - m(eta)| < 15 MeV/c^2 with a 1C mass-constrained fit.")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m13.with_decay_card(decay_card_m13).apply(sel_m13)
alg_m13.execute_on(scan_points + incMC_points + exMC_m13)

# =============================================================================
# Mode 14 : eta_c -> 2 (pi+ pi-) eta
# =============================================================================
decay_card_m14 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi+ pi- eta   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m14 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_2pipieta"
  config.events        = 40_000
  config.decay_card    = decay_card_m14
  config.cross_section = :default
end

alg_m14 = Algorithm.new("EtaHcEtac2PiPiEta")
alg_m14.set_header(["EtaHcEtac2PiPiEtaAlg/EtaHcEtac2PiPiEta.h"])
sel_m14 = Selection.new
sel_m14.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=2"
             nChrn     ">=2"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=5"   # E1 gamma + the two photons of each of the two eta
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :pion, against: [:kaon, :proton]
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=2"             # primary eta and the eta from eta_c
           }
          .kinematic_fit([:gamma, :eta, :pip, :pim, :pip, :pim, :eta]) {
             nominal
             vertex_fit([2, 3, 4, 5])
             constrain_four_momentum
             chi2_cut 200
           }
alg_m14.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex; resolved in BOSS by the minimum chi2_4C.")
      .note(:vertex_fit, "Common production-vertex fit on all charged tracks.")
      .note(:two_eta_ambiguity, "The final state contains two eta mesons (the primary eta and the eta from eta_c); both are reconstructed from gamma gamma pairs with a 1C mass constraint and the pairing is resolved by the combined chi2.")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:eta_mass_window, "eta candidates: |M(gamma gamma) - m(eta)| < 15 MeV/c^2 with a 1C mass-constrained fit.")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m14.with_decay_card(decay_card_m14).apply(sel_m14)
alg_m14.execute_on(scan_points + incMC_points + exMC_m14)

# =============================================================================
# Mode 15 : eta_c -> pi+ pi- pi0 pi0
# =============================================================================
decay_card_m15 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi0 pi0   PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m15 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_pipi2pi0"
  config.events        = 40_000
  config.decay_card    = decay_card_m15
  config.cross_section = :default
end

alg_m15 = Algorithm.new("EtaHcEtacPiPi2Pi0")
alg_m15.set_header(["EtaHcEtacPiPi2Pi0Alg/EtaHcEtacPiPi2Pi0.h"])
sel_m15 = Selection.new
sel_m15.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=1"
             nChrn     ">=1"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=7"   # E1 gamma + primary eta + the four photons of the two pi0
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :pion, against: [:kaon, :proton]
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 25
             npi0 ">=2"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=1"
           }
          .kinematic_fit([:gamma, :eta, :pip, :pim, :pi0, :pi0]) {
             nominal
             vertex_fit([2, 3])
             constrain_four_momentum
             chi2_cut 200
           }
alg_m15.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex, where chi2_1C is the mass constraint of the two photons of each pi0/eta; resolved in BOSS by the minimum chi2_4C.")
      .note(:vertex_fit, "Common production-vertex fit on the two charged tracks.")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:pi0_mass_window, "pi0 candidates: |M(gamma gamma) - m(pi0)| < 15 MeV/c^2 with a 1C mass-constrained fit.")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m15.with_decay_card(decay_card_m15).apply(sel_m15)
alg_m15.execute_on(scan_points + incMC_points + exMC_m15)

# =============================================================================
# Mode 16 : eta_c -> 2 (pi+ pi-) pi0 pi0
# =============================================================================
decay_card_m16 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta h_c   PHSP;
  Enddecay
  Decay h_c
  1.0000 gamma eta_c   PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma   PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi+ pi- pi0 pi0   PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay
  End
DECAYCARD

exMC_m16 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_etahc_etac_2pipi2pi0"
  config.events        = 40_000
  config.decay_card    = decay_card_m16
  config.cross_section = :default
end

alg_m16 = Algorithm.new("EtaHcEtac2PiPi2Pi0")
alg_m16.set_header(["EtaHcEtac2PiPi2Pi0Alg/EtaHcEtac2PiPi2Pi0.h"])
sel_m16 = Selection.new
sel_m16.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     ">=2"
             nChrn     ">=2"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    10.0
             nGam              ">=7"   # E1 gamma + primary eta + the four photons of the two pi0
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :pion, against: [:kaon, :proton]
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 25
             npi0 ">=2"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             chi2_cut 25
             neta ">=1"
           }
          .kinematic_fit([:gamma, :eta, :pip, :pim, :pip, :pim, :pi0, :pi0]) {
             nominal
             vertex_fit([2, 3, 4, 5])
             constrain_four_momentum
             chi2_cut 200
           }
alg_m16.note(:kinematic_fit_chi2, "The paper requires chi2_4C < 25; the loose default (200) is applied in BOSS and the published cut is applied in ROOT.")
      .note(:combination_selection, "Candidate minimising chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex; resolved in BOSS by the minimum chi2_4C.")
      .note(:vertex_fit, "Common production-vertex fit on all charged tracks.")
      .note(:e1_photon_energy, "E1 photon energy in (400, 600) MeV, enlarged to (350, 650) MeV for sqrt(s) > 4.416 GeV (ROOT level).")
      .note(:eta_recoil_mass, "M_recoil(eta) in (3480, 3600) MeV/c^2 (ROOT level).")
      .note(:eta_c_mass_window, "Hadronic invariant mass in (2940, 3020) MeV/c^2; the eta candidate giving the eta_c mass closest to the nominal value is selected (ROOT level).")
      .note(:pi0_mass_window, "pi0 candidates: |M(gamma gamma) - m(pi0)| < 15 MeV/c^2 with a 1C mass-constrained fit.")
      .note(:pid_selection, "PID confidence levels from TOF and dE/dx combined with kinematic-fit information.")
      .note(:born_cross_section, "sigma^Born = N_obs / (L (1+delta) |1+Pi|^2 B(eta->gamma gamma) B(h_c->gamma eta_c) sum_i eps_i B_i), with (1+delta) and |1+Pi|^2 obtained iteratively.")
      .note(:fit_model, "Simultaneous unbinned maximum-likelihood fit of the 16 eta recoil mass spectra (ROOT level).")
alg_m16.with_decay_card(decay_card_m16).apply(sel_m16)
alg_m16.execute_on(scan_points + incMC_points + exMC_m16)
