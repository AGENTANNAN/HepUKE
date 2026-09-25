# BESIII analysis: Search for eta_c(2S) -> K+ K- eta' and updated chi_cJ -> K+ K- eta'
# arXiv: 2407.17184v1
# psi(3686) -> gamma eta_c(2S)/chi_cJ with eta' -> pi+ pi- gamma / pi+ pi- eta (eta -> gamma gamma)
# 2.712e9 psi(3686) events + continuum at 3.650 GeV

### Datasets ###
data_3686 = DatasetManager.real_data.find("709_3686")
data_3650 = DatasetManager.real_data.find("709_3650")

incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")
incMC_3650 = DatasetManager.inclusive_mc.find("709_3650")

all_datasets = [data_3686, data_3650]
all_incMC = [incMC_3686, incMC_3650]

### Decay card for psi(3686) -> gamma eta_c(2S)/chi_cJ, eta' -> pi+ pi- gamma ###
decay_card_etap_pipig = <<~DECAYCARD
    Decay psi(3686)
    1.000 gamma eta_c(2S) VSP_PWAVE 1.0 0.0 1.0 0.0 1.0;
    Enddecay

    Decay eta_c(2S)
    1.000 K+ K- eta' PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- gamma PHSP;
    Enddecay

    End
DECAYCARD

### Decay card for psi(3686) -> gamma chi_cJ, eta' -> pi+ pi- eta (eta -> gamma gamma) ###
decay_card_etap_pipieta = <<~DECAYCARD
    Decay psi(3686)
    1.000 gamma chi_c1 VSP_PWAVE 1.0 -0.333 1.0 0.0 1.0;
    Enddecay

    Decay chi_c1
    1.000 K+ K- eta' PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC for Mode I: eta' -> pi+ pi- gamma ###
exMC_etap_pipig = DatasetManager.create_exclusive_mc do |config|
  config.sample_name   = "psip_gamma_etac_etap_pipig"
  config.events        = 500_000
  config.decay_card    = decay_card_etap_pipig
  config.related_dataset = data_3686
  config.cross_section = :default
end

### Exclusive MC for Mode II: eta' -> pi+ pi- eta, eta -> gamma gamma ###
exMC_etap_pipieta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name   = "psip_gamma_chic_etap_pipieta"
  config.events        = 500_000
  config.decay_card    = decay_card_etap_pipieta
  config.related_dataset = data_3686
  config.cross_section = :default
end

### Algorithm I: eta' -> pi+ pi- gamma ###
alg_modeI = Algorithm.new("PsipGammaEtacEtapPiPiGam")
alg_modeI.set_header(["PsipGammaEtacEtapPiPiGamAlg/PsipGammaEtacEtapPiPiGam.h"])
  .set_constant({"ECMS" => [:double, 3.686]})
  .note(:psip_sample, "2.712e9 psi(3686) events. Continuum data at 3.650 GeV (401 pb^-1) used for background estimation.")
  .note(:kinematic_fit_4c, "4C kinematic fit constraining final state to initial e+e- system. chi2_4C < 25 for mode I.")
  .note(:pid_combined, "Combined PID using dE/dx and TOF. chi2_PID minimized together with chi2_4C to select best track species assignment.")
  .note(:etap_window, "eta' signal region: M(pi+ pi- gamma) in (0.94, 0.97) GeV/c^2. Sidebands: [0.8,0.9] and [1.0,1.1] for bg estimation.")
  .note(:mass_vetoes_modeI, "pi0 veto: M(gamma gamma) outside (0.122,0.146); eta veto: M(gamma gamma) outside (0.526,0.566); J/psi veto: M(K+K-pi+pi-) outside (3.08,3.12); chi_cJ vetoes: M(K+K-pi+pi-) outside (3.37,3.43)/(3.48,3.51)/(3.53,3.54); pi+pi-J/psi veto: M_recoil(pi+pi-) outside (3.091,3.103); phi veto: M(K+K-) < 1.03 GeV/c^2")
  .note(:kinematic_fit_3c, "3C kinematic fit applied after selection: radiative photon energy not used as input. M(K+K-eta')^3C used for signal extraction via simultaneous unbinned ML fit.")
  .note(:signal_extraction, "Simultaneous fit to M(K+K-eta')^3C distributions from both eta' modes. Components: eta_c(2S)/chi_cJ signals (E_gamma^3 * damping * BW * resolution * efficiency), psi(3686)->K+K-eta' peaking bg, non-eta' bg (sideband), continuum, remaining bg (2nd-order polynomial).")
  .note(:fsr_correction, "FSR correction factor f_FSR = 1.39 +/- 0.08 +/- 0.04 applied to psi(3686)->K+K-eta' MC shape.")
  .note(:continuum_shift, "Continuum M(K+K-eta')^3C shifted from 3.650 to 3.686 GeV via m_shifted -> a(m - m0) + m0 with a=1.021, m0=1.945 GeV/c^2.")

sel_modeI = Selection.new
sel_modeI.select_track {
            cos_theta 0.93
            Vz 10.0
            Vr 1.0
            nChrp ">=2"
            nChrn ">=2"
          }
         .select_photon {
            tdc_emc_start 0
            tdc_emc_end 700
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam ">=2"
         }
         .kinematic_fit([:kp, :km, :pip, :pim, :gamma]) {
            nominal
            constrain_four_momentum
            chi2_cut 25
         }

alg_modeI.with_decay_card(decay_card_etap_pipig).apply(sel_modeI)
alg_modeI.execute_on(all_datasets + all_incMC + [exMC_etap_pipig])

### Algorithm II: eta' -> pi+ pi- eta, eta -> gamma gamma ###
alg_modeII = Algorithm.new("PsipGammaChicEtapPiPiEta")
alg_modeII.set_header(["PsipGammaChicEtapPiPiEtaAlg/PsipGammaChicEtapPiPiEta.h"])
  .set_constant({"ECMS" => [:double, 3.686]})
  .note(:psip_sample, "2.712e9 psi(3686) events. Continuum data at 3.650 GeV.")
  .note(:kinematic_fit_4c, "4C kinematic fit. chi2_4C < 20 for mode II.")
  .note(:kalman_eta, "1C kinematic fit for eta->gamma gamma candidate with M(gamma gamma) constrained to nominal eta mass. Best eta candidate from all photon pairs.")
  .note(:pid_combined, "Combined PID: chi2_tot = chi2_1C + chi2_4C + sum chi2_PID minimized to select best track species.")
  .note(:etap_window, "eta' signal region: M(pi+ pi- eta) in (0.93, 0.98) GeV/c^2. Sidebands: [0.85,0.90] and [1.0,1.1].")
  .note(:eta_window, "eta -> gamma gamma signal region: M(gamma gamma) in (0.52, 0.57) GeV/c^2.")
  .note(:mass_vetoes_modeII, "pi0 veto: M(gamma gamma_eta) outside (0.118,0.150); eta J/psi veto: M_recoil(eta) outside (3.070,3.133); phi veto: M(K+K-) < 1.03 GeV/c^2")
  .note(:kinematic_fit_3c, "3C kinematic fit applied after selection. M(K+K-eta')^3C used for simultaneous fit.")
  .note(:signal_extraction, "Simultaneous unbinned ML fit to M(K+K-eta')^3C for both eta' modes. See mode I notes for details.")

sel_modeII = Selection.new
sel_modeII.select_track {
             cos_theta 0.93
             Vz 10.0
             Vr 1.0
             nChrp ">=2"
             nChrn ">=2"
           }
          .select_photon {
             tdc_emc_start 0
             tdc_emc_end 700
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             nGam ">=3"
           }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
             neta ">=1"
           }
          .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {
             nominal
             constrain_four_momentum
             chi2_cut 20
           }

alg_modeII.with_decay_card(decay_card_etap_pipieta).apply(sel_modeII)
alg_modeII.execute_on(all_datasets + all_incMC + [exMC_etap_pipieta])