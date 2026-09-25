# BOSS DSL specification for psi(3686) -> gamma eta', gamma eta, gamma pi0
# arXiv: 1708.03103v1 (BESIII Collaboration)

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# =============================================================================
# Decay cards
# =============================================================================

# psi(3686) -> gamma eta', eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_gamma_etap_pipiEta = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta'                 HELAMP 1.0 0.0 -1.0 0.0;
  Enddecay

  Decay eta'
  1.000 pi+ pi- eta                ETA_DALITZ;
  Enddecay

  Decay eta
  1.000 gamma gamma                PHSP;
  Enddecay

  End
DECAYCARD

# psi(3686) -> gamma eta', eta' -> pi0 pi0 eta, eta -> gamma gamma
decay_card_gamma_etap_pi0pi0Eta = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta'                 HELAMP 1.0 0.0 -1.0 0.0;
  Enddecay

  Decay eta'
  1.000 pi0 pi0 eta                PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma                PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma                PHSP;
  Enddecay

  End
DECAYCARD

# psi(3686) -> gamma eta, eta -> pi+ pi- pi0
decay_card_gamma_eta_pipipi0 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta                  HELAMP 1.0 0.0 -1.0 0.0;
  Enddecay

  Decay eta
  1.000 pi+ pi- pi0                ETA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma                PHSP;
  Enddecay

  End
DECAYCARD

# psi(3686) -> gamma eta, eta -> pi0 pi0 pi0
decay_card_gamma_eta_3pi0 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta                  HELAMP 1.0 0.0 -1.0 0.0;
  Enddecay

  Decay eta
  1.000 pi0 pi0 pi0                PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma                PHSP;
  Enddecay

  End
DECAYCARD

# psi(3686) -> gamma pi0, pi0 -> gamma gamma
decay_card_gamma_pi0 = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma pi0                  HELAMP 1.0 0.0 -1.0 0.0;
  Enddecay

  Decay pi0
  1.000 gamma gamma                PHSP;
  Enddecay

  End
DECAYCARD

# =============================================================================
# Exclusive MC samples
# =============================================================================

exMC_gamma_etap_pipiEta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "gamma_etap_pipiEta"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_gamma_etap_pipiEta
  config.cross_section   = :default
end

exMC_gamma_etap_pi0pi0Eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "gamma_etap_pi0pi0Eta"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_gamma_etap_pi0pi0Eta
  config.cross_section   = :default
end

exMC_gamma_eta_pipipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "gamma_eta_pipipi0"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_gamma_eta_pipipi0
  config.cross_section   = :default
end

exMC_gamma_eta_3pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "gamma_eta_3pi0"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_gamma_eta_3pi0
  config.cross_section   = :default
end

exMC_gamma_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "gamma_pi0"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_gamma_pi0
  config.cross_section   = :default
end

# =============================================================================
# A1) psi(3686) -> gamma eta', eta' -> pi+ pi- eta, eta -> gamma gamma
# Final state: gamma pi+ pi- gamma gamma  (3 photons + 2 charged pions)
# =============================================================================
alg_gamma_etap_pipiEta = Algorithm.new("PsipGammaEtapPipPimEta")
alg_gamma_etap_pipiEta.set_header(["PsipGammaEtapPipPimEtaAlg/PsipGammaEtapPipPimEta.h"])
                     .set_constant({"ECMS" => [:double, 3.686]})

sel_gamma_etap_pipiEta = Selection.new
sel_gamma_etap_pipiEta.select_track {
                        cos_theta 0.93
                        Vz        10.0
                        Vr        1.0
                        nChrp     "==1"
                        nChrn     "==1"
                        nNet      "==0"
                      }
                     .select_photon {
                        energyThreshold_b  0.025
                        energyThreshold_e  0.050
                        angle_to_track     10.0
                        tdc_emc_start      0
                        tdc_emc_end        14
                        nGam               ">=3"
                      }
                     .pid(method: :probability) {
                        prob_cut 0.001
                        identify :pion, against: [:kaon, :proton]
                      }
                     .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {
                        nominal
                        constrain_four_momentum
                        chi2_cut 80
                      }

alg_gamma_etap_pipiEta.with_decay_card(decay_card_gamma_etap_pipiEta).apply(sel_gamma_etap_pipiEta)
alg_gamma_etap_pipiEta.execute_on([psip_data, psip_incMC, exMC_gamma_etap_pipiEta])

# =============================================================================
# A2) psi(3686) -> gamma eta', eta' -> pi0 pi0 eta, eta -> gamma gamma
# Final state: 7 photons, no charged tracks
# =============================================================================
alg_gamma_etap_pi0pi0Eta = Algorithm.new("PsipGammaEtapPi0Pi0Eta")
alg_gamma_etap_pi0pi0Eta.set_header(["PsipGammaEtapPi0Pi0EtaAlg/PsipGammaEtapPi0Pi0Eta.h"])
                       .set_constant({"ECMS" => [:double, 3.686]})

sel_gamma_etap_pi0pi0Eta = Selection.new
sel_gamma_etap_pi0pi0Eta.select_track {
                          cos_theta 0.93
                          Vz        10.0
                          Vr        1.0
                          nChrp     "==0"
                          nChrn     "==0"
                          nNet      "==0"
                        }
                       .select_photon {
                          energyThreshold_b  0.025
                          energyThreshold_e  0.050
                          angle_to_track     10.0
                          nGam               ">=7"
                        }
                       .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) {
                          nominal
                          constrain_four_momentum
                          chi2_cut 80
                        }
                        .note(:best_combination_chi2M,
                             "For eta' -> pi0 pi0 eta with multiple photon combinations, the " \
                             "combination minimising chi2_M = sum_i (M(gamma_i gamma_j) - M_pi0)^2/sigma_pi0^2 " \
                             "for the two pi0 pairs plus (M(gamma_5 gamma_6) - M_eta)^2/sigma_eta^2 " \
                             "is retained. sigma_pi0 = 4.8 MeV/c^2, sigma_eta = 8.7 MeV/c^2.")

alg_gamma_etap_pi0pi0Eta.with_decay_card(decay_card_gamma_etap_pi0pi0Eta).apply(sel_gamma_etap_pi0pi0Eta)
alg_gamma_etap_pi0pi0Eta.execute_on([psip_data, psip_incMC, exMC_gamma_etap_pi0pi0Eta])

# =============================================================================
# B1) psi(3686) -> gamma eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma
# Final state: gamma pi+ pi- gamma gamma  (3 photons + 2 charged pions)
# No angle-to-track cut for photons (higher-momentum eta -> pi0 photons close to charged pions)
# =============================================================================
alg_gamma_eta_pipipi0 = Algorithm.new("PsipGammaEtaPipPimPi0")
alg_gamma_eta_pipipi0.set_header(["PsipGammaEtaPipPimPi0Alg/PsipGammaEtaPipPimPi0.h"])
                    .set_constant({"ECMS" => [:double, 3.686]})

sel_gamma_eta_pipipi0 = Selection.new
sel_gamma_eta_pipipi0.select_track {
                       cos_theta 0.93
                       Vz        10.0
                       Vr        1.0
                       nChrp     "==1"
                       nChrn     "==1"
                       nNet      "==0"
                     }
                    .select_photon {
                       energyThreshold_b  0.025
                       energyThreshold_e  0.050
                       tdc_emc_start      0
                       tdc_emc_end        14
                       nGam               ">=3"
                     }
                    .pid(method: :probability) {
                       prob_cut 0.001
                       identify :pion, against: [:kaon, :proton]
                     }
                    .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {
                       nominal
                       constrain_four_momentum
                       chi2_cut 80
                     }

alg_gamma_eta_pipipi0.with_decay_card(decay_card_gamma_eta_pipipi0).apply(sel_gamma_eta_pipipi0)
alg_gamma_eta_pipipi0.execute_on([psip_data, psip_incMC, exMC_gamma_eta_pipipi0])

# =============================================================================
# B2) psi(3686) -> gamma eta, eta -> pi0 pi0 pi0
# Final state: 7 photons, no charged tracks
# =============================================================================
alg_gamma_eta_3pi0 = Algorithm.new("PsipGammaEta3Pi0")
alg_gamma_eta_3pi0.set_header(["PsipGammaEta3Pi0Alg/PsipGammaEta3Pi0.h"])
                 .set_constant({"ECMS" => [:double, 3.686]})

sel_gamma_eta_3pi0 = Selection.new
sel_gamma_eta_3pi0.select_track {
                    cos_theta 0.93
                    Vz        10.0
                    Vr        1.0
                    nChrp     "==0"
                    nChrn     "==0"
                    nNet      "==0"
                  }
                 .select_photon {
                    energyThreshold_b  0.025
                    energyThreshold_e  0.050
                    nGam               ">=7"
                  }
                 .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) {
                    nominal
                    constrain_four_momentum
                    chi2_cut 80
                  }

alg_gamma_eta_3pi0.with_decay_card(decay_card_gamma_eta_3pi0).apply(sel_gamma_eta_3pi0)
alg_gamma_eta_3pi0.execute_on([psip_data, psip_incMC, exMC_gamma_eta_3pi0])

# =============================================================================
# C) psi(3686) -> gamma pi0, pi0 -> gamma gamma
# Final state: 3 photons in barrel EMC, no charged tracks
# =============================================================================
alg_gamma_pi0 = Algorithm.new("PsipGammaPi0")
alg_gamma_pi0.set_header(["PsipGammaPi0Alg/PsipGammaPi0.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

sel_gamma_pi0 = Selection.new
sel_gamma_pi0.select_track {
               cos_theta 0.93
               Vz        10.0
               Vr        1.0
               nChrp     "==0"
               nChrn     "==0"
               nNet      "==0"
             }
            .select_photon {
               energyThreshold_b  0.025
               energyThreshold_e  0.050
               nGam               "==3"
             }
            .kinematic_fit([:gamma, :gamma, :gamma]) {
               nominal
               constrain_four_momentum
               chi2_cut 40
             }
            .note(:barrel_only_photons,
                 "Only photons in the EMC barrel region (|cos(theta)| < 0.8) are accepted, " \
                 "to suppress the QED background e+e- -> gamma gamma (gamma_ISR).")
            .note(:pi0_helicity_angle_cut,
                 "|cos(theta_hel)| < 0.7 for the pi0, where theta_hel is the angle between " \
                 "the more energetic photon momentum in the pi0 rest frame and the pi0 " \
                 "momentum in the psi(3686) rest frame. Suppresses QED background.")
            .note(:mdc_hits_gamma_conversion_veto,
                 "Fewer than 8 MDC hits are required in the region between the two radial " \
                 "lines connecting the IP and the two shower positions in the EMC, to reject " \
                 "gamma conversion background from e+e- -> gamma gamma (gamma_ISR) with one " \
                 "photon converting to an e+e- pair.")

alg_gamma_pi0.with_decay_card(decay_card_gamma_pi0).apply(sel_gamma_pi0)
alg_gamma_pi0.execute_on([psip_data, psip_incMC, exMC_gamma_pi0])
