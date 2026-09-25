### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# --- Decay cards ---
decay_card_phiphi_eta = <<~DECAYCARD
  Alias phi_a phi
  Alias phi_b phi

  Decay psi(3686)
  1.000 gamma chi_c1                       PHSP;
  Enddecay

  Decay chi_c1
  1.000 phi_a phi_b eta                    PHSP;
  Enddecay

  Decay phi_a
  1.000 K+ K-                              PHSP;
  Enddecay

  Decay phi_b
  1.000 K+ K-                              PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma                        PHSP;
  Enddecay

  End
DECAYCARD

decay_card_phiphi_etap_mode1 = <<~DECAYCARD
  Alias phi_a phi
  Alias phi_b phi

  Decay psi(3686)
  1.000 gamma chi_c1                       PHSP;
  Enddecay

  Decay chi_c1
  1.000 phi_a phi_b eta'                   PHSP;
  Enddecay

  Decay phi_a
  1.000 K+ K-                              PHSP;
  Enddecay

  Decay phi_b
  1.000 K+ K-                              PHSP;
  Enddecay

  Decay eta'
  1.000 pi+ pi- gamma                      PHSP;
  Enddecay

  End
DECAYCARD

decay_card_phiphi_etap_mode2 = <<~DECAYCARD
  Alias phi_a phi
  Alias phi_b phi

  Decay psi(3686)
  1.000 gamma chi_c1                       PHSP;
  Enddecay

  Decay chi_c1
  1.000 phi_a phi_b eta'                   PHSP;
  Enddecay

  Decay phi_a
  1.000 K+ K-                              PHSP;
  Enddecay

  Decay phi_b
  1.000 K+ K-                              PHSP;
  Enddecay

  Decay eta'
  1.000 pi+ pi- eta                        PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma                        PHSP;
  Enddecay

  End
DECAYCARD

decay_card_phi_KKeta = <<~DECAYCARD
  Decay psi(3686)
  1.000 gamma chi_c1                       PHSP;
  Enddecay

  Decay chi_c1
  1.000 phi K+ K- eta                      PHSP;
  Enddecay

  Decay phi
  1.000 K+ K-                              PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma                        PHSP;
  Enddecay

  End
DECAYCARD

# --- Exclusive MC samples ---
exMC_phiphi_eta = DatasetManager.create_exclusive_mc do |c|
  c.sample_name    = "chi_cJ_phiphi_eta"
  c.related_dataset = psip_data
  c.events         = 1_000_000
  c.decay_card     = decay_card_phiphi_eta
  c.cross_section  = :default
end
exMC_phiphi_eta.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_phiphi_etap_m1 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name    = "chi_cJ_phiphi_etap_pipiG"
  c.related_dataset = psip_data
  c.events         = 1_000_000
  c.decay_card     = decay_card_phiphi_etap_mode1
  c.cross_section  = :default
end
exMC_phiphi_etap_m1.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_phiphi_etap_m2 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name    = "chi_cJ_phiphi_etap_pipiEta"
  c.related_dataset = psip_data
  c.events         = 1_000_000
  c.decay_card     = decay_card_phiphi_etap_mode2
  c.cross_section  = :default
end
exMC_phiphi_etap_m2.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_phi_KKeta = DatasetManager.create_exclusive_mc do |c|
  c.sample_name    = "chi_cJ_phi_KKeta"
  c.related_dataset = psip_data
  c.events         = 1_000_000
  c.decay_card     = decay_card_phi_KKeta
  c.cross_section  = :default
end
exMC_phi_KKeta.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection ###

# --------- Channel 1: chi_cJ -> phi phi eta ---------
alg1 = Algorithm.new("ChicJPhiPhiEta")
alg1.set_header(["ChicJPhiPhiEtaAlg/ChicJPhiPhiEta.h"])
    .set_constant({ "ECMS" => [:double, 3.686] })
    .set_alias({ "std::vector<double>" => "Vdouble" })

sel1 = Selection.new
sel1.select_track {
       cos_theta 0.93
       Vz 10.0
       Vr 1.0
       nChrp "==2"
       nChrn "==2"
       nNet  "==0"
     }
    .select_photon {
       tdc_emc_start 0
       tdc_emc_end 14
       angle_to_track 10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam ">=3"
     }
    .pid(method: :probability) {
       prob_cut 0.001
       identify :kaon,  against: [:pion, :proton]
       identify :pion,  against: [:kaon, :proton]
       nkp  ">=2"
       nkm  ">=2"
       npip "==0"
       npim "==0"
     }
    .kinematic_fit([:gamma, :gamma, :gamma, :kp, :kp, :km, :km]) {
       nominal
       constrain_four_momentum
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
       chi2_cut 200
    }

alg1.note(:phi_mass_window,
          "Two phi candidates selected by minimising DeltaM^2 = (M(K+K-)i - m_phi)^2 + " \
          "(M(K+K-)j - m_phi)^2; phi signal region 1.005 <= M(K+K-) <= 1.035 GeV/c^2 " \
          "(applied at ROOT stage).")
    .note(:eta_signal_region,
          "eta signal region 0.51 <= M(gg) <= 0.57 GeV/c^2 with the gamma-gamma pair " \
          "closest to nominal eta mass; applied at ROOT stage.")
    .note(:chi2_4c_cut,
          "chi2_4C < 38 required in the paper; loose 200 used here (Rule T3), tight cut in ROOT.")
alg1.with_decay_card(decay_card_phiphi_eta).apply(sel1)
alg1.execute_on([psip_data, psip_incMC, exMC_phiphi_eta])

# --------- Channel 2: chi_cJ -> phi phi eta' (Mode I: eta' -> pi+ pi- gamma) ---------
alg2 = Algorithm.new("ChicJPhiPhiEtapModeI")
alg2.set_header(["ChicJPhiPhiEtapModeIAlg/ChicJPhiPhiEtapModeI.h"])
    .set_constant({ "ECMS" => [:double, 3.686] })
    .set_alias({ "std::vector<double>" => "Vdouble" })

sel2 = Selection.new
sel2.select_track {
       cos_theta 0.93
       Vz 10.0
       Vr 1.0
       nChrp "==3"
       nChrn "==3"
       nNet  "==0"
     }
    .select_photon {
       tdc_emc_start 0
       tdc_emc_end 14
       angle_to_track 10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam ">=2"
     }
    .pid(method: :probability) {
       prob_cut 0.001
       identify :kaon,  against: [:pion, :proton]
       identify :pion,  against: [:kaon, :proton]
       nkp  ">=2"
       nkm  ">=2"
       npip ">=1"
       npim ">=1"
     }
    .kinematic_fit([:gamma, :gamma, :kp, :kp, :km, :km, :pip, :pim]) {
       nominal
       constrain_four_momentum
       chi2_cut 200
    }

alg2.note(:etap_signal_region,
          "eta' signal region |M(pi+ pi- gamma) - 0.958| < 0.015 GeV/c^2 (ROOT stage).")
    .note(:pi0_veto,
          "|M(gg) - m_pi0| > 0.015 GeV/c^2 to suppress psi(3686) -> K+K+K-K- pi+ pi- pi0 (ROOT stage).")
    .note(:phi_mass_window,
          "Two phi from K+K- pairs with 3-sigma window 1.005-1.035 GeV/c^2 (ROOT stage).")
alg2.with_decay_card(decay_card_phiphi_etap_mode1).apply(sel2)
alg2.execute_on([psip_data, psip_incMC, exMC_phiphi_etap_m1])

# --------- Channel 3: chi_cJ -> phi phi eta' (Mode II: eta' -> pi+ pi- eta, eta -> gg) ---------
alg3 = Algorithm.new("ChicJPhiPhiEtapModeII")
alg3.set_header(["ChicJPhiPhiEtapModeIIAlg/ChicJPhiPhiEtapModeII.h"])
    .set_constant({ "ECMS" => [:double, 3.686] })
    .set_alias({ "std::vector<double>" => "Vdouble" })

sel3 = Selection.new
sel3.select_track {
       cos_theta 0.93
       Vz 10.0
       Vr 1.0
       nChrp "==3"
       nChrn "==3"
       nNet  "==0"
     }
    .select_photon {
       tdc_emc_start 0
       tdc_emc_end 14
       angle_to_track 10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam ">=3"
     }
    .pid(method: :probability) {
       prob_cut 0.001
       identify :kaon,  against: [:pion, :proton]
       identify :pion,  against: [:kaon, :proton]
       nkp  ">=2"
       nkm  ">=2"
       npip ">=1"
       npim ">=1"
     }
    .kinematic_fit([:gamma, :gamma, :gamma, :kp, :kp, :km, :km, :pip, :pim]) {
       nominal
       constrain_four_momentum
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
       chi2_cut 200
    }

alg3.note(:etap_signal_region,
          "eta' signal region |M(pi+ pi- eta) - 0.958| < 0.015 GeV/c^2 (ROOT stage).")
    .note(:phi_mass_window,
          "Two phi from K+K- pairs with 3-sigma window 1.005-1.035 GeV/c^2 (ROOT stage).")
alg3.with_decay_card(decay_card_phiphi_etap_mode2).apply(sel3)
alg3.execute_on([psip_data, psip_incMC, exMC_phiphi_etap_m2])

# --------- Channel 4: chi_cJ -> phi K+ K- eta (non-phi K+K-) ---------
alg4 = Algorithm.new("ChicJPhiKKEta")
alg4.set_header(["ChicJPhiKKEtaAlg/ChicJPhiKKEta.h"])
    .set_constant({ "ECMS" => [:double, 3.686] })
    .set_alias({ "std::vector<double>" => "Vdouble" })

sel4 = Selection.new
sel4.select_track {
       cos_theta 0.93
       Vz 10.0
       Vr 1.0
       nChrp "==2"
       nChrn "==2"
       nNet  "==0"
     }
    .select_photon {
       tdc_emc_start 0
       tdc_emc_end 14
       angle_to_track 10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam ">=3"
     }
    .pid(method: :probability) {
       prob_cut 0.001
       identify :kaon,  against: [:pion, :proton]
       nkp  ">=2"
       nkm  ">=2"
     }
    .kinematic_fit([:gamma, :gamma, :gamma, :kp, :kp, :km, :km]) {
       nominal
       constrain_four_momentum
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
       chi2_cut 200
    }

alg4.note(:phi_signal_region,
          "One K+K- pair defines phi (1.005-1.035 GeV/c^2), the other K+K- has M > 1.09 GeV/c^2 " \
          "(non-phi region); applied at ROOT stage.")
    .note(:jpsi_veto,
          "|M(K+ K+ K- K-) - M_Jpsi| > 0.015 GeV/c^2 to suppress psi(3686) -> eta J/psi, " \
          "J/psi -> K+K+K-K- (ROOT stage).")
alg4.with_decay_card(decay_card_phi_KKeta).apply(sel4)
alg4.execute_on([psip_data, psip_incMC, exMC_phi_KKeta])
