# ============================================================================
# Search for lepton number violating decays of D_s+ -> h- h0 e+ e+
# at c.m. energies 4.128-4.226 GeV (7.33 fb^-1). ST method.
# 2410.02421v2  —  BESIII
# ============================================================================

# Shared datasets — representative energy points (largest: 4.178 GeV, BOSS 703)
ds_data_4180 = DatasetManager.real_data.find("703_4180")
ds_incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")

# ============================================================================
# Mode I: D_s+ -> phi pi- e+ e+  (phi -> K+ K-)
# ============================================================================

decay_card_phi_pi = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.000 phi pi- e+ e+ PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

exMC_phi_pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Ds_phi_pi_ee"
  config.related_dataset = ds_data_4180
  config.events = 100_000
  config.decay_card = decay_card_phi_pi
  config.cross_section = :default
end

alg_phi_pi = Algorithm.new("DsToPhiPiEE")
alg_phi_pi.set_header(["DsToPhiPiEEAlg/DsToPhiPiEE.h"])
          .set_constant({ "ECMS" => [:double, 4.178] })

sel_phi_pi = Selection.new
sel_phi_pi.select_track {
            cos_theta 0.93
            Vz 10.0
            Vxy 1.0
            nChrp ">=3"
            nChrn ">=2"
            nNet "==1"
          }
          .select_photon {
            tdc_emc_start 0
            tdc_emc_end 700
            angle_to_track 10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion]
            identify :pion, against: [:kaon]
            identify :positron, against: [:kaon, :pion]
            nkp "==1"
            nkm "==1"
            npim "==1"
            nep "==2"
          }
          .kinematic_fit([:kp, :km, :pim, :ep]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_phi_pi
  .note(:signal_mode, "D_s+ -> phi pi- e+ e+ with phi -> K+ K-")
  .note(:production, "e+e- -> D_s*+- D_s-+ at c.m. energies 4.128-4.226 GeV; total 7.33 fb^-1")
  .note(:multi_energy, "Data at 8 energy points: 705_4130(4.128), 705_4160(4.157), 703_4180(4.178), 703_4190(4.189), 703_4200(4.199), 703_4210(4.209), 703_4220(4.219), 703_4230(4.226); representative MC at 4.178 GeV")
  .note(:e_pid, "e+: L_e/(L_e+L_pi+L_K) > 0.8, L_e > 0.001, E>25MeV in EMC, E/p > 0.8c (p>400MeV) or > 0.7c (p<400MeV), recovery of photons within 5deg cone. For phi pi mode: both e+ must satisfy PID criterion. For phi K mode: at least one e+ must satisfy PID criterion")
  .note(:phi_selection, "phi: M_K+K- in [1.00, 1.05] GeV/c^2; K+K- pair closest to nominal phi mass selected if multiple combinations")
  .note(:background_veto, "cos(theta_e+e-) < 0.95 (gamma conv., qqbar, Bhabha); cos(theta_e+pi-) < 0.98 (e/pi mis-ID); L_pi/sigma_L_pi < 3 (K_S0 veto for modes with pi- and without K_S0); L/sigma_L > 2 (fake K_S0 veto for modes with K_S0); E(pi0) > 0.17 GeV (soft pi0 veto)")
  .note(:signal_region, "M_rec and Delta_M 2D signal region optimized per energy point and decay mode via Punzi FOM (see Table 3); D_s+ -> phi K- e+ e+ mode exempt from M_rec/Delta_M requirement")
  .note(:e_ep_dedx, "Optimized E/p and chi2_dE/dx requirements for modes without phi determined per mode via Punzi FOM (see Table 2)")
  .note(:photon_selection, "Barrel: E>25MeV |cos(theta)|<0.80; Endcap: E>50MeV 0.86<|cos(theta)|<0.92; TDC [0,700]ns; angle to nearest track>10deg")
  .note(:ks0_reconstruction, "K_S0: two oppositely charged tracks |Vz|<20, vertex fit, mass [0.487,0.511], decay length L/sigma_L>2")
  .note(:pi0_selection, "pi0: M_gg in [0.115,0.150], 1C kinematic fit constraining to pi0 mass")
  .note(:signal_yield, "Signal yield from unbinned ML fit to M(h- h0 e+ e+) in [1.86, 2.04] GeV/c^2; signal shape: double-sided CB + bifurcated Gaussian; background shape: RooKeysPdf from inclusive MC")
  .note(:majorana_search, "Majorana neutrino nu_m searched in D_s+ -> phi e+ nu_m(->pi- e+) with m_nu_m in [0.20,0.80] GeV/c^2; M_pi-e+ within [m-5sigma, m+4sigma]; counting method with profile likelihood ULs")
  .with_decay_card(decay_card_phi_pi)
  .apply(sel_phi_pi)

alg_phi_pi.execute_on([ds_data_4180, ds_incMC_4180, exMC_phi_pi])

# ============================================================================
# Mode II: D_s+ -> phi K- e+ e+  (phi -> K+ K-)
# ============================================================================

decay_card_phi_K = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.000 phi K- e+ e+ PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

exMC_phi_K = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Ds_phi_K_ee"
  config.related_dataset = ds_data_4180
  config.events = 100_000
  config.decay_card = decay_card_phi_K
  config.cross_section = :default
end

alg_phi_K = Algorithm.new("DsToPhiKEE")
alg_phi_K.set_header(["DsToPhiKEEAlg/DsToPhiKEE.h"])
          .set_constant({ "ECMS" => [:double, 4.178] })

sel_phi_K = Selection.new
sel_phi_K.select_track {
            cos_theta 0.93
            Vz 10.0
            Vxy 1.0
            nChrp ">=3"
            nChrn ">=2"
            nNet "==1"
          }
          .select_photon {
            tdc_emc_start 0
            tdc_emc_end 700
            angle_to_track 10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion]
            identify :pion, against: [:kaon]
            identify :positron, against: [:kaon, :pion]
            nkp "==1"
            nkm "==2"
            nep "==2"
          }
          .kinematic_fit([:kp, :km, :ep]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_phi_K
  .note(:signal_mode, "D_s+ -> phi K- e+ e+ with phi -> K+ K-")
  .note(:production, "e+e- -> D_s*+- D_s-+ at c.m. energies 4.128-4.226 GeV; total 7.33 fb^-1")
  .note(:multi_energy, "Data at 8 energy points; representative MC at 4.178 GeV (see Mode I for full energy list)")
  .note(:e_pid, "At least one e+ must satisfy L_e/(L_e+L_pi+L_K) > 0.8 (looser than other modes); E/p > 0.70c (looser), chi2_dE/dx < 5.6 (optimized)")
  .note(:phi_selection, "phi: M_K+K- in [1.00, 1.05] GeV/c^2; best K+K- pair from two possible combinations selected by mass closest to nominal phi mass")
  .note(:background_veto, "See Mode I notes for full background veto details")
  .note(:signal_region, "M_rec/Delta_M requirement NOT applied for this mode (background already low)")
  .with_decay_card(decay_card_phi_K)
  .apply(sel_phi_K)

alg_phi_K.execute_on([ds_data_4180, ds_incMC_4180, exMC_phi_K])

# ============================================================================
# Mode III: D_s+ -> K_S0 pi- e+ e+  (K_S0 -> pi+ pi-)
# ============================================================================

decay_card_ks_pi = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.000 K_S0 pi- e+ e+ PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_ks_pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Ds_Ks_pi_ee"
  config.related_dataset = ds_data_4180
  config.events = 100_000
  config.decay_card = decay_card_ks_pi
  config.cross_section = :default
end

alg_ks_pi = Algorithm.new("DsToKsPiEE")
alg_ks_pi.set_header(["DsToKsPiEEAlg/DsToKsPiEE.h"])
          .set_constant({ "ECMS" => [:double, 4.178] })

sel_ks_pi = Selection.new
sel_ks_pi.select_track {
            cos_theta 0.93
            Vz 10.0
            Vxy 1.0
            nChrp ">=2"
            nChrn ">=2"
            nNet "==1"
          }
          .select_photon {
            tdc_emc_start 0
            tdc_emc_end 700
            angle_to_track 10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :pion, against: [:kaon]
            identify :positron, against: [:kaon, :pion]
            npip ">=1"
            npim ">=2"
            nep "==2"
          }
          .secondary_vertex_fit([:pip, :pim]) do
            build_virtual_particle(:K_S0).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          end
          .kinematic_fit([:K_S0, :pim, :ep]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_ks_pi
  .note(:signal_mode, "D_s+ -> K_S0 pi- e+ e+ with K_S0 -> pi+ pi-")
  .note(:production, "e+e- -> D_s*+- D_s-+ at c.m. energies 4.128-4.226 GeV; total 7.33 fb^-1")
  .note(:multi_energy, "Data at 8 energy points; representative MC at 4.178 GeV (see Mode I for full energy list)")
  .note(:e_pid, "Both e+ must satisfy L_e/(L_e+L_pi+L_K) > 0.8; E/p > 0.81c (optimized), chi2_dE/dx < 4.3 (optimized)")
  .note(:ks0_selection, "K_S0: |Vz|<20cm, vertex fit, M_pi+pi- in [0.487, 0.511] GeV/c^2, L/sigma_L > 2")
  .note(:background_veto, "cos(theta_e+e-) < 0.95; cos(theta_e+pi-) < 0.98; L/sigma_L > 2 (fake K_S0 veto)")
  .note(:signal_region, "M_rec and Delta_M signal region optimized per energy point (see Table 3)")
  .with_decay_card(decay_card_ks_pi)
  .apply(sel_ks_pi)

alg_ks_pi.execute_on([ds_data_4180, ds_incMC_4180, exMC_ks_pi])

# ============================================================================
# Mode IV: D_s+ -> K_S0 K- e+ e+  (K_S0 -> pi+ pi-)
# ============================================================================

decay_card_ks_K = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.000 K_S0 K- e+ e+ PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_ks_K = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Ds_Ks_K_ee"
  config.related_dataset = ds_data_4180
  config.events = 100_000
  config.decay_card = decay_card_ks_K
  config.cross_section = :default
end

alg_ks_K = Algorithm.new("DsToKsKEE")
alg_ks_K.set_header(["DsToKsKEEAlg/DsToKsKEE.h"])
          .set_constant({ "ECMS" => [:double, 4.178] })

sel_ks_K = Selection.new
sel_ks_K.select_track {
            cos_theta 0.93
            Vz 10.0
            Vxy 1.0
            nChrp ">=2"
            nChrn ">=2"
            nNet "==1"
          }
          .select_photon {
            tdc_emc_start 0
            tdc_emc_end 700
            angle_to_track 10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion]
            identify :pion, against: [:kaon]
            identify :positron, against: [:kaon, :pion]
            nkp ">=0"
            nkm ">=1"
            npip ">=1"
            npim ">=1"
            nep "==2"
          }
          .secondary_vertex_fit([:pip, :pim]) do
            build_virtual_particle(:K_S0).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          end
          .kinematic_fit([:K_S0, :km, :ep]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_ks_K
  .note(:signal_mode, "D_s+ -> K_S0 K- e+ e+ with K_S0 -> pi+ pi-")
  .note(:production, "e+e- -> D_s*+- D_s-+ at c.m. energies 4.128-4.226 GeV; total 7.33 fb^-1")
  .note(:multi_energy, "Data at 8 energy points; representative MC at 4.178 GeV (see Mode I for full energy list)")
  .note(:e_pid, "Both e+ must satisfy L_e/(L_e+L_pi+L_K) > 0.8; E/p > 0.70c (optimized), chi2_dE/dx < 5.6 (optimized)")
  .note(:ks0_selection, "K_S0: |Vz|<20cm, vertex fit, M_pi+pi- in [0.487, 0.511] GeV/c^2, L/sigma_L > 2")
  .note(:background_veto, "cos(theta_e+e-) < 0.95; L/sigma_L > 2 (fake K_S0 veto)")
  .note(:signal_region, "M_rec and Delta_M signal region optimized per energy point (see Table 3)")
  .with_decay_card(decay_card_ks_K)
  .apply(sel_ks_K)

alg_ks_K.execute_on([ds_data_4180, ds_incMC_4180, exMC_ks_K])

# ============================================================================
# Mode V: D_s+ -> pi- pi0 e+ e+  (pi0 -> gamma gamma)
# ============================================================================

decay_card_pi_pi0 = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.000 pi- pi0 e+ e+ PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_pi_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Ds_pi_pi0_ee"
  config.related_dataset = ds_data_4180
  config.events = 100_000
  config.decay_card = decay_card_pi_pi0
  config.cross_section = :default
end

alg_pi_pi0 = Algorithm.new("DsToPiPi0EE")
alg_pi_pi0.set_header(["DsToPiPi0EEAlg/DsToPiPi0EE.h"])
          .set_constant({ "ECMS" => [:double, 4.178] })

sel_pi_pi0 = Selection.new
sel_pi_pi0.select_track {
            cos_theta 0.93
            Vz 10.0
            Vxy 1.0
            nChrp ">=2"
            nChrn ">=1"
            nNet "==1"
          }
          .select_photon {
            tdc_emc_start 0
            tdc_emc_end 700
            angle_to_track 10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam ">=2"
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :pion, against: [:kaon]
            identify :positron, against: [:kaon, :pion]
            npim ">=1"
            nep "==2"
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 50
            npi0 ">=1"
          }
          .kinematic_fit([:pim, :pi0, :ep]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_pi_pi0
  .note(:signal_mode, "D_s+ -> pi- pi0 e+ e+ with pi0 -> gamma gamma")
  .note(:production, "e+e- -> D_s*+- D_s-+ at c.m. energies 4.128-4.226 GeV; total 7.33 fb^-1")
  .note(:multi_energy, "Data at 8 energy points; representative MC at 4.178 GeV (see Mode I for full energy list)")
  .note(:e_pid, "Both e+ must satisfy L_e/(L_e+L_pi+L_K) > 0.8; E/p > 0.81c (optimized), chi2_dE/dx < 4.3 (optimized)")
  .note(:pi0_selection, "pi0: M_gg in [0.115, 0.150] GeV/c^2, 1C kinematic fit constraining to pi0 mass; E(pi0) > 0.17 GeV (soft pi0 from D* veto)")
  .note(:background_veto, "cos(theta_e+e-) < 0.95; cos(theta_e+pi-) < 0.98; L_pi/sigma_L_pi < 3 (K_S0 veto)")
  .note(:signal_region, "M_rec and Delta_M signal region optimized per energy point (see Table 3)")
  .with_decay_card(decay_card_pi_pi0)
  .apply(sel_pi_pi0)

alg_pi_pi0.execute_on([ds_data_4180, ds_incMC_4180, exMC_pi_pi0])

# ============================================================================
# Mode VI: D_s+ -> K- pi0 e+ e+  (pi0 -> gamma gamma)
# ============================================================================

decay_card_K_pi0 = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.000 K- pi0 e+ e+ PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_K_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Ds_K_pi0_ee"
  config.related_dataset = ds_data_4180
  config.events = 100_000
  config.decay_card = decay_card_K_pi0
  config.cross_section = :default
end

alg_K_pi0 = Algorithm.new("DsToKPi0EE")
alg_K_pi0.set_header(["DsToKPi0EEAlg/DsToKPi0EE.h"])
          .set_constant({ "ECMS" => [:double, 4.178] })

sel_K_pi0 = Selection.new
sel_K_pi0.select_track {
            cos_theta 0.93
            Vz 10.0
            Vxy 1.0
            nChrp ">=2"
            nChrn ">=1"
            nNet "==1"
          }
          .select_photon {
            tdc_emc_start 0
            tdc_emc_end 700
            angle_to_track 10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam ">=2"
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion]
            identify :pion, against: [:kaon]
            identify :positron, against: [:kaon, :pion]
            nkm ">=1"
            nep "==2"
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 50
            npi0 ">=1"
          }
          .kinematic_fit([:km, :pi0, :ep]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
          }

alg_K_pi0
  .note(:signal_mode, "D_s+ -> K- pi0 e+ e+ with pi0 -> gamma gamma")
  .note(:production, "e+e- -> D_s*+- D_s-+ at c.m. energies 4.128-4.226 GeV; total 7.33 fb^-1")
  .note(:multi_energy, "Data at 8 energy points; representative MC at 4.178 GeV (see Mode I for full energy list)")
  .note(:e_pid, "Both e+ must satisfy L_e/(L_e+L_pi+L_K) > 0.8; E/p > 0.75c (optimized), chi2_dE/dx < 5.0 (optimized)")
  .note(:pi0_selection, "pi0: M_gg in [0.115, 0.150] GeV/c^2, 1C kinematic fit constraining to pi0 mass; E(pi0) > 0.17 GeV (soft pi0 veto)")
  .note(:background_veto, "cos(theta_e+e-) < 0.95; L_pi/sigma_L_pi < 3 (K_S0 veto); no cos(theta_e+pi-) veto for this mode (K- instead of pi-)")
  .note(:signal_region, "M_rec and Delta_M signal region optimized per energy point (see Table 3)")
  .with_decay_card(decay_card_K_pi0)
  .apply(sel_K_pi0)

alg_K_pi0.execute_on([ds_data_4180, ds_incMC_4180, exMC_K_pi0])