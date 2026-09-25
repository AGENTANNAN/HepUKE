# Paper: 2307.09266v3
# Analysis: Lambda_c+ -> p eta and Lambda_c+ -> p omega branching fractions
# Energy: 7 points from 4.600 to 4.699 GeV
# Method: Ordinary analysis (single-tag method — direct reconstruction, not DTagAlg)

# --- Datasets (7 energy points) ---
data_4600 = DatasetManager.real_data.find("703_4600")
data_4612 = DatasetManager.real_data.find("706_4610")
data_4628 = DatasetManager.real_data.find("706_4620")
data_4641 = DatasetManager.real_data.find("706_4640")
data_4661 = DatasetManager.real_data.find("706_4660")
data_4682 = DatasetManager.real_data.find("706_4680")
data_4699 = DatasetManager.real_data.find("706_4700")

all_data = [data_4600, data_4612, data_4628, data_4641, data_4661, data_4682, data_4699]

incMC_4682 = DatasetManager.inclusive_mc.find("706_4680")

# ===================================================================
# Signal mode 1: Lambda_c+ -> p eta, eta -> gamma gamma
# ===================================================================

decay_card_p_eta_gg = <<~DECAYCARD
  Decay anti-Lambda_c-
  1.000 p+ eta PHSP;
  Enddecay
  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

sig_mc_p_eta_gg = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name = "sig_Lc_p_eta_gg"
  config.events = 100_000
  config.decay_card = decay_card_p_eta_gg
  config.cross_section = :default
end

algo_p_eta_gg = Algorithm.new("LcpEtaGG")
algo_p_eta_gg.set_header(["LcpEtaGGAlg/LcpEtaGG.h"])
  .set_constant({ "ECMS" => [:double, 4.600] })
  .with_decay_card(decay_card_p_eta_gg)
  .note(:analysis, "Lambda_c+ -> p eta, eta -> gamma gamma; ST method with direct reconstruction; M_BC vs Delta_E 2D simultaneous fit across 7 energy points")
  .note(:eta_reco, "eta -> gamma gamma: 1C kinematic fit constraining M_gg to nominal eta mass, chi2 < 20; |cos(theta_decay)| < 0.9 for helicity angle; lateral_moment < 0.4; E3x3/E5x5 > 0.85")
  .note(:proton_vr, "Additional proton Vr < 0.2 cm to reject beam-gas and beam-pipe backgrounds")
  .note(:signal_extraction, "2D unbinned ML fit on M_BC vs Delta_E; signal = MC shape convolved with 2D Gaussian; background = ARGUS x 2nd-order Chebyshev polynomial; BF shared across 7 energy points")
  .note(:bf_normalization, "BF = N_sig / (2 * N_LcLc * epsilon * BF_inter); N_LcLc = L_int * sigma from Ref.[53]")
  .note(:result, "BF(Lambda_c+ -> p eta) = (1.57 +/- 0.11_stat +/- 0.04_syst) x 10^{-3}")

sel_p_eta_gg = Selection.new
sel_p_eta_gg.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==1"
  nChrn "==0"
  nNet "==1"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 8.0
  nGam ">=2"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:pion, :kaon]
  nprp "==1"
end
.for_each(:prp) do
  where { vr < 0.2 }
end
.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 20
  neta 1
end
.kinematic_fit([:prp, :eta]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

algo_p_eta_gg.apply(sel_p_eta_gg)

# ===================================================================
# Signal mode 2: Lambda_c+ -> p eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma
# ===================================================================

decay_card_p_eta_3pi = <<~DECAYCARD
  Decay anti-Lambda_c-
  1.000 p+ eta PHSP;
  Enddecay
  Decay eta
  1.000 pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

sig_mc_p_eta_3pi = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name = "sig_Lc_p_eta_3pi"
  config.events = 100_000
  config.decay_card = decay_card_p_eta_3pi
  config.cross_section = :default
end

algo_p_eta_3pi = Algorithm.new("LcpEta3Pi")
algo_p_eta_3pi.set_header(["LcpEta3PiAlg/LcpEta3Pi.h"])
  .set_constant({ "ECMS" => [:double, 4.600] })
  .with_decay_card(decay_card_p_eta_3pi)
  .note(:analysis, "Lambda_c+ -> p eta, eta -> pi+ pi- pi0; M_BC vs Delta_E 2D simultaneous fit")
  .note(:eta_reco, "eta -> pi+ pi- pi0: invariant mass [0.536, 0.560] GeV/c^2; M(p pi0) veto [1.17, 1.20] GeV for Sigma+ suppression")
  .note(:pi0_reco, "pi0 -> gamma gamma: 1C kinematic fit chi2 < 50; M_gg [0.115, 0.150] GeV/c^2")
  .note(:proton_vr, "Additional proton Vr < 0.2 cm")
  .note(:bf_combined, "BF from p_eta_gg and p_eta_3pi consistent and combined")
  .note(:signal_extraction, "2D unbinned ML fit on M_BC vs Delta_E; signal = MC shape convolved with 2D Gaussian; background = ARGUS x Chebyshev")

sel_p_eta_3pi = Selection.new
sel_p_eta_3pi.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==1"
  nNet "==1"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 8.0
  nGam ">=2"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:pion, :kaon]
  identify :pion, against: [:proton, :kaon]
  nprp "==1"
  npip "==1"
  npim "==1"
end
.for_each(:prp) do
  where { vr < 0.2 }
end
.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 50
  npi0 1
end
.for_each(:pi0) do
  where { mass >= 0.115 && mass <= 0.150 }
end
.build_virtual_particle(:eta, from: [:pip, :pim, :pi0])
.for_each(:eta) do
  where { mass >= 0.536 && mass <= 0.560 }
end
.kinematic_fit([:prp, :eta]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

algo_p_eta_3pi.apply(sel_p_eta_3pi)

# ===================================================================
# Signal mode 3: Lambda_c+ -> p omega, omega -> pi+ pi- pi0, pi0 -> gamma gamma
# ===================================================================

decay_card_p_omega = <<~DECAYCARD
  Decay anti-Lambda_c-
  1.000 p+ omega PHSP;
  Enddecay
  Decay omega
  1.000 pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

sig_mc_p_omega = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name = "sig_Lc_p_omega"
  config.events = 100_000
  config.decay_card = decay_card_p_omega
  config.cross_section = :default
end

algo_p_omega = Algorithm.new("LcpOmega")
algo_p_omega.set_header(["LcpOmegaAlg/LcpOmega.h"])
  .set_constant({ "ECMS" => [:double, 4.600] })
  .with_decay_card(decay_card_p_omega)
  .note(:analysis, "Lambda_c+ -> p omega, omega -> pi+ pi- pi0; 1D M_BC fit with tight Delta_E window [-0.03, 0.02] GeV")
  .note(:omega_reco, "omega -> pi+ pi- pi0: invariant mass [0.750, 0.810] GeV/c^2; Dalitz R value < 0.9 to suppress non-omega background")
  .note(:vertex_fit, "Vertex fit applied to proton and two charged pions before omega reconstruction")
  .note(:vetoes, "Vetoes applied at ROOT level: M(p pi0) [1.17,1.20] (Sigma+), M(p pi-) [1.10,1.12] (Lambda), M(pi+pi-) [0.47,0.51] (K_S0)")
  .note(:signal_extraction, "1D unbinned ML simultaneous fit on M_BC for omega signal and sideband regions; signal = MC shape convolved with Gaussian; combinatorial background = ARGUS; non-omega peaking = Lambda_c+ -> p 3pi MC shape")
  .note(:result, "BF(Lambda_c+ -> p omega) = (1.11 +/- 0.20_stat +/- 0.07_syst) x 10^{-3}; significance 5.7 sigma")

sel_p_omega = Selection.new
sel_p_omega.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==1"
  nNet "==1"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 8.0
  nGam ">=2"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:pion, :kaon]
  identify :pion, against: [:proton, :kaon]
  nprp "==1"
  npip "==1"
  npim "==1"
end
.for_each(:prp) do
  where { vr < 0.2 }
end
.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 50
  npi0 1
end
.for_each(:pi0) do
  where { mass >= 0.115 && mass <= 0.150 }
end
.build_virtual_particle(:omega, from: [:pip, :pim, :pi0])
.for_each(:omega) do
  where { mass >= 0.750 && mass <= 0.810 }
end
.kinematic_fit([:prp, :omega]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

algo_p_omega.apply(sel_p_omega)

# ===================================================================
# Execute all three algorithms
# ===================================================================

algo_p_eta_gg.execute_on(all_data + [incMC_4682] + sig_mc_p_eta_gg)
algo_p_eta_3pi.execute_on(all_data + [incMC_4682] + sig_mc_p_eta_3pi)
algo_p_omega.execute_on(all_data + [incMC_4682] + sig_mc_p_omega)