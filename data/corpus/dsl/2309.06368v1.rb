# BESIII Analysis: Measurements of absolute branching fractions of Omega- decays
# Paper: 2309.06368v1
# Data: psi(2S) at 3.686 GeV, BOSS 709

# ============================================================
# Datasets
# ============================================================
psip_data   = DatasetManager.real_data.find("709_3686")
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")

# ============================================================
# Signal Mode A: Omega- → Xi0 pi-, Xi0 → Lambda pi0
# (ST: Omega+ → Lambda K+)
# ============================================================
decay_card_Omega_Xi0pi = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Omega- Omega+  JPIPI;
  Enddecay
  Decay Omega-
  1.0000 Xi0 pi-  PHSP;
  Enddecay
  Decay Omega+
  1.0000 Lambda K+  PHSP;
  Enddecay
  Decay Xi0
  1.0000 Lambda pi0  PHSP;
  Enddecay
  Decay Lambda
  1.0000 p+ pi-  PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

signal_mc_A = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Omega_Xi0pi"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_Omega_Xi0pi
  config.cross_section   = :default
end

alg_A = Algorithm.new("Omega_Xi0pi_Analysis")
alg_A
  .set_header(["Omega_Xi0pi_Analysis/Omega_Xi0pi.h"])

sel_A = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=3"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
  end
  .remove([:prp <= :chrgp])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp ">=1"
  end
  .assign({chrgp: :pip, chrgn: :pim})
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=2"
  end
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  end
  .kinematic_fit([:Lambda, :kp, :pim, :Lambda, :pi0, :pim]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

alg_A
  .note(:st_omega_plus, "ST side: Omega+ reconstructed via Lambda K+ vertex fit and M(Lambda K+) mass window [1.664, 1.680] GeV")
  .note(:lambda_mass_window, "Lambda mass window |M(p pi-) - m_Lambda| < 2*sigma (~11 MeV)")
  .note(:highest_energy_pi, "Signal pi- selected as highest-energy pi- among remaining candidates")
  .note(:pi0_mass_window, "pi0 mass window [0.115, 0.150] GeV after 1C kinematic fit")
  .note(:recoil_mass_veto, "Various background vetoes using recoil mass of reconstructed particles")
  .note(:dt_yield_extraction, "DT yield extracted from M(Lambda K+) vs signal-side recoil mass 2D fit")
  .with_decay_card(decay_card_Omega_Xi0pi)
  .apply(sel_A)

alg_A.execute_on([psip_data, psip_incMC, signal_mc_A])

# ============================================================
# Signal Mode B: Omega- → Xi- pi0, Xi- → Lambda pi-
# (ST: Omega+ → Lambda K+)
# ============================================================
decay_card_Omega_Xipi0 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Omega- Omega+  JPIPI;
  Enddecay
  Decay Omega-
  1.0000 Xi- pi0  PHSP;
  Enddecay
  Decay Omega+
  1.0000 Lambda K+  PHSP;
  Enddecay
  Decay Xi-
  1.0000 Lambda pi-  PHSP;
  Enddecay
  Decay Lambda
  1.0000 p+ pi-  PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

signal_mc_B = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Omega_Xipi0"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_Omega_Xipi0
  config.cross_section   = :default
end

alg_B = Algorithm.new("Omega_Xipi0_Analysis")
alg_B
  .set_header(["Omega_Xipi0_Analysis/Omega_Xipi0.h"])

sel_B = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=3"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
  end
  .remove([:prp <= :chrgp])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp ">=1"
  end
  .assign({chrgp: :pip, chrgn: :pim})
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=2"
  end
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  end
  .kinematic_fit([:Lambda, :kp, :pim, :Lambda, :pi0]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

alg_B
  .note(:st_omega_plus, "ST side: Omega+ reconstructed via Lambda K+ vertex fit and mass window [1.664, 1.680] GeV")
  .note(:lambda_mass_window, "Lambda mass window |M(p pi-) - m_Lambda| < 11 MeV")
  .note(:highest_energy_pi, "Signal pi- selected as highest-energy pi-")
  .note(:pi0_mass_window, "pi0 mass window [0.115, 0.150] GeV after 1C kinematic fit")
  .note(:recoil_mass_veto, "Background vetoes using recoil mass")
  .note(:dt_yield_extraction, "DT yield from 2D mass fit")
  .with_decay_card(decay_card_Omega_Xipi0)
  .apply(sel_B)

alg_B.execute_on([psip_data, psip_incMC, signal_mc_B])

# ============================================================
# Signal Mode C: Omega- → Lambda K-
# (ST: Omega+ → Lambda K+)
# ============================================================
decay_card_Omega_LambdaK = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Omega- Omega+  JPIPI;
  Enddecay
  Decay Omega-
  1.0000 Lambda K-  PHSP;
  Enddecay
  Decay Omega+
  1.0000 Lambda K+  PHSP;
  Enddecay
  Decay Lambda
  1.0000 p+ pi-  PHSP;
  Enddecay
  End
DECAYCARD

signal_mc_C = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Omega_LambdaK"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_Omega_LambdaK
  config.cross_section   = :default
end

alg_C = Algorithm.new("Omega_LambdaK_Analysis")
alg_C
  .set_header(["Omega_LambdaK_Analysis/Omega_LambdaK.h"])

sel_C = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=3"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
  end
  .remove([:prp <= :chrgp])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp ">=1"
    nkm ">=1"
  end
  .assign({chrgp: :pip, chrgn: :pim})
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:Lambda, :kp, :Lambda, :km]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

alg_C
  .note(:st_omega_plus, "ST side: Omega+ reconstructed via Lambda K+ vertex fit and mass window [1.664, 1.680] GeV")
  .note(:lambda_mass_window, "Lambda mass window |M(p pi-) - m_Lambda| < 11 MeV")
  .note(:recoil_mass_veto, "Background vetoes using recoil mass of reconstructed Lambda K- system")
  .note(:dt_yield_extraction, "DT yield from 2D mass fit of M(Lambda K+) vs M(Lambda K-)")
  .note(:no_pi0, "Mode C has no pi0 in final state; no photon selection needed on signal side")
  .with_decay_card(decay_card_Omega_LambdaK)
  .apply(sel_C)

alg_C.execute_on([psip_data, psip_incMC, signal_mc_C])