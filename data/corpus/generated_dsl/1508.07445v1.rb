# =====================================================================
# Datasets — ψ(3770), ψ(3686), J/ψ and the 3.65 GeV continuum
# =====================================================================
data_3773 = DatasetManager.real_data.find("712_3773")   # ψ(3770), 2.92 fb⁻¹
data_3686 = DatasetManager.real_data.find("709_3686")   # ψ(3686), 106 M
data_3097 = DatasetManager.real_data.find("708_3097")   # J/ψ, 1.31 B
data_3650 = DatasetManager.real_data.find("709_3650")   # 3.65 GeV continuum, 42 pb⁻¹ (QED bkg)

incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")
incMC_3097 = DatasetManager.inclusive_mc.find("708_3097")
incMC_3650 = DatasetManager.inclusive_mc.find("709_3650")

# =====================================================================
# Decay cards — ψ(3770): seven baryonic channels
# =====================================================================
dc_LLpipi = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Lambda0 anti-Lambda0 pi+ pi- PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay
  End
DECAYCARD

dc_LLpi0 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Lambda0 anti-Lambda0 pi0 PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

dc_LLeta = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Lambda0 anti-Lambda0 eta PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

dc_SpSm = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Sigma+ anti-Sigma- PHSP;
  Enddecay
  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay
  Decay anti-Sigma-
  1.0000 anti-p- pi0 PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

dc_S0S0bar = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Sigma0 anti-Sigma0 PHSP;
  Enddecay
  Decay Sigma0
  1.0000 Lambda0 gamma PHSP;
  Enddecay
  Decay anti-Sigma0
  1.0000 anti-Lambda0 gamma PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay
  End
DECAYCARD

dc_XmXp = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Xi- anti-Xi+ PHSP;
  Enddecay
  Decay Xi-
  1.0000 Lambda0 pi- PHSP;
  Enddecay
  Decay anti-Xi+
  1.0000 anti-Lambda0 pi+ PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay
  End
DECAYCARD

dc_X0X0bar = <<~DECAYCARD
  Decay psi(3770)
  1.0000 Xi0 anti-Xi0 PHSP;
  Enddecay
  Decay Xi0
  1.0000 Lambda0 pi0 PHSP;
  Enddecay
  Decay anti-Xi0
  1.0000 anti-Lambda0 pi0 PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# ---- ψ(3770): radiative η_c(2S) and χ_cJ → γ J/ψ ----
dc_etac2s = <<~DECAYCARD
  Decay psi(3770)
  1.0000 gamma eta_c(2S) PHSP;
  Enddecay
  Decay eta_c(2S)
  1.0000 K_S0 K+ pi- PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

dc_chic1_ee = <<~DECAYCARD
  Decay psi(3770)
  1.0000 gamma chi_c1 PHSP;
  Enddecay
  Decay chi_c1
  1.0000 gamma J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay
  End
DECAYCARD

dc_chic1_mumu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 gamma chi_c1 PHSP;
  Enddecay
  Decay chi_c1
  1.0000 gamma J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0000 mu+ mu- PHOTOS VLL;
  Enddecay
  End
DECAYCARD

dc_chic2_ee = <<~DECAYCARD
  Decay psi(3770)
  1.0000 gamma chi_c2 PHSP;
  Enddecay
  Decay chi_c2
  1.0000 gamma J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay
  End
DECAYCARD

dc_chic2_mumu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 gamma chi_c2 PHSP;
  Enddecay
  Decay chi_c2
  1.0000 gamma J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0000 mu+ mu- PHOTOS VLL;
  Enddecay
  End
DECAYCARD

# ---- ψ(3686): γχ_c0,2 → π⁰ η_c ; π⁺π⁻J/ψ ----
dc_chic0_pi0etac = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c0 PHSP;
  Enddecay
  Decay chi_c0
  1.0000 pi0 eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 K_S0 K+ pi- PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

dc_chic2_pi0etac = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c2 PHSP;
  Enddecay
  Decay chi_c2
  1.0000 pi0 eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 K_S0 K+ pi- PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

dc_pipiJpsi_gg = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi+ pi- J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

dc_pipiJpsi_gphi = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi+ pi- J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0000 gamma phi PHSP;
  Enddecay
  Decay phi
  1.0000 K+ K- VSS;
  Enddecay
  End
DECAYCARD

# ---- J/ψ → π⁰ φ ----
dc_pi0phi = <<~DECAYCARD
  Decay J/psi
  1.0000 pi0 phi PHSP;
  Enddecay
  Decay phi
  1.0000 K+ K- VSS;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# =====================================================================
# Exclusive MC — 200k per mode, 500k for J/ψ → π⁰φ
# =====================================================================
exMC_LLpipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_LLpipi"
  config.related_dataset = data_3773
  config.events = 200_000
  config.decay_card = dc_LLpipi
  config.cross_section = :default
end

exMC_LLpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_LLpi0"
  config.related_dataset = data_3773
  config.events = 200_000
  config.decay_card = dc_LLpi0
  config.cross_section = :default
end

exMC_LLeta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_LLeta"
  config.related_dataset = data_3773
  config.events = 200_000
  config.decay_card = dc_LLeta
  config.cross_section = :default
end

exMC_SpSm = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_SpSm"
  config.related_dataset = data_3773
  config.events = 200_000
  config.decay_card = dc_SpSm
  config.cross_section = :default
end

exMC_S0S0bar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_S0S0bar"
  config.related_dataset = data_3773
  config.events = 200_000
  config.decay_card = dc_S0S0bar
  config.cross_section = :default
end

exMC_XmXp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_XmXp"
  config.related_dataset = data_3773
  config.events = 200_000
  config.decay_card = dc_XmXp
  config.cross_section = :default
end

exMC_X0X0bar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_X0X0bar"
  config.related_dataset = data_3773
  config.events = 200_000
  config.decay_card = dc_X0X0bar
  config.cross_section = :default
end

exMC_etac2s = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_etac2s"
  config.related_dataset = data_3773
  config.events = 200_000
  config.decay_card = dc_etac2s
  config.cross_section = :default
end

exMC_chic1_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_chic1_ee"
  config.related_dataset = data_3773
  config.events = 200_000
  config.decay_card = dc_chic1_ee
  config.cross_section = :default
end

exMC_chic1_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_chic1_mumu"
  config.related_dataset = data_3773
  config.events = 200_000
  config.decay_card = dc_chic1_mumu
  config.cross_section = :default
end

exMC_chic2_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_chic2_ee"
  config.related_dataset = data_3773
  config.events = 200_000
  config.decay_card = dc_chic2_ee
  config.cross_section = :default
end

exMC_chic2_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_chic2_mumu"
  config.related_dataset = data_3773
  config.events = 200_000
  config.decay_card = dc_chic2_mumu
  config.cross_section = :default
end

exMC_chic0_pi0etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3686_chic0_pi0etac"
  config.related_dataset = data_3686
  config.events = 200_000
  config.decay_card = dc_chic0_pi0etac
  config.cross_section = :default
end

exMC_chic2_pi0etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3686_chic2_pi0etac"
  config.related_dataset = data_3686
  config.events = 200_000
  config.decay_card = dc_chic2_pi0etac
  config.cross_section = :default
end

exMC_pipiJpsi_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3686_pipiJpsi_gg"
  config.related_dataset = data_3686
  config.events = 200_000
  config.decay_card = dc_pipiJpsi_gg
  config.cross_section = :default
end

exMC_pipiJpsi_gphi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3686_pipiJpsi_gphi"
  config.related_dataset = data_3686
  config.events = 200_000
  config.decay_card = dc_pipiJpsi_gphi
  config.cross_section = :default
end

exMC_pi0phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3097_pi0phi"
  config.related_dataset = data_3097
  config.events = 500_000
  config.decay_card = dc_pi0phi
  config.cross_section = :default
end

# =====================================================================
# ψ(3770) — baryonic channels
# =====================================================================

# --- ψ(3770) → Λ Λ̄ π⁺ π⁻ ---
alg_LLpipi = Algorithm.new("LLpipi")
alg_LLpipi.set_header(["LLpipiAlg/LLpipi.h"]).set_constant({"ECMS" => [:double, 3.773]})
sel_LLpipi = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==3"; nChrn "==3"; nNet "==0" }
  .pid(method: :probability) { prob_cut 0.001; identify :proton, against: [:kaon, :pion]; nprp "==1"; nprm "==1" }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:prp, :pim]) { build_virtual_particle(:Lambda).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .secondary_vertex_fit([:prm, :pip]) { build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim]) { nominal; constrain_four_momentum; chi2_cut 200 }
alg_LLpipi.with_decay_card(dc_LLpipi).apply(sel_LLpipi)

# --- ψ(3770) → Λ Λ̄ π⁰ ---
alg_LLpi0 = Algorithm.new("LLpi0")
alg_LLpi0.set_header(["LLpi0Alg/LLpi0.h"]).set_constant({"ECMS" => [:double, 3.773]})
sel_LLpi0 = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==2"; nChrn "==2"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=2" }
  .pid(method: :probability) { prob_cut 0.001; identify :proton, against: [:kaon, :pion]; nprp "==1"; nprm "==1" }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:prp, :pim]) { build_virtual_particle(:Lambda).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .secondary_vertex_fit([:prm, :pip]) { build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .kalman_kinematic_fit([:gamma, :gamma]) { invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0); chi2_cut 25; npi0 ">=1" }
  .kinematic_fit([:Lambda, :Lambda_bar, :pi0]) { nominal; constrain_four_momentum; chi2_cut 200 }
alg_LLpi0.with_decay_card(dc_LLpi0).apply(sel_LLpi0)

# --- ψ(3770) → Λ Λ̄ η ---
alg_LLeta = Algorithm.new("LLeta")
alg_LLeta.set_header(["LLetaAlg/LLeta.h"]).set_constant({"ECMS" => [:double, 3.773]})
sel_LLeta = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==2"; nChrn "==2"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=2" }
  .pid(method: :probability) { prob_cut 0.001; identify :proton, against: [:kaon, :pion]; nprp "==1"; nprm "==1" }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:prp, :pim]) { build_virtual_particle(:Lambda).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .secondary_vertex_fit([:prm, :pip]) { build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .kalman_kinematic_fit([:gamma, :gamma]) { invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta); chi2_cut 25; neta ">=1" }
  .kinematic_fit([:Lambda, :Lambda_bar, :eta]) { nominal; constrain_four_momentum; chi2_cut 200 }
alg_LLeta.with_decay_card(dc_LLeta).apply(sel_LLeta)

# --- ψ(3770) → Σ⁺ Σ̄⁻ ---
alg_SpSm = Algorithm.new("SpSm")
alg_SpSm.set_header(["SpSmAlg/SpSm.h"]).set_constant({"ECMS" => [:double, 3.773]})
sel_SpSm = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==1"; nChrn "==1"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=4" }
  .pid(method: :probability) { prob_cut 0.001; identify :proton, against: [:kaon, :pion]; nprp "==1"; nprm "==1" }
  .kalman_kinematic_fit([:gamma, :gamma]) { invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0); chi2_cut 25; npi0 ">=2" }
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) { nominal; constrain_four_momentum; chi2_cut 200 }
alg_SpSm.with_decay_card(dc_SpSm).apply(sel_SpSm)

# --- ψ(3770) → Σ⁰ Σ̄⁰ ---
alg_S0S0bar = Algorithm.new("S0S0bar")
alg_S0S0bar.set_header(["S0S0barAlg/S0S0bar.h"]).set_constant({"ECMS" => [:double, 3.773]})
sel_S0S0bar = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==2"; nChrn "==2"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=2" }
  .pid(method: :probability) { prob_cut 0.001; identify :proton, against: [:kaon, :pion]; nprp "==1"; nprm "==1" }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:prp, :pim]) { build_virtual_particle(:Lambda).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .secondary_vertex_fit([:prm, :pip]) { build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) { nominal; constrain_four_momentum; chi2_cut 200 }
alg_S0S0bar.with_decay_card(dc_S0S0bar).apply(sel_S0S0bar)

# --- ψ(3770) → Ξ⁻ Ξ̄⁺ ---
alg_XmXp = Algorithm.new("XmXp")
alg_XmXp.set_header(["XmXpAlg/XmXp.h"]).set_constant({"ECMS" => [:double, 3.773]})
sel_XmXp = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==3"; nChrn "==3"; nNet "==0" }
  .pid(method: :probability) { prob_cut 0.001; identify :proton, against: [:kaon, :pion]; nprp "==1"; nprm "==1" }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:prp, :pim]) { build_virtual_particle(:Lambda).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .secondary_vertex_fit([:prm, :pip]) { build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .kinematic_fit([:Lambda, :Lambda_bar, :pim, :pip]) { nominal; constrain_four_momentum; chi2_cut 200 }
alg_XmXp.with_decay_card(dc_XmXp).apply(sel_XmXp)

# --- ψ(3770) → Ξ⁰ Ξ̄⁰ ---
alg_X0X0bar = Algorithm.new("X0X0bar")
alg_X0X0bar.set_header(["X0X0barAlg/X0X0bar.h"]).set_constant({"ECMS" => [:double, 3.773]})
sel_X0X0bar = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==2"; nChrn "==2"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=4" }
  .pid(method: :probability) { prob_cut 0.001; identify :proton, against: [:kaon, :pion]; nprp "==1"; nprm "==1" }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:prp, :pim]) { build_virtual_particle(:Lambda).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .secondary_vertex_fit([:prm, :pip]) { build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .kalman_kinematic_fit([:gamma, :gamma]) { invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0); chi2_cut 25; npi0 ">=2" }
  .kinematic_fit([:Lambda, :Lambda_bar, :pi0, :pi0]) { nominal; constrain_four_momentum; chi2_cut 200 }
alg_X0X0bar.with_decay_card(dc_X0X0bar).apply(sel_X0X0bar)

# --- ψ(3770) → γ η_c(2S), η_c(2S) → K_S⁰ K⁺ π⁻ ---
alg_etac2s = Algorithm.new("EtaC2S")
alg_etac2s.set_header(["EtaC2SAlg/EtaC2S.h"]).set_constant({"ECMS" => [:double, 3.773]})
sel_etac2s = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==2"; nChrn "==2"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=1" }
  .pid(method: :probability) { prob_cut 0.001; identify :kaon, against: [:pion, :proton]; identify :pion, against: [:kaon, :proton]; nkp "==1" }
  .secondary_vertex_fit([:pip, :pim]) { build_virtual_particle(:K_S0).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .kinematic_fit([:gamma, :K_S0, :kp, :pim]) { nominal; constrain_four_momentum; chi2_cut 200 }
alg_etac2s.with_decay_card(dc_etac2s).apply(sel_etac2s)

# --- ψ(3770) → γ χ_cJ → γ γ ℓ⁺ℓ⁻ (e, μ for χ_c1 and χ_c2) ---
alg_chic1_ee = Algorithm.new("ChiC1GamEE")
alg_chic1_ee.set_header(["ChiC1GamEEAlg/ChiC1GamEE.h"]).set_constant({"ECMS" => [:double, 3.773]})
sel_chic1_ee = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==1"; nChrn "==1"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=2" }
  .pid(method: :probability) { prob_cut 0.001; identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0, treat_as_electron_if_energy_above: 0.6; nlp "==1"; nlm "==1" }
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) { nominal; constrain_four_momentum; invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi); chi2_cut 200 }
alg_chic1_ee.with_decay_card(dc_chic1_ee).apply(sel_chic1_ee)

alg_chic1_mumu = Algorithm.new("ChiC1GamMuMu")
alg_chic1_mumu.set_header(["ChiC1GamMuMuAlg/ChiC1GamMuMu.h"]).set_constant({"ECMS" => [:double, 3.773]})
sel_chic1_mumu = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==1"; nChrn "==1"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=2" }
  .pid(method: :probability) { prob_cut 0.001; identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0, treat_as_electron_if_energy_above: 0.6; nlp "==1"; nlm "==1" }
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) { nominal; constrain_four_momentum; invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi); chi2_cut 200 }
alg_chic1_mumu.with_decay_card(dc_chic1_mumu).apply(sel_chic1_mumu)

alg_chic2_ee = Algorithm.new("ChiC2GamEE")
alg_chic2_ee.set_header(["ChiC2GamEEAlg/ChiC2GamEE.h"]).set_constant({"ECMS" => [:double, 3.773]})
sel_chic2_ee = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==1"; nChrn "==1"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=2" }
  .pid(method: :probability) { prob_cut 0.001; identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0, treat_as_electron_if_energy_above: 0.6; nlp "==1"; nlm "==1" }
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) { nominal; constrain_four_momentum; invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi); chi2_cut 200 }
alg_chic2_ee.with_decay_card(dc_chic2_ee).apply(sel_chic2_ee)

alg_chic2_mumu = Algorithm.new("ChiC2GamMuMu")
alg_chic2_mumu.set_header(["ChiC2GamMuMuAlg/ChiC2GamMuMu.h"]).set_constant({"ECMS" => [:double, 3.773]})
sel_chic2_mumu = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==1"; nChrn "==1"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=2" }
  .pid(method: :probability) { prob_cut 0.001; identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0, treat_as_electron_if_energy_above: 0.6; nlp "==1"; nlm "==1" }
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) { nominal; constrain_four_momentum; invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi); chi2_cut 200 }
alg_chic2_mumu.with_decay_card(dc_chic2_mumu).apply(sel_chic2_mumu)

# =====================================================================
# ψ(3686) — γχ_cJ → π⁰ η_c  and  π⁺π⁻ J/ψ
# =====================================================================

# --- ψ(3686) → γ χ_c0, χ_c0 → π⁰ η_c, η_c → K_S⁰ K⁺ π⁻ ---
alg_chic0_pi0etac = Algorithm.new("ChiC0Pi0EtaC")
alg_chic0_pi0etac.set_header(["ChiC0Pi0EtaCAlg/ChiC0Pi0EtaC.h"]).set_constant({"ECMS" => [:double, 3.686]})
sel_chic0_pi0etac = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==2"; nChrn "==2"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=3" }
  .pid(method: :probability) { prob_cut 0.001; identify :kaon, against: [:pion, :proton]; identify :pion, against: [:kaon, :proton]; nkp "==1" }
  .kalman_kinematic_fit([:gamma, :gamma]) { invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0); chi2_cut 25; npi0 ">=1" }
  .secondary_vertex_fit([:pip, :pim]) { build_virtual_particle(:K_S0).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .kinematic_fit([:gamma, :pi0, :K_S0, :kp, :pim]) { nominal; constrain_four_momentum; chi2_cut 200 }
alg_chic0_pi0etac.with_decay_card(dc_chic0_pi0etac).apply(sel_chic0_pi0etac)

# --- ψ(3686) → γ χ_c2, χ_c2 → π⁰ η_c, η_c → K_S⁰ K⁺ π⁻ ---
alg_chic2_pi0etac = Algorithm.new("ChiC2Pi0EtaC")
alg_chic2_pi0etac.set_header(["ChiC2Pi0EtaCAlg/ChiC2Pi0EtaC.h"]).set_constant({"ECMS" => [:double, 3.686]})
sel_chic2_pi0etac = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==2"; nChrn "==2"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=3" }
  .pid(method: :probability) { prob_cut 0.001; identify :kaon, against: [:pion, :proton]; identify :pion, against: [:kaon, :proton]; nkp "==1" }
  .kalman_kinematic_fit([:gamma, :gamma]) { invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0); chi2_cut 25; npi0 ">=1" }
  .secondary_vertex_fit([:pip, :pim]) { build_virtual_particle(:K_S0).by_minimizing_mass_difference; remove_used_particle_from_candidate_list }
  .kinematic_fit([:gamma, :pi0, :K_S0, :kp, :pim]) { nominal; constrain_four_momentum; chi2_cut 200 }
alg_chic2_pi0etac.with_decay_card(dc_chic2_pi0etac).apply(sel_chic2_pi0etac)

# --- ψ(3686) → π⁺π⁻ J/ψ, J/ψ → γ γ ---
alg_pipiJpsi_gg = Algorithm.new("PipiJpsiGG")
alg_pipiJpsi_gg.set_header(["PipiJpsiGGAlg/PipiJpsiGG.h"]).set_constant({"ECMS" => [:double, 3.686]})
sel_pipiJpsi_gg = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==1"; nChrn "==1"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=2" }
  .pid(method: :probability) { prob_cut 0.001; identify :pion, against: [:kaon, :proton]; identify :kaon, against: [:pion, :proton]; npip "==1"; npim "==1" }
  .kinematic_fit([:pip, :pim, :gamma, :gamma]) { nominal; constrain_four_momentum; chi2_cut 200 }
alg_pipiJpsi_gg.with_decay_card(dc_pipiJpsi_gg).apply(sel_pipiJpsi_gg)

# --- ψ(3686) → π⁺π⁻ J/ψ, J/ψ → γ φ, φ → K⁺K⁻ ---
alg_pipiJpsi_gphi = Algorithm.new("PipiJpsiGPhi")
alg_pipiJpsi_gphi.set_header(["PipiJpsiGPhiAlg/PipiJpsiGPhi.h"]).set_constant({"ECMS" => [:double, 3.686]})
sel_pipiJpsi_gphi = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==2"; nChrn "==2"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=1" }
  .pid(method: :probability) { prob_cut 0.001; identify :pion, against: [:kaon, :proton]; identify :kaon, against: [:pion, :proton]; npip "==1"; npim "==1"; nkp "==1"; nkm "==1" }
  .kinematic_fit([:pip, :pim, :gamma, :kp, :km]) { nominal; constrain_four_momentum; chi2_cut 200 }
alg_pipiJpsi_gphi.with_decay_card(dc_pipiJpsi_gphi).apply(sel_pipiJpsi_gphi)

# =====================================================================
# J/ψ → π⁰ φ, φ → K⁺K⁻
# =====================================================================
alg_pi0phi = Algorithm.new("Pi0Phi")
alg_pi0phi.set_header(["Pi0PhiAlg/Pi0Phi.h"]).set_constant({"ECMS" => [:double, 3.097]})
sel_pi0phi = Selection.new
  .select_track { cos_theta 0.93; Vz 10.0; Vr 1.0; nChrp "==1"; nChrn "==1"; nNet "==0" }
  .select_photon { tdc_emc_start 0; tdc_emc_end 14; angle_to_track 20.0; energyThreshold_b 0.025; energyThreshold_e 0.050; nGam ">=2" }
  .pid(method: :probability) { prob_cut 0.001; identify :kaon, against: [:pion, :proton]; nkp "==1"; nkm "==1" }
  .kalman_kinematic_fit([:gamma, :gamma]) { invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0); chi2_cut 25; npi0 ">=1" }
  .kinematic_fit([:pi0, :kp, :km]) { nominal; constrain_four_momentum; chi2_cut 200 }
alg_pi0phi.with_decay_card(dc_pi0phi).apply(sel_pi0phi)

# =====================================================================
# Execute on real data + inclusive MC + 3.65 GeV continuum + exclusive MC
# =====================================================================
root_files_LLpipi = alg_LLpipi.execute_on([data_3773, incMC_3773, data_3650, incMC_3650, exMC_LLpipi])
root_files_LLpi0 = alg_LLpi0.execute_on([data_3773, incMC_3773, data_3650, incMC_3650, exMC_LLpi0])
root_files_LLeta = alg_LLeta.execute_on([data_3773, incMC_3773, data_3650, incMC_3650, exMC_LLeta])
root_files_SpSm = alg_SpSm.execute_on([data_3773, incMC_3773, data_3650, incMC_3650, exMC_SpSm])
root_files_S0S0bar = alg_S0S0bar.execute_on([data_3773, incMC_3773, data_3650, incMC_3650, exMC_S0S0bar])
root_files_XmXp = alg_XmXp.execute_on([data_3773, incMC_3773, data_3650, incMC_3650, exMC_XmXp])
root_files_X0X0bar = alg_X0X0bar.execute_on([data_3773, incMC_3773, data_3650, incMC_3650, exMC_X0X0bar])
root_files_etac2s = alg_etac2s.execute_on([data_3773, incMC_3773, data_3650, incMC_3650, exMC_etac2s])
root_files_chic1_ee = alg_chic1_ee.execute_on([data_3773, incMC_3773, data_3650, incMC_3650, exMC_chic1_ee])
root_files_chic1_mumu = alg_chic1_mumu.execute_on([data_3773, incMC_3773, data_3650, incMC_3650, exMC_chic1_mumu])
root_files_chic2_ee = alg_chic2_ee.execute_on([data_3773, incMC_3773, data_3650, incMC_3650, exMC_chic2_ee])
root_files_chic2_mumu = alg_chic2_mumu.execute_on([data_3773, incMC_3773, data_3650, incMC_3650, exMC_chic2_mumu])
root_files_chic0_pi0etac = alg_chic0_pi0etac.execute_on([data_3686, incMC_3686, data_3650, incMC_3650, exMC_chic0_pi0etac])
root_files_chic2_pi0etac = alg_chic2_pi0etac.execute_on([data_3686, incMC_3686, data_3650, incMC_3650, exMC_chic2_pi0etac])
root_files_pipiJpsi_gg = alg_pipiJpsi_gg.execute_on([data_3686, incMC_3686, data_3650, incMC_3650, exMC_pipiJpsi_gg])
root_files_pipiJpsi_gphi = alg_pipiJpsi_gphi.execute_on([data_3686, incMC_3686, data_3650, incMC_3650, exMC_pipiJpsi_gphi])
root_files_pi0phi = alg_pi0phi.execute_on([data_3097, incMC_3097, data_3650, incMC_3650, exMC_pi0phi])