# DSL for 2212.13361v2: D_s^{*+} → D_s^+ π^0 / D_s^{*+} → D_s^+ γ BF ratio
# Ordinary analysis: D_s pairs, 3 decay modes, 2C kinematic fit
# M_miss^2 separation between D_sγ and D_sπ^0 channels
# 8 energy points: 4.128-4.226 GeV, 7.33 fb^{-1} total

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# All 8 energy points (in order of increasing E_cm)
data_4130 = DatasetManager.real_data.find("705_4130")
data_4160 = DatasetManager.real_data.find("705_4160")
data_4180 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")

incMC_4130 = DatasetManager.inclusive_mc.find("705_4130")
incMC_4160 = DatasetManager.inclusive_mc.find("705_4160")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

all_data = [data_4130, data_4160, data_4180, data_4190, data_4200, data_4210, data_4220, data_4230]
all_incMC = [incMC_4130, incMC_4160, incMC_4180, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4230]

energy_points = [4.130, 4.160, 4.178, 4.190, 4.200, 4.210, 4.220, 4.230]

# ============================================================
# Mode I: D_s^+ → K^+ K^- π^+  vs  D_s^- → K^+ K^- π^-
# ============================================================

decay_card_mode1 = <<~DECAY
Decay e+ e-
  1.0  D_s*+  D_s*-  VSS;
Enddecay
Decay D_s*+
  1.0  D_s+  gamma  PHSP;
Enddecay
Decay D_s*-
  1.0  D_s-  pi0  PHSP;
Enddecay
Decay D_s+
  1.0  K+  K-  pi+  PHSP;
Enddecay
Decay D_s-
  1.0  K+  K-  pi-  PHSP;
Enddecay
Decay pi0
  1.0  gamma  gamma  PHSP;
Enddecay
DECAY

sigMC_mode1 = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name     = "sig_DsStar_DsKKpi_DsKKpi"
  config.events          = 500_000
  config.decay_card      = decay_card_mode1
  config.cross_section   = :default
end

alg_mode1 = Algorithm.new("DsStarDsKKpiModeI", "00-00-01")
alg_mode1.set_header(["DsStarDsKKpiAlg/DsStarDsKKpi.h"])
          .set_constant({ "ECMS" => [:double, 4.178] })
          .with_decay_card(decay_card_mode1)

sel_mode1 = Selection.new

sel_mode1.select_track do |t|
  t.nChrp 2
  t.nChrn 2
  t.nTot 4
end

sel_mode1.pid(method: :probability) do |pid|
  pid.identify :kp, :km, against: [:pip, :ep]
  pid.identify :pip, :pim, against: [:kp, :ep]
  pid.prob_cut 0.001
end

sel_mode1.select_photon do |p|
  p.min_energy 0.025
end

sel_mode1.build_virtual_particle(:DsP_cand, from: [:kp, :km, :pip])
sel_mode1.build_virtual_particle(:DsM_cand, from: [:kp, :km, :pim])

sel_mode1.kinematic_fit([:kp, :km, :pip, :kp, :km, :pim]) do |fit|
  fit.constrain_four_momentum
  fit.invariant_mass_of(:kp, :km, :pip).constrain_to_nominal_mass_of(:Ds)
  fit.invariant_mass_of(:kp, :km, :pim).constrain_to_nominal_mass_of(:Ds)
  fit.chi2_cut 200
  fit.nominal
end

alg_mode1.note(:mmiss2, "M_miss^2 = (E_cm - E_Ds+ - E_Ds-)^2 - |p_Ds+ + p_Ds-|^2 used to separate D_s^* → D_sγ (M_miss^2 near 0) from D_s^* → D_sπ^0 (M_miss^2 near M_π0^2)")
alg_mode1.note(:kinematic_fit, "2C kinematic fit constraining both D_s masses to nominal. Fit performed with all photon candidates 0..N")
alg_mode1.note(:photon_pair, "D_s^* → D_sπ^0 channel: recoil π^0 reconstructed from γγ pair. D_s^* → D_sγ channel: single photon from transition")
alg_mode1.note(:ratio, "BF ratio R = B(D_s^*+→D_s+π0) / B(D_s^*+→D_s+γ) extracted from simultaneous fit to M_miss^2 at all 8 energies")
alg_mode1.note(:cross_feed, "Cross-feed between γ and π^0 channels modeled with MC-derived transfer factors")
alg_mode1.note(:isospin, "Isospin symmetry: B(D_s^*+→D_s+π^0) = B(D_s^*+→D_s0π^+). Combined signal yields from both channels")
alg_mode1.note(:background, "Background from non-D_sD_s events estimated from D_s sidebands in M_BC and ΔE")
alg_mode1.note(:energy_scan, "Data taken at 8 c.m. energies: 4.130, 4.160, 4.178, 4.190, 4.200, 4.210, 4.220, 4.230 GeV. BOSS versions: 705 for 4.130/4.160, 703 for 4.178-4.230")

alg_mode1.apply(sel_mode1)
alg_mode1.execute_on(all_data + all_incMC + sigMC_mode1.to_a)

# ============================================================
# Mode II: D_s^+ → K^+ K^- π^+  vs  D_s^- → K_S^0 K^-
# ============================================================

decay_card_mode2 = <<~DECAY
Decay e+ e-
  1.0  D_s*+  D_s*-  VSS;
Enddecay
Decay D_s*+
  1.0  D_s+  gamma  PHSP;
Enddecay
Decay D_s*-
  1.0  D_s-  pi0  PHSP;
Enddecay
Decay D_s+
  1.0  K+  K-  pi+  PHSP;
Enddecay
Decay D_s-
  1.0  K_S0  K-  PHSP;
Enddecay
Decay K_S0
  1.0  pi+  pi-  PHSP;
Enddecay
Decay pi0
  1.0  gamma  gamma  PHSP;
Enddecay
DECAY

sigMC_mode2 = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name     = "sig_DsStar_DsKKpi_DsKsK"
  config.events          = 500_000
  config.decay_card      = decay_card_mode2
  config.cross_section   = :default
end

alg_mode2 = Algorithm.new("DsStarDsKKpiModeII", "00-00-01")
alg_mode2.set_header(["DsStarDsKKpiAlg/DsStarDsKKpi.h"])
          .set_constant({ "ECMS" => [:double, 4.178] })
          .with_decay_card(decay_card_mode2)

sel_mode2 = Selection.new

sel_mode2.select_track do |t|
  t.nChrp 2
  t.nChrn 3
  t.nTot 5
end

sel_mode2.pid(method: :probability) do |pid|
  pid.identify :kp, :km, against: [:pip, :ep]
  pid.identify :pip, :pim, against: [:kp, :ep]
  pid.prob_cut 0.001
end

sel_mode2.select_photon do |p|
  p.min_energy 0.025
end

sel_mode2.secondary_vertex_fit(:K_S0, daughters: [:pip, :pim]) do |v|
  v.mass_window [0.485, 0.510]
  v.flight_significance 2.0
end

sel_mode2.build_virtual_particle(:DsP_cand, from: [:kp, :km, :pip])
sel_mode2.build_virtual_particle(:DsM_cand, from: [:K_S0, :km])

sel_mode2.kinematic_fit([:kp, :km, :pip, :K_S0, :km]) do |fit|
  fit.constrain_four_momentum
  fit.invariant_mass_of(:kp, :km, :pip).constrain_to_nominal_mass_of(:Ds)
  fit.invariant_mass_of(:K_S0, :km).constrain_to_nominal_mass_of(:Ds)
  fit.chi2_cut 200
  fit.nominal
end

alg_mode2.note(:mmiss2, "M_miss^2 = (E_cm - E_Ds+ - E_Ds-)^2 - |p_Ds+ + p_Ds-|^2 used to separate D_s^* → D_sγ from D_s^* → D_sπ^0")
alg_mode2.note(:kinematic_fit, "2C kinematic fit constraining both D_s masses to nominal")
alg_mode2.note(:ks_reconstruction, "K_S^0 candidates: secondary vertex fit, flight significance > 2, mass window [0.485,0.510] GeV/c^2")
alg_mode2.note(:ratio, "BF ratio from simultaneous fit to M_miss^2 at all 8 energies")
alg_mode2.note(:cross_feed, "Cross-feed between γ and π^0 channels modeled with MC-derived transfer factors")
alg_mode2.note(:background, "Background from non-D_sD_s events estimated from D_s sidebands in M_BC and ΔE")
alg_mode2.note(:energy_scan, "Data at 8 c.m. energies: 4.130-4.230 GeV")

alg_mode2.apply(sel_mode2)
alg_mode2.execute_on(all_data + all_incMC + sigMC_mode2.to_a)

# ============================================================
# Mode III: D_s^+ → K_S^0 K^+  vs  D_s^- → K_S^0 K^-
# ============================================================

decay_card_mode3 = <<~DECAY
Decay e+ e-
  1.0  D_s*+  D_s*-  VSS;
Enddecay
Decay D_s*+
  1.0  D_s+  gamma  PHSP;
Enddecay
Decay D_s*-
  1.0  D_s-  pi0  PHSP;
Enddecay
Decay D_s+
  1.0  K_S0  K+  PHSP;
Enddecay
Decay D_s-
  1.0  K_S0  K-  PHSP;
Enddecay
Decay K_S0
  1.0  pi+  pi-  PHSP;
Enddecay
Decay pi0
  1.0  gamma  gamma  PHSP;
Enddecay
DECAY

sigMC_mode3 = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name     = "sig_DsStar_DsKsK_DsKsK"
  config.events          = 500_000
  config.decay_card      = decay_card_mode3
  config.cross_section   = :default
end

alg_mode3 = Algorithm.new("DsStarDsKsKModeIII", "00-00-01")
alg_mode3.set_header(["DsStarDsKsKAlg/DsStarDsKsK.h"])
          .set_constant({ "ECMS" => [:double, 4.178] })
          .with_decay_card(decay_card_mode3)

sel_mode3 = Selection.new

sel_mode3.select_track do |t|
  t.nChrp 2
  t.nChrn 2
  t.nTot 4
end

sel_mode3.pid(method: :probability) do |pid|
  pid.identify :kp, :km, against: [:pip, :ep]
  pid.identify :pip, :pim, against: [:kp, :ep]
  pid.prob_cut 0.001
end

sel_mode3.select_photon do |p|
  p.min_energy 0.025
end

sel_mode3.secondary_vertex_fit(:K_S0_a, daughters: [:pip, :pim]) do |v|
  v.mass_window [0.485, 0.510]
  v.flight_significance 2.0
end

sel_mode3.secondary_vertex_fit(:K_S0_b, daughters: [:pip, :pim]) do |v|
  v.mass_window [0.485, 0.510]
  v.flight_significance 2.0
end

sel_mode3.build_virtual_particle(:DsP_cand, from: [:K_S0_a, :kp])
sel_mode3.build_virtual_particle(:DsM_cand, from: [:K_S0_b, :km])

sel_mode3.kinematic_fit([:K_S0_a, :kp, :K_S0_b, :km]) do |fit|
  fit.constrain_four_momentum
  fit.invariant_mass_of(:K_S0_a, :kp).constrain_to_nominal_mass_of(:Ds)
  fit.invariant_mass_of(:K_S0_b, :km).constrain_to_nominal_mass_of(:Ds)
  fit.chi2_cut 200
  fit.nominal
end

alg_mode3.note(:mmiss2, "M_miss^2 = (E_cm - E_Ds+ - E_Ds-)^2 - |p_Ds+ + p_Ds-|^2 used to separate D_s^* → D_sγ from D_s^* → D_sπ^0")
alg_mode3.note(:kinematic_fit, "2C kinematic fit constraining both D_s masses to nominal")
alg_mode3.note(:ks_reconstruction, "Two K_S^0 candidates: each via secondary vertex fit, flight significance > 2, mass window [0.485,0.510] GeV/c^2. Ambiguity resolved by pairing K_S^0 closest in vertex")
alg_mode3.note(:ratio, "BF ratio from simultaneous fit to M_miss^2 at all 8 energies")
alg_mode3.note(:cross_feed, "Cross-feed between γ and π^0 channels modeled with MC-derived transfer factors")
alg_mode3.note(:background, "Background from non-D_sD_s events estimated from D_s sidebands in M_BC and ΔE")
alg_mode3.note(:energy_scan, "Data at 8 c.m. energies: 4.130-4.230 GeV")

alg_mode3.apply(sel_mode3)
alg_mode3.execute_on(all_data + all_incMC + sigMC_mode3.to_a)