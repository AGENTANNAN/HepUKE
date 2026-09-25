# BESIII Analysis: Measurement of e+e- → KS KL pi0 cross sections
# Paper: 2309.13883v2
# Data: BOSS 713 R-scan, 19 energy points from 2.000 to 3.080 GeV
# ConExc generator for continuum production

# ============================================================
# Datasets — all 19 R-scan energy points
# ============================================================
data_points = [
  DatasetManager.real_data.find("713_Rscan_2000"),
  DatasetManager.real_data.find("713_Rscan_2050"),
  DatasetManager.real_data.find("713_Rscan_2100"),
  DatasetManager.real_data.find("713_Rscan_2125"),
  DatasetManager.real_data.find("713_Rscan_2150"),
  DatasetManager.real_data.find("713_Rscan_2175"),
  DatasetManager.real_data.find("713_Rscan_2200"),
  DatasetManager.real_data.find("713_Rscan_2232"),
  DatasetManager.real_data.find("713_Rscan_2309"),
  DatasetManager.real_data.find("713_Rscan_2386"),
  DatasetManager.real_data.find("713_Rscan_2396"),
  DatasetManager.real_data.find("713_Rscan_2644"),
  DatasetManager.real_data.find("713_Rscan_2646"),
  DatasetManager.real_data.find("713_Rscan_2900"),
  DatasetManager.real_data.find("713_Rscan_2950"),
  DatasetManager.real_data.find("713_Rscan_2981"),
  DatasetManager.real_data.find("713_Rscan_3000"),
  DatasetManager.real_data.find("713_Rscan_3020"),
  DatasetManager.real_data.find("713_Rscan_3080")
]

incMC_points = [
  DatasetManager.inclusive_mc.find("713_Rscan_2000"),
  DatasetManager.inclusive_mc.find("713_Rscan_2050"),
  DatasetManager.inclusive_mc.find("713_Rscan_2100"),
  DatasetManager.inclusive_mc.find("713_Rscan_2125"),
  DatasetManager.inclusive_mc.find("713_Rscan_2150"),
  DatasetManager.inclusive_mc.find("713_Rscan_2175"),
  DatasetManager.inclusive_mc.find("713_Rscan_2200"),
  DatasetManager.inclusive_mc.find("713_Rscan_2232"),
  DatasetManager.inclusive_mc.find("713_Rscan_2309"),
  DatasetManager.inclusive_mc.find("713_Rscan_2386"),
  DatasetManager.inclusive_mc.find("713_Rscan_2396"),
  DatasetManager.inclusive_mc.find("713_Rscan_2644"),
  DatasetManager.inclusive_mc.find("713_Rscan_2646"),
  DatasetManager.inclusive_mc.find("713_Rscan_2900"),
  DatasetManager.inclusive_mc.find("713_Rscan_2950"),
  DatasetManager.inclusive_mc.find("713_Rscan_2981"),
  DatasetManager.inclusive_mc.find("713_Rscan_3000"),
  DatasetManager.inclusive_mc.find("713_Rscan_3020"),
  DatasetManager.inclusive_mc.find("713_Rscan_3080")
]

# ============================================================
# ConExc Decay Card for e+e- → KS KL pi0
# Using vpho (no hard-coded energy — DSL injects per scan point)
# ============================================================
decay_card_KSKLpi0 = <<~DECAYCARD
  Decay vpho
  1 ConExc 9;
  Enddecay
  Decay vhdr
  1 K_S0 K_S0 pi0 PHSP;
  Enddecay
  Decay K_S0
  1 pi+ pi- PHSP;
  Enddecay
  Decay pi0
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Batch signal MC across all 19 scan points
signal_mc_scan = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_KSKLpi0_ConExc"
  config.events        = 100_000
  config.decay_card    = decay_card_KSKLpi0
  config.cross_section = :default
end

# ============================================================
# Analysis Algorithm — shared across all energy points
# KS → pi+ pi-, pi0 → gamma gamma, KL as missing particle
# ============================================================
alg = Algorithm.new("KSKLpi0_Analysis")
alg
  .set_header(["KSKLpi0_Analysis/KSKLpi0.h"])

sel = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==1"
    nChrn "==1"
    nNet "==0"
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
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 30
    npi0 ">=1"
  end
  # ISR veto: 1C fit under gamma KS KL hypothesis
  .kinematic_fit([:K_S0, :gamma]) do
    miss_track_of :K_S0
    constrain_four_momentum
  end
  # Nominal fit: 3C (energy-momentum + pi0 mass + KS mass constraints)
  .kinematic_fit([:K_S0, :pi0]) do
    miss_track_of :K_L0
    invariant_mass_of(:K_S0).constrain_to_nominal_mass_of(:K_S0)
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

alg
  .note(:ks_mass_window, "|M(pi+ pi-) - m_KS| < 12 MeV, decay length > 2*sigma")
  .note(:pi0_mass_window, "|M(gamma gamma) - m_pi0| < 0.015 GeV")
  .note(:isr_veto, "ISR veto using 1C kinematic fit under gamma KS KL hypothesis")
  .note(:kl_reconstruction, "KL reconstructed as missing particle in 3C kinematic fit")
  .note(:pwa, "Partial wave analysis performed in ROOT on selected events")
  .note(:isr_correction, "ISR correction factor f_ISR from ConExc generator log")
  .note(:vp_correction, "Vacuum polarization correction from ConExc generator log")
  .note(:conexc_mode, "ConExc signal MC approximation: e+e- → KS KS pi0 (mode 9), KS treated as KS+KL for cross-section extraction")
  .with_decay_card(decay_card_KSKLpi0)
  .apply(sel)

alg.execute_on(data_points + incMC_points + signal_mc_scan.flatten)