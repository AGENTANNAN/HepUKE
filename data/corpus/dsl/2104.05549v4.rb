# DSL for: e+e- → φη at √s = 2.00–3.08 GeV
# Paper: 2104.05549v4
# Classification: Ordinary (Algorithm + Selection), ConExc generator
# Data: 713 R-scan data, 22 energy points
# Process: e+e- → φη, φ → K+K-, η → γγ

# All 22 R-scan energy points
data_2000 = DatasetManager.real_data.find("713_2000")
data_2050 = DatasetManager.real_data.find("713_2050")
data_2100 = DatasetManager.real_data.find("713_2100")
data_2125 = DatasetManager.real_data.find("713_2125")
data_2150 = DatasetManager.real_data.find("713_2150")
data_2175 = DatasetManager.real_data.find("713_2175")
data_2200 = DatasetManager.real_data.find("713_2200")
data_2232 = DatasetManager.real_data.find("713_2232")
data_2309 = DatasetManager.real_data.find("713_2309")
data_2386 = DatasetManager.real_data.find("713_2386")
data_2396 = DatasetManager.real_data.find("713_2396")
data_2400 = DatasetManager.real_data.find("713_2400")
data_2644 = DatasetManager.real_data.find("713_2644")
data_2646 = DatasetManager.real_data.find("713_2646")
data_2900 = DatasetManager.real_data.find("713_2900")
data_2950 = DatasetManager.real_data.find("713_2950")
data_2981 = DatasetManager.real_data.find("713_2981")
data_3000 = DatasetManager.real_data.find("713_3000")
data_3020 = DatasetManager.real_data.find("713_3020")
data_3050 = DatasetManager.real_data.find("713_3050")
data_3060 = DatasetManager.real_data.find("713_3060")
data_3080 = DatasetManager.real_data.find("713_3080")

incMC_2000 = DatasetManager.inclusive_mc.find("713_2000")
incMC_2050 = DatasetManager.inclusive_mc.find("713_2050")
incMC_2100 = DatasetManager.inclusive_mc.find("713_2100")
incMC_2125 = DatasetManager.inclusive_mc.find("713_2125")
incMC_2150 = DatasetManager.inclusive_mc.find("713_2150")
incMC_2175 = DatasetManager.inclusive_mc.find("713_2175")
incMC_2200 = DatasetManager.inclusive_mc.find("713_2200")
incMC_2232 = DatasetManager.inclusive_mc.find("713_2232")
incMC_2309 = DatasetManager.inclusive_mc.find("713_2309")
incMC_2386 = DatasetManager.inclusive_mc.find("713_2386")
incMC_2396 = DatasetManager.inclusive_mc.find("713_2396")
incMC_2400 = DatasetManager.inclusive_mc.find("713_2400")
incMC_2644 = DatasetManager.inclusive_mc.find("713_2644")
incMC_2646 = DatasetManager.inclusive_mc.find("713_2646")
incMC_2900 = DatasetManager.inclusive_mc.find("713_2900")
incMC_2950 = DatasetManager.inclusive_mc.find("713_2950")
incMC_2981 = DatasetManager.inclusive_mc.find("713_2981")
incMC_3000 = DatasetManager.inclusive_mc.find("713_3000")
incMC_3020 = DatasetManager.inclusive_mc.find("713_3020")
incMC_3050 = DatasetManager.inclusive_mc.find("713_3050")
incMC_3060 = DatasetManager.inclusive_mc.find("713_3060")
incMC_3080 = DatasetManager.inclusive_mc.find("713_3080")

all_datasets = [data_2000, data_2050, data_2100, data_2125, data_2150, data_2175,
                data_2200, data_2232, data_2309, data_2386, data_2396, data_2400,
                data_2644, data_2646, data_2900, data_2950, data_2981, data_3000,
                data_3020, data_3050, data_3060, data_3080]

all_incMC = [incMC_2000, incMC_2050, incMC_2100, incMC_2125, incMC_2150, incMC_2175,
             incMC_2200, incMC_2232, incMC_2309, incMC_2386, incMC_2396, incMC_2400,
             incMC_2644, incMC_2646, incMC_2900, incMC_2950, incMC_2981, incMC_3000,
             incMC_3020, incMC_3050, incMC_3060, incMC_3080]

# ConExc signal MC: mode 23 = φη, φ → K+K-, η → γγ
# create_exclusive_mc_for handles multi-energy scan points
# Note: omit "Particle vpho" line for multi-energy ConExc — DSL auto-injects per energy
sig_mc_phi_eta = DatasetManager.create_exclusive_mc_for(all_datasets) do |config|
  config.sample_name     = "sig_eemm_to_phi_eta"
  config.events          = 2_500_000
  config.decay_card      = <<~DECAYCARD
    ConExc 23
    End
  DECAYCARD
  config.cross_section   = :default
end

# ===========================================================================
# Algorithm: e+e- → φη
# ===========================================================================

alg = Algorithm.new("eeToPhiEta", version: "00-00-01")
alg.set_header(["eeToPhiEtaAlg/eeToPhiEta.h"])
  .set_constant({ "ECMS" => [:double, 2.125] })

sel = Selection.new("PhiEtaSelection")

# --- Track selection: exactly 2 charged tracks, opposite charge ---
sel.select_track do |t|
  t.cos_theta 0.93
  t.Vz 10.0
  t.Vr 1.0
  t.nTot 2
  t.nNet 0
end

# --- Kaon PID: both tracks identified as kaons ---
sel.pid(method: :probability) do |p|
  p.identify :kp, :km, against: [:pip, :pim]
  p.prob_cut 0.0
end

# --- Photon selection: at least 2 photons, E > 70 MeV ---
sel.select_photon do |g|
  g.energyThreshold_b 0.070
  g.energyThreshold_e 0.070
  g.nGam 2
end

# Photon timing: EMC time within 700 ns
sel.note(:photon_timing,
  "EMC time within [0, 700] ns of event start time required for all photon candidates.")

# --- 4C kinematic fit: e+e- → K+K-γγ ---
# For events with >2 photons, combination with smallest χ² is kept
sel.kinematic_fit([:kp, :km, :gamma, :gamma]) do |k|
  k.constrain_four_momentum
  k.chi2_cut 100
  k.nominal
end

# Best photon combination selection
sel.note(:best_photon_combination,
  "For events with >2 good photon candidates, the combination with smallest χ²_4C is retained. Implemented via photon multiplicity and best-combination logic in BOSS algorithm.")

# --- Competing hypothesis veto: e+e- → K+K-γ (ISR φ) ---
# Rule T2: non-nominal kinematic fit for veto
sel.kinematic_fit([:kp, :km, :gamma]) do |k|
  k.constrain_four_momentum
  # χ²(K+K-γ) compared to χ²(K+K-γγ) in ROOT; events with lower χ²(K+K-γ) are rejected
end

sel.note(:isr_phi_veto,
  "Events rejected if χ²(e+e-→K+K-γ) < χ²(e+e-→K+K-γγ) to suppress ISR φ background. Competing χ² comparison done in ROOT.")

# --- Post-fit signal region (applied in ROOT) ---
sel.note(:signal_region,
  "|M(γγ) - m_η| < 30 MeV/c² and 0.98 < M(K+K-) < 1.08 GeV/c² applied in ROOT analysis. Signal yield extracted via unbinned maximum likelihood fit to M(K+K-) with P-wave relativistic Breit-Wigner ⊗ Gaussian.")

# --- Inexpressible: fit range and MC shape details ---
sel.note(:yield_extraction,
  "Cross section from Born cross section formula: σBorn = N_sig / (L · ε · (1+δ) · 1/(1-Π)² · B). Efficiency and radiative corrections determined iteratively using ConExc MC with sampled line-shape parameters.")

alg.with_decay_card(nil)  # ConExc card embedded in create_exclusive_mc_for
alg.apply(sel)
alg.execute_on(all_datasets + all_incMC + [sig_mc_phi_eta].flatten)