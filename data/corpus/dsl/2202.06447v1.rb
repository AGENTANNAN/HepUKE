# Paper: 2202.06447v1 — e+e- → K+K-π0 ConExc analysis (R-scan, Born cross section)
# 19 R-scan energy points (2.000–3.080 GeV, BOSS 713, 648 pb−1 total)
# Key features: ConExc mode 8, kaon PID, 4C kinematic fit χ2<65, π0 mass window, ISR φ veto

### Dataset preparation — 19 R-scan points ###
data_2000  = DatasetManager.real_data.find("713_Rscan_2000")
data_2050  = DatasetManager.real_data.find("713_Rscan_2050")
data_2100  = DatasetManager.real_data.find("713_Rscan_2100")
data_2125  = DatasetManager.real_data.find("713_Rscan_2125")
data_2150  = DatasetManager.real_data.find("713_Rscan_2150")
data_2175  = DatasetManager.real_data.find("713_Rscan_2175")
data_2200  = DatasetManager.real_data.find("713_Rscan_2200")
data_2232  = DatasetManager.real_data.find("713_Rscan_2232")
data_2309  = DatasetManager.real_data.find("713_Rscan_2309")
data_2386  = DatasetManager.real_data.find("713_Rscan_2386")
data_2396  = DatasetManager.real_data.find("713_Rscan_2396")
data_2644  = DatasetManager.real_data.find("713_Rscan_2644")
data_2646  = DatasetManager.real_data.find("713_Rscan_2646")
data_2900  = DatasetManager.real_data.find("713_Rscan_2900")
data_2950  = DatasetManager.real_data.find("713_Rscan_2950")
data_2981  = DatasetManager.real_data.find("713_Rscan_2981")
data_3000  = DatasetManager.real_data.find("713_Rscan_3000")
data_3020  = DatasetManager.real_data.find("713_Rscan_3020")
data_3080  = DatasetManager.real_data.find("713_Rscan_3080")

data_points = [data_2000, data_2050, data_2100, data_2125, data_2150, data_2175,
               data_2200, data_2232, data_2309, data_2386, data_2396, data_2644,
               data_2646, data_2900, data_2950, data_2981, data_3000, data_3020, data_3080]

incMC_2000  = DatasetManager.inclusive_mc.find("713_Rscan_2000")
incMC_2050  = DatasetManager.inclusive_mc.find("713_Rscan_2050")
incMC_2100  = DatasetManager.inclusive_mc.find("713_Rscan_2100")
incMC_2125  = DatasetManager.inclusive_mc.find("713_Rscan_2125")
incMC_2150  = DatasetManager.inclusive_mc.find("713_Rscan_2150")
incMC_2175  = DatasetManager.inclusive_mc.find("713_Rscan_2175")
incMC_2200  = DatasetManager.inclusive_mc.find("713_Rscan_2200")
incMC_2232  = DatasetManager.inclusive_mc.find("713_Rscan_2232")
incMC_2309  = DatasetManager.inclusive_mc.find("713_Rscan_2309")
incMC_2386  = DatasetManager.inclusive_mc.find("713_Rscan_2386")
incMC_2396  = DatasetManager.inclusive_mc.find("713_Rscan_2396")
incMC_2644  = DatasetManager.inclusive_mc.find("713_Rscan_2644")
incMC_2646  = DatasetManager.inclusive_mc.find("713_Rscan_2646")
incMC_2900  = DatasetManager.inclusive_mc.find("713_Rscan_2900")
incMC_2950  = DatasetManager.inclusive_mc.find("713_Rscan_2950")
incMC_2981  = DatasetManager.inclusive_mc.find("713_Rscan_2981")
incMC_3000  = DatasetManager.inclusive_mc.find("713_Rscan_3000")
incMC_3020  = DatasetManager.inclusive_mc.find("713_Rscan_3020")
incMC_3080  = DatasetManager.inclusive_mc.find("713_Rscan_3080")

incMC_points = [incMC_2000, incMC_2050, incMC_2100, incMC_2125, incMC_2150, incMC_2175,
                incMC_2200, incMC_2232, incMC_2309, incMC_2386, incMC_2396, incMC_2644,
                incMC_2646, incMC_2900, incMC_2950, incMC_2981, incMC_3000, incMC_3020, incMC_3080]

### ConExc decay card — mode 8 = K+K−π0 (1.34–4.68 GeV) ###
# Particle vpho is OMITTED — DSL auto-injects per-point √s for multi-energy scans
decay_card = <<~DECAYCARD
  Decay vpho
  1 ConExc 8;
  Enddecay
  Decay vhdr
  1 K+ K- pi0 PHSP;
  Enddecay
  End
DECAYCARD

### Exclusive MC — one card over all 19 scan points ###
sig_mc = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_conexc_KKPi0"
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "KKPi0ConExc"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.080]})

event_selection = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=1"
    nChrn ">=1"
    nNet "==0"
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    tdc_emc_start 0
    tdc_emc_end 14
    nGam ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkp ">=1"
    nkm ">=1"
  }
  .kinematic_fit([:kp, :km, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 65
  }
  # ISR φ veto: second 4C fit under K+K−γ hypothesis (Rule T2 — no nominal, no chi2_cut)
  .kinematic_fit([:kp, :km, :gamma]) {
    constrain_four_momentum
  }

my_Algorithm.with_decay_card(decay_card).apply(event_selection)
root_files = my_Algorithm.execute_on(data_points + incMC_points + sig_mc)

# Note: π0 mass window [0.120, 0.150] GeV and ISR φ veto comparison
# χ2(4C_K+K−π0) < χ2(4C_K+K−γ) are applied in the ROOT analysis stage.