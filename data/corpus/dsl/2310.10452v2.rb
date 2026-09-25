# e+e- -> eta pi+ pi- cross section measurement at sqrt(s)=2.00-3.08 GeV
# BESIII, arXiv:2310.10452v2
# ConExc mode 37 = eta pi+ pi-

### Dataset description ###

# R-scan data points: 19 energies from 2.000 to 3.080 GeV
scan_datasets = [
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
  DatasetManager.real_data.find("713_Rscan_3080"),
]

# Inclusive MC for background studies at key energy points
incMC_datasets = [
  DatasetManager.inclusive_mc.find("713_Rscan_2125"),
  DatasetManager.inclusive_mc.find("713_Rscan_2396"),
  DatasetManager.inclusive_mc.find("713_Rscan_2900"),
]

# ConExc decay card for e+e- -> eta pi+ pi- (mode 37)
# Particle vpho is OMITTED — DSL auto-injects per energy point
decay_card_eta_pipi = <<~DECAYCARD
  Decay vpho
  1 ConExc 37;
  Enddecay
  Decay vhdr
  1 pi+ pi- eta PHSP;
  Enddecay
  Decay eta
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive MC for the energy scan (one MC per energy point)
exMC_signal = DatasetManager.create_exclusive_mc_for(scan_datasets) do |config|
  config.sample_name   = "sig_eta_pipi_conexc"
  config.events        = 100000
  config.decay_card    = decay_card_eta_pipi
  config.cross_section = :default
end

# Save configs for each MC point
exMC_signal.each { |m| m.save_to_config(format: :yaml, file_path: "temp_for_test") }

### Event selection (BOSS) ###

alg = Algorithm.new("EtaPiPiScan")
alg.set_header(["EtaPiPiScanAlg/EtaPiPiScan.h"])
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  # Charged track selection: one oppositely charged pion pair
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  end
  # Photon selection: at least 2 photons for eta -> gamma gamma
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  # PID: identify pions
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  end
  # Bhabha removal: reject charged pions with E/p >= 0.8
  .for_each(:pip) do
    define(:e_over_p) { eraw / p }
    where { e_over_p >= 0.8 }
    remove
  end
  .for_each(:pim) do
    define(:e_over_p) { eraw / p }
    where { e_over_p >= 0.8 }
    remove
  end
  # Reconstruct eta -> gamma gamma (1-C Kalman kinematic fit)
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  end
  # 4C kinematic fit with vertex constraint on charged pions
  .kinematic_fit([:pip, :pim, :eta]) do
    nominal
    vertex_fit([0, 1])
    constrain_four_momentum
    chi2_cut 100
  end

alg
  .note(:eta_helicity_cut, "cos(theta_gamma) < 0.95 required for photon in eta helicity frame to suppress photon mis-combination background; applied on corrected four-momenta from nominal kinematic fit")
  .note(:eta_mass_window, "eta mass signal region [0.523, 0.573] GeV/c^2 applied to M(gamma gamma); sideband regions [0.478, 0.503] and [0.593, 0.618] used for background estimation")
  .note(:pwa_fit, "partial-wave analysis (GPUPWA) performed to extract intermediate contributions (rho eta, a2(1320) pi, etc.); fit fractions used to determine signal yields for sub-processes")
  .note(:conexc_range_note, "ConExc mode 37 supports sqrt(s) 1.025-2.975 GeV; points at 3.000, 3.020, 3.080 GeV may exceed the table range")
  .with_decay_card(decay_card_eta_pipi)
  .apply(event_selection)

root_files = alg.execute_on(scan_datasets + incMC_datasets + exMC_signal)