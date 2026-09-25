# Paper: 2310.10452v2 - e+e- -> eta pi+ pi- at 19 energies 2.000-3.080 GeV
# ConExc generator for continuum; eta -> gamma gamma
# 4C kinematic fit; PWA analysis (ROOT side); BOSS DSL covers selection only

DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

# 19 energy scan points from 2.000 to 3.080 GeV (R-scan)
scan_energies = [
  "713_Rscan_2000", "713_Rscan_2050", "713_Rscan_2100", "713_Rscan_2125",
  "713_Rscan_2150", "713_Rscan_2175", "713_Rscan_2200", "713_Rscan_2232",
  "713_Rscan_2309", "713_Rscan_2386", "713_Rscan_2396",
  "713_Rscan_2644", "713_Rscan_2646",
  "713_Rscan_2900",
  "713_Rscan_2950", "713_Rscan_2981",
  "713_Rscan_3000", "713_Rscan_3020", "713_Rscan_3080"
]

scan_data = scan_energies.map { |s| DatasetManager.real_data.find(s) }
scan_incMC = scan_energies.map { |s| DatasetManager.inclusive_mc.find(s) }

# ConExc decay card for e+e- -> eta pi+ pi- (mode 37)
decay_card_eta_pipi = <<~DECAYCARD
  Decay vpho
  1.0000 eta pi+ pi- ConExc 37;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

signal_mc = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name     = "sig_eta_pipi"
  config.events          = 100_000
  config.decay_card      = decay_card_eta_pipi
  config.cross_section   = :default
end

algorithm = Algorithm.new("EtaPiPi", "00-00-01")

event_selection = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nTot        "==2"      # two oppositely charged pions
    nNet        "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .kinematic_fit([:pip, :pim, :gamma, :gamma]) do
    constrain_four_momentum
    chi2_cut 100
    nominal
  end
  .for_each(:gamma) do
    define(:cos_theta_hel) { cos_theta }
    where { cos_theta_hel < 0.95 }
  end

algorithm
  .set_header(["EtaPiPiAlg/EtaPiPiAlg.h"])
  .note(:generator, "ConExc mode 37 (eta pi+ pi-) for continuum e+e- annihilation. ISR correction incorporated via ConExc.")
  .note(:eta_reco, "eta -> gamma gamma. M(gamma gamma) signal region [0.523, 0.573] GeV/c2. Sideband regions [0.478, 0.503] and [0.593, 0.618] for background estimation (applied in ROOT).")
  .note(:kinematic_fit, "4C kinematic fit under e+e- -> gamma gamma pi+ pi- hypothesis, chi2_4C < 100.")
  .note(:Bhabha_veto, "Bhabha events removed by E/p < 0.8 requirement on charged pions (applied in ROOT).")
  .note(:helicity_cut, "cos(theta_gamma) < 0.95 in helicity frame of eta to suppress mis-combined photons.")
  .note(:pwa, "Partial-wave analysis performed using GPUPWA framework on surviving candidates. Intermediate processes: e+e- -> rho eta, e+e- -> a2(1320) pi, and others. Born cross sections extracted for total and intermediate processes.")
  .note(:cross_section, "Cross section = N_obs / (L * epsilon * B * (1+delta_gamma)). Iterative ISR correction using input cross sections until convergence.")
  .note(:scan_points, "19 energy points from 2.000 to 3.080 GeV, total integrated luminosity 648 pb^-1. R-scan datasets (713).")
  .note(:upper_limit, "Resonant structure in e+e- -> a2(1320)pi cross section observed with 5.5 sigma significance at M = 2044 MeV/c2.")
  .with_decay_card(decay_card_eta_pipi)
  .apply(event_selection)

algorithm.execute_on(scan_data + scan_incMC + signal_mc)