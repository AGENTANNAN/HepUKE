# DSL for 2003.13064v2: e+e- -> phi eta' Born cross section measurement
# 20 energy points from 2.05 to 3.08 GeV
# phi -> K+ K-, eta' -> pi+ pi- gamma
# ConExc mode 48 = phi eta'

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# 20 scan energy points (2015 Rscan, BOSS 713)
scan_points = [
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
  DatasetManager.real_data.find("713_Rscan_2500"),
  DatasetManager.real_data.find("713_Rscan_2644"),
  DatasetManager.real_data.find("713_Rscan_2646"),
  DatasetManager.real_data.find("713_Rscan_2800"),
  DatasetManager.real_data.find("713_Rscan_2900"),
  DatasetManager.real_data.find("713_Rscan_2950"),
  DatasetManager.real_data.find("713_Rscan_2981"),
  DatasetManager.real_data.find("713_Rscan_3000"),
  DatasetManager.real_data.find("713_Rscan_3020"),
  DatasetManager.real_data.find("713_Rscan_3080"),
]

incMC_points = scan_points.map do |dp|
  DatasetManager.inclusive_mc.find("#{dp.boss}_#{dp.sample_name}")
end

# ConExc decay card: mode 48 = phi eta'
# No Particle vpho line - DSL auto-injects per energy point
decay_card_phi_etap = <<~DECAYCARD
  Decay vpho
  1 ConExc 48;
  Enddecay
  Decay phi
  1 K+ K- VSS;
  Enddecay
  Decay eta'
  1 pi+ pi- gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_scan = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_phi_etap"
  config.events        = 100_000
  config.decay_card    = decay_card_phi_etap
  config.cross_section = :default
end

# Algorithm for e+e- -> phi eta', phi -> K+K-, eta' -> pi+pi-gamma
alg = Algorithm.new("PhiEtapConExc")
alg.set_header(["PhiEtapConExcAlg/PhiEtapConExc.h"])

# Event selection
event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==2"
    nChrn "==2"
    nTot "==4"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=1"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    identify :kaon, against: [:pion, :proton]
    npip "==1"
    npim "==1"
    nkp "==1"
    nkm "==1"
  end
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg
  .with_decay_card(decay_card_phi_etap)
  .note(:partial_reconstruction, "Paper also uses 3-track + missing-kaon partial reconstruction: 1C kinematic fit constraining K K_miss pi+ pi- gamma with missing mass = m_K, chi2_1C < 20, and phi mass window |M(K K_miss) - m_phi| < 3 sigma. Only the 4-track (both kaons found) case is expressible in DSL.")
  .note(:photon_isr_suppression, "Photon energy > 70 MeV required to suppress ISR background. Not expressible in select_photon; applied at ROOT level.")
  .note(:eta_prime_signal_extraction, "Signal yield extracted from unbinned fit to M(pi+ pi- gamma) distribution. Signal shape from MC convolved with Gaussian resolution. Background parametrized by 2nd-order polynomial.")
  .note(:born_cross_section, "Born cross section = N_obs / (L * (1+delta) * epsilon * B). ISR correction factor (1+delta) and detection efficiency from MC. Iterative procedure until convergence < 1.0%.")
  .note(:resonance_fit, "Line shape fit with coherent sum of phase-space modified Breit-Wigner and phase-space term. Resonance mass and width extracted: M = 2177.5 +/- 4.8 MeV/c^2, Gamma = 149.0 +/- 15.6 MeV.")
  .note(:phi_mass_window, "phi signal region: |M(K K_miss) - m_phi| < 3 sigma. Sideband: 1.050 < M(K K_miss) < 1.130 GeV/c^2 for non-phi background estimation.")
  .note(:interaction_vertex, "Event interaction vertex reconstructed from two pions and one kaon. Not expressible in DSL.")
  .apply(event_selection)

alg.execute_on(scan_points + incMC_points + exMC_scan)