### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi(3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding inclusive MC

# Decay card for the signal J/psi -> gamma pi0 pi0 pi0, pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 gamma pi0 pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal process (100k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_gamma3pi0"
  config.related_dataset = jpsi_data
  config.events         = 100000
  config.decay_card     = decay_card_signal
  config.cross_section  = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "Gamma3Pi0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # J/psi CM energy
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  # No charged tracks in the final state
  .select_track {
     nChrp "==0"
     nChrn "==0"
  }
  # Photon selection: EMC energy 25 MeV (barrel) / 50 MeV (endcap), TDC in [0,14], >=7 photons
  .select_photon {
     tdc_emc_start 0
     tdc_emc_end 14
     angle_to_track 10.0
     energyThreshold_b 0.025
     energyThreshold_e 0.050
     nGam ">=7"
  }
  # 1C fit: pair photons to the nominal pi0 mass, chi2 < 10, at least three pi0 candidates
  .kalman_kinematic_fit([:gamma, :gamma]) {
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
     chi2_cut 10
     npi0 ">=3"
  }
  # Reconstruct eta -> gamma gamma, needed for the competing hypothesis below
  .kalman_kinematic_fit([:gamma, :gamma]) {
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
     chi2_cut 25
     neta ">=1"
  }
  # Nominal 7C fit to gamma pi0 pi0 pi0 (4C four-momentum + nominal pi0 masses), chi2 < 40
  .kinematic_fit([:gamma, :pi0, :pi0, :pi0]) {
     nominal
     constrain_four_momentum
     chi2_cut 40
  }
  # Competing gamma eta pi0 pi0 hypothesis: four-momentum constraint only, chi2 stored for veto
  .kinematic_fit([:gamma, :eta, :pi0, :pi0]) {
     constrain_four_momentum
  }

# Capture inexpressible (post-fit) BOSS/analysis-side procedures
my_algorithm
  .note(:background_veto,
        "veto |M(gamma pi0) - M(omega)| < 0.06 GeV/c^2 to suppress the omega -> gamma pi0 background")
  .note(:radiative_photon_veto,
        "veto |M(gamma_r gamma) - M(pi0)| < 0.02 GeV/c^2 against radiative-photon miscombination")
  .note(:f0_980_window,
        "require at least one pi0 pi0 pair with invariant mass in the f0(980) window [0.89, 1.09] GeV/c^2")
  .note(:pwa,
        "partial-wave analysis for M(pi0 pi0 pi0) < 1.6 GeV/c^2 with eta(1405), f1(1285), f1(1420), f1(1510) and nonresonant 0^-+ / 1^++ PHSP components")

# Attach decay card and render the BOSS selection
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])