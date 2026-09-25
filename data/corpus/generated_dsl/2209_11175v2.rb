### Dataset description ###
# J/psi (3.097 GeV) data and inclusive MC: sample name "708_3097"
jpsi_data  = DatasetManager.real_data.find("708_3097")     # ~10.09e9 J/psi events (combined)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding inclusive MC

# Signal decay card: J/psi -> gamma K_S0 K_S0 pi0, K_S0 -> pi+ pi-, pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma K_S0 K_S0 pi0    PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Generate 500k exclusive-MC events for J/psi -> gamma K_S0 K_S0 pi0
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "jpsi_gamma_KS0_KS0_pi0_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiGammaKS0KS0pi0"
jpsi_alg = Algorithm.new(alg_name)
jpsi_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})   # J/psi center-of-mass energy (GeV)

event_selection = Selection.new
event_selection
  .select_track {                                   # Charged track selection
    cos_theta 0.93                                  # |cos(theta)| < 0.93
    Vz        20.0                                  # |Vz| < 20 cm
    Vr        1.0                                   # Vr < 1 cm
    nChrp     "==2"                                 # exactly two positive tracks
    nChrn     "==2"                                 # exactly two negative tracks
    nNet      "==0"                                 # net charge zero
  }
  .select_photon {                                  # Photon selection
    tdc_emc_start     0                             # EMC timing start (50 ns units)
    tdc_emc_end       14                            # EMC timing end
    angle_to_track    10.0                          # angle to nearest track > 10 deg
    energyThreshold_b 0.025                         # barrel energy > 25 MeV
    energyThreshold_e 0.050                         # endcap energy > 50 MeV
    nGam              ">=3"                         # at least three photons
  }
  .pid(method: :probability) {                      # Particle ID (probability method)
    prob_cut 0.001                                  # PID probability > 0.001
    identify :pion, against: [:kaon, :proton]       # identify pi+ and pi- against K, p
    npip ">=2"                                      # at least two pi+
    npim ">=2"                                      # at least two pi-
  }
  .secondary_vertex_fit([:pip, :pim]) {             # First K_S0 -> pi+ pi-
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list        # remove used tracks
  }
  .secondary_vertex_fit([:pip, :pim]) {             # Second K_S0 -> pi+ pi-
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {         # Reconstruct pi0 -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # constraint to pi0 mass
    chi2_cut 25                                     # chi2 < 25
    npi0 ">=1"                                      # at least one pi0 candidate
  }
  .kinematic_fit([:gamma, :K_S0, :K_S0, :pi0]) {    # Final fit: 4C + pi0/K_S0 mass constraints (7C)
    nominal                                         # nominal fit
    constrain_four_momentum                         # four-momentum conservation
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # pi0 mass constraint
    chi2_cut 200                                    # chi2 < 200 (loose; tight cut in ROOT)
  }

# Attach the decay card and render the selection into the BOSS algorithm
jpsi_alg.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC, and the exclusive signal MC
root_files = jpsi_alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])