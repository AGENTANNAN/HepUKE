# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(3686) real data (448.1M events) at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # Corresponding inclusive MC sample

# Decay card for the signal process psi(3686) -> gamma chi_cJ, chi_cJ -> K_S0 K_S0 K_S0 K_S0, K_S0 -> pi+ pi-
decay_card_for_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1           PHSP;
    Enddecay

    Decay chi_c1
    1.0000 K_S0 K_S0 K_S0 K_S0    PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC sample: 500k events of psi(3686) -> gamma chi_cJ -> 4K_S0 -> 8pi
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi2S_to_gamma_chicJ_to_4KS0"
  config.related_dataset = psip_data            # Associate with the real dataset for simulation conditions
  config.events          = 500000               # 500k events
  config.decay_card      = decay_card_for_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')   # Save configuration (optional)

### Event selection (BOSS) ###
alg_name = "GammaChicJTo4KS0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})           # CMS energy in GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# BOSS-side procedures that cannot be expressed in formal DSL syntax
my_algorithm
  .note(:helix_correction,
        "MC helix-parameter correction applied to all charged tracks (data and MC) before the 4C kinematic fit")
  .note(:ks_mass_window,
        "K_S0 candidates are required to have |M(pi+pi-) - m_K_S0| < 12 MeV/c^2 after the secondary vertex fit; \
the fit itself only selects the pair closest to the nominal K_S0 mass (by_minimizing_mass_difference)")
  .note(:ks_decay_length,
        "K_S0 candidates are required to have a decay-length significance L/sigma_L > 2 (flight length from the \
secondary vertex relative to its uncertainty)")

# Build the event selection chain
event_selection = Selection.new
  .select_track {                 # Charged track selection
    cos_theta  0.93               # |cos(theta)| < 0.93
    Vz         20.0               # |Vz| < 20 cm
    Vr         1.0                # Vr < 1 cm in the transverse plane
    nChrp      "==4"              # Exactly 4 positively charged tracks (from 4 K_S0 -> pi+ pi-)
    nChrn      "==4"              # Exactly 4 negatively charged tracks
    nNet       "==0"              # Net charge zero
  }
  .select_photon {                # Photon selection (at least one photon from psi(2S) -> gamma chi_cJ)
    tdc_emc_start     0           # EMC TDC start time
    tdc_emc_end       14          # EMC TDC end time
    energyThreshold_b 0.025       # E > 25 MeV in the barrel region
    energyThreshold_e 0.050       # E > 50 MeV in the endcap region
    angle_to_track    10.0        # Angle to nearest charged track > 10 degrees
    nGam              ">=1"       # At least one good photon
  }
  # No PID: all charged tracks are assumed to be pions
  .assign({:chrgp => :pip, :chrgn => :pim})   # Treat all positive tracks as pi+, all negative as pi-

# Form four K_S0 candidates by secondary-vertex fits to pi+ pi- pairs
  .secondary_vertex_fit([:pip, :pim]) {       # First K_S0
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list  # Do not reuse the tracks of this K_S0
  }
  .secondary_vertex_fit([:pip, :pim]) {       # Second K_S0
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:pip, :pim]) {       # Third K_S0
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:pip, :pim]) {       # Fourth K_S0
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Nominal 4C kinematic fit of gamma + 4 K_S0 to the CMS four-momentum
  .kinematic_fit([:gamma, :K_S0, :K_S0, :K_S0, :K_S0]) {
    nominal                     # Mark this as the nominal fit
    constrain_four_momentum     # 4C energy-momentum constraint
    chi2_cut 200                # Loose chi^2 < 200 (tight cut optimised later in ROOT)
  }

# Generate the complete algorithm for the process defined in the decay card
my_algorithm.with_decay_card(decay_card_for_signal).apply(event_selection)

# Execute the same selection on data, inclusive MC and exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])