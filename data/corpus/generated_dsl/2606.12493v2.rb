# =============================================================================
# e+e- -> K_S K+ pi-  (K_S -> pi+ pi-)  on the J/psi R-scan
# BOSS part: dataset preparation + event selection up to the final kinematic fit
# =============================================================================

### Dataset description ###
# Scan points currently available in the dataset tables (BOSS 713 R-scan).
# The full 26-point scan (3000.00 - 3119.88 MeV, 440.7 pb-1) will be added once
# the remaining points are available; only the five points below can be used now.
data_points = [
  DatasetManager.real_data.find("713_3080"),
  DatasetManager.real_data.find("713_3020"),
  DatasetManager.real_data.find("713_3000"),
  DatasetManager.real_data.find("713_2981"),
  DatasetManager.real_data.find("713_2950")
]

# Matching inclusive MC for each scan point
inclusive_mc_points = [
  DatasetManager.inclusive_mc.find("713_3080"),
  DatasetManager.inclusive_mc.find("713_3020"),
  DatasetManager.inclusive_mc.find("713_3000"),
  DatasetManager.inclusive_mc.find("713_2981"),
  DatasetManager.inclusive_mc.find("713_2950")
]

# Decay card for the signal process: continuum production e+e- -> K_S K+ pi-
# modelled with the ConExc generator (ISR / vacuum-polarisation corrected,
# mode 9).  The virtual-photon energy is injected automatically per scan point,
# so `Particle vpho` is deliberately omitted (multi-energy scan).
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 9;
    Enddecay

    Decay K_S0
    1.0 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive MC per scan point (one ExclusiveMC per energy point)
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ks_kpi_scan_exclusive_mc"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "KsKPi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.080]})  # representative scan energy (per-point value set by the framework, see note)

event_selection = Selection.new
event_selection
  .select_track {                       # charged track selection
      cos_theta   0.93                  # |cos(theta)| < 0.93
      Vz          10.0                  # |Vz| < 10 cm
      Vr          1.0                   # Vr < 1 cm
      nChrp       ">=2"                 # at least two positive tracks (K+ and K_S pi+)
      nChrn       ">=2"                 # at least two negative tracks (primary pi- and K_S pi-)
      nNet        "==0"                 # net charge zero
  }
  .select_photon {                      # photon selection (quality cut / pi0 veto)
      tdc_emc_start     0               # TDC start time
      tdc_emc_end       14              # TDC end time
      angle_to_track    10.0            # angle to nearest charged track > 10 deg
      energyThreshold_b 0.025           # 25 MeV barrel threshold
      energyThreshold_e 0.050           # 50 MeV endcap threshold
  }
  .pid(method: :probability) {          # PID with the probability method
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]   # separate kaons from pions/protons
      identify :pion, against: [:kaon, :proton]   # separate pions from kaons/protons
  }
  .secondary_vertex_fit([:pip, :pim]) { # common vertex fit for K_S0 -> pi+ pi-
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  # Nominal 4C kinematic fit to K_S0 K+ pi-
  .kinematic_fit([:K_S0, :kp, :pim]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
  }
  # 5C fit: same tracks/combination, adding the K_S0 mass constraint
  # (stored for the downstream PWA performed in ROOT)
  .kinematic_fit([:K_S0, :kp, :pim]) {
      use_track_index_from_nominal_kmfit
      constrain_four_momentum
      invariant_mass_of(:K_S0).constrain_to_nominal_mass_of(:K_S0)
  }

# Inexpressible BOSS-side procedures, captured for the systematic-uncertainty step
my_algorithm
  .note(:per_point_ecms, "Five R-scan points span 2950-3080 MeV; the beam energy must be "
        "set per energy point (3.080 GeV used here as a placeholder). The full 26-point "
        "scan 3000.00-3119.88 MeV (440.7 pb-1) will replace this once the remaining "
        "samples are available.")
  .note(:k_s_daughter_selection, "The two K_S0 daughter pions are used without a PID "
        "requirement and with a looser |Vz| < 20 cm cut; the standard |Vz| < 10 cm cut "
        "applies only to the primary K+ and pi- tracks.")
  .note(:k_s_selection, "After the common-vertex fit (chi2 < 100, mass-difference "
        "minimisation) the K_S0 is required to satisfy |M(pi+pi-) - m_K_S0| < 15 MeV and a "
        "decay length significance > 2 sigma; among surviving candidates the one with the "
        "longest decay length is kept.")
  .note(:background_veto, "gamma-conversion veto: K_S0 -> pi+pi- candidates are rejected "
        "when the opening angle between the two pions is below 30 degrees.")

# Generate the algorithm for the signal decay card and run on all datasets
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on(data_points + inclusive_mc_points + exMCs_signal)