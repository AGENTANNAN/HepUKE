### Dataset description ###
# J/psi scan region — resolved 713 R-scan real-data points used here
rscan_3080 = DatasetManager.real_data.find("713_3080")
rscan_3020 = DatasetManager.real_data.find("713_3020")
rscan_3000 = DatasetManager.real_data.find("713_3000")
rscan_2981 = DatasetManager.real_data.find("713_2981")
rscan_2950 = DatasetManager.real_data.find("713_2950")
rscan_points = [rscan_3080, rscan_3020, rscan_3000, rscan_2981, rscan_2950]

# Matched inclusive MC at each scan point
incMC_3080 = DatasetManager.inclusive_mc.find("713_3080")
incMC_3020 = DatasetManager.inclusive_mc.find("713_3020")
incMC_3000 = DatasetManager.inclusive_mc.find("713_3000")
incMC_2981 = DatasetManager.inclusive_mc.find("713_2981")
incMC_2950 = DatasetManager.inclusive_mc.find("713_2950")
incMCs = [incMC_3080, incMC_3020, incMC_3000, incMC_2981, incMC_2950]

# Decay card: ConExc continuum / J/psi lineshape (ISR up to second order) for
# e+e- -> K_S0 K+ pi-, K_S0 -> pi+ pi-.
# ConExc is auto-detected from the literal token; the `Particle vpho <ECMS>`
# line is injected per scan energy point, so it is omitted here (multi-energy).
decay_card_signal = <<~DECAYCARD
  Decay vpho
  1.0000 ConExc 8;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# 200k-event exclusive MC at each scan point (same card, one MC per energy point)
exMCs = DatasetManager.create_exclusive_mc_for(rscan_points) do |config|
  config.sample_name   = "exmc_kskpi_scan"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMCs.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) ###
alg_name = "KsKPi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.080]})   # representative scan energy
            .set_alias({"std::vector<double>" => "Vdouble"})
            .note(:multi_energy_ecms, "R-scan over 3080/3020/3000/2981/2950 MeV: a single ECMS constant is set to a representative 3.080 GeV; per-energy beam energy is injected by the ConExc simulation template")

event_selection = Selection.new
event_selection
  .select_track {                       # charged track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     ">=2"                     # at least two positive tracks
    nChrn     ">=2"                     # at least two negative tracks
    nNet      "==0"                     # net charge zero
  }
  .select_photon {                      # photon candidate selection (no multiplicity cut)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0              # >= 10 deg from any charged track
    energyThreshold_b 0.025            # > 25 MeV in the barrel
    energyThreshold_e 0.050            # > 50 MeV in the endcap
  }
  .pid(method: :probability) {          # PID: probability method, 0.001 cut
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]   # K+ / K- vs pi, p
    identify :pion, against: [:kaon, :proton]   # pi+ / pi- vs K, p
  }
  .secondary_vertex_fit([:pip, :pim]) { # K_S0 -> pi+ pi-
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list   # remove used pions from candidate list
  }
  # Nominal 4C kinematic fit: e+e- -> K_S0 K+ pi-
  .kinematic_fit([:K_S0, :kp, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200                         # loose cut; tight cut applied in ROOT (Rule T3)
  }
  # 5C fit: adds the constraint M(pi+ pi-) = m(K_S0); stores its chi2
  .kinematic_fit([:K_S0, :kp, :pim]) {
    constrain_four_momentum
    invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
    chi2_cut 200
  }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on the scan points, their matched inclusive MC, and the per-point signal MC
root_files = my_algorithm.execute_on(rscan_points + incMCs + exMCs)