### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC sample

# Decay card: ψ(3686) → γ χ_cJ, χ_cJ → 3K_S0 K+ π-, K_S0 → π+π-  (J = 0, 1, 2)
# One card covers all three χ_cJ signal modes, which share the same decay chain.
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    0.3333 gamma chi_c0 P2GC0;
    0.3333 gamma chi_c1 P2GC1;
    0.3334 gamma chi_c2 P2GC2;
    Enddecay

    Decay chi_c0
    1.000 K_S0 K_S0 K_S0 K+ pi- PHSP;
    Enddecay

    Decay chi_c1
    1.000 K_S0 K_S0 K_S0 K+ pi- PHSP;
    Enddecay

    Decay chi_c2
    1.000 K_S0 K_S0 K_S0 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Single 500k-event exclusive MC sample serving the three χ_cJ modes
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_chicJ_3ks_kpi"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "ChiCJ3KsKPi"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  .select_track {                       # Charged track selection
      cos_theta 0.93                    # |cos(theta)| < 0.93
      Vz        10.0                    # |Vz| < 10 cm
      Vr        1.0                     # Vr < 1 cm
      nTot      ">=8"                   # at least 8 charged tracks
  }
  .select_photon {                      # Photon selection
      tdc_emc_start     0               # EMC TDC start time
      tdc_emc_end       14              # EMC TDC end time
      angle_to_track    10.0            # at least 10 degrees from any charged track
      energyThreshold_b 0.025           # 25 MeV in the barrel
      energyThreshold_e 0.050           # 50 MeV in the endcap
      nGam              ">=1"           # at least one photon
  }
  .pid(method: :probability) {          # Probability-based PID, zero probability cut
      prob_cut 0.0                      # no PID probability cut
      identify :kaon, against: [:pion]  # K+ against π+ and K- against π-
      identify :pion, against: [:kaon]  # π+ against K+ and π- against K-
  }
  # K_S0 → π+π- from oppositely charged tracks; three K_S0 are needed.
  # The pions are not PID-constrained (prob_cut = 0), i.e. the K_S0 daughters are
  # used without any PID requirement, as the description asks.
  .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference  # keep the combination closest to the K_S0 mass
      remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  # Nominal 4C kinematic fit to γ K_S0 K_S0 K_S0 K+ π-;
  # the candidate combination with the smallest chi2 (and hence the radiative,
  # highest-energy photon) is selected automatically by the fit.
  .kinematic_fit([:gamma, :K_S0, :K_S0, :K_S0, :kp, :pim]) {
      nominal
      constrain_four_momentum          # 4C energy-momentum constraint
      chi2_cut 50                      # chi2 < 50
  }

# K_S0 quality requirements not exposed by the secondary_vertex_fit DSL surface
algorithm.note(:ks0_selection, "K_S0 candidates are built from oppositely charged track pairs without applying PID: "
                               "|Vz| < 20 cm and a common-vertex fit with chi2 < 200; the decay length from the "
                               "interaction point must exceed 2 sigma; the pi+pi- invariant mass is required to lie "
                               "within |M - 0.498| < 0.012 GeV/c^2 (3 sigma). Mass sidebands "
                               "0.020 < |M - 0.498| < 0.044 GeV/c^2 are retained for background subtraction at the ROOT stage.")

# Generate the algorithm for the three χ_cJ modes described by the (shared) decay card
algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the exclusive signal MC
root_files = algorithm.execute_on([psip_data, psip_incMC, exMC_signal])