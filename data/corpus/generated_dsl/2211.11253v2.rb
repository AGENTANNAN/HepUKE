# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset preparation ###
# Six energy points; the closest available real-data sample and its matching
# inclusive MC are used at each point. No exclusive signal MC is generated:
#   2.2324 GeV -> 2232.4 MeV  (713 R-scan)
#   2.4000 GeV -> 2396.0 MeV  (713 R-scan)
#   2.8000 GeV -> 2800.0 MeV  (713 R-scan)
#   3.0500 GeV -> 3020.0 MeV  (713 R-scan)
#   3.4000 GeV -> 3080.0 MeV  (713 R-scan)
#   3.6710 GeV -> 3686.0 MeV  (709 psi(2S))
data_2232  = DatasetManager.real_data.find("713_Rscan_2232")
incMC_2232 = DatasetManager.inclusive_mc.find("713_Rscan_2232")

data_2400  = DatasetManager.real_data.find("713_Rscan_2396")
incMC_2400 = DatasetManager.inclusive_mc.find("713_Rscan_2396")

data_2800  = DatasetManager.real_data.find("713_Rscan_2800")
incMC_2800 = DatasetManager.inclusive_mc.find("713_Rscan_2800")

data_3050  = DatasetManager.real_data.find("713_Rscan_3020")
incMC_3050 = DatasetManager.inclusive_mc.find("713_Rscan_3020")

data_3400  = DatasetManager.real_data.find("713_Rscan_3080")
incMC_3400 = DatasetManager.inclusive_mc.find("713_Rscan_3080")

data_3671  = DatasetManager.real_data.find("709_3686")
incMC_3671 = DatasetManager.inclusive_mc.find("709_3686")

all_data  = [data_2232, data_2400, data_2800, data_3050, data_3400, data_3671]
all_incMC = [incMC_2232, incMC_2400, incMC_2800, incMC_3050, incMC_3400, incMC_3671]

# Decay card for the inclusive signal topology e+e- -> pi0 / K_S0 + X as it is
# produced in the inclusive hadronic sample generated with the luarlw generator.
decay_card_inclusive = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi0 K_S0 luarlw;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# QED background decay card (events generated with babayaga3.5)
decay_card_qed = <<~DECAYCARD
    Decay vpho
    1.0000 e+ e- PHSP;
    Enddecay
    End
DECAYCARD

# tau+tau- background decay card (events generated with KKMC)
decay_card_tautau = <<~DECAYCARD
    Decay psi(4260)
    1.0000 tau+ tau- PHSP;
    Enddecay

    Decay tau+
    1.0000 e+ nu_e anti-nu_tau PHSP;
    Enddecay

    Decay tau-
    1.0000 pi- pi0 anti-nu_tau PHSP;
    Enddecay

    End
DECAYCARD

# QED background exclusive MC at every energy point
exMC_qed = DatasetManager.create_exclusive_mc_for(all_data) do |config|
    config.sample_name   = "exmc_qed_inclusive"
    config.events        = 100000
    config.decay_card    = decay_card_qed
    config.cross_section = :default
end

# tau+tau- background exclusive MC at the 3.671 GeV point
exMC_tautau = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_tautau_3671"
    config.related_dataset = data_3671
    config.events          = 100000
    config.decay_card      = decay_card_tautau
    config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "InclusivePi0KS0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 2.2324]})   # representative beam energy; each point is run on its own dataset
            .set_alias({"std::vector<double>" => "Vdouble"})

# Steps that cannot be expressed as BOSS/DSL selection code are preserved as notes.
my_algorithm
  .note(:inclusive_signal, "Signal e+e- -> pi0 / K_S0 + X is embedded in the inclusive hadronic sample produced with the luarlw generator; no exclusive signal MC is generated.")
  .note(:background_simulation, "QED background simulated with babayaga3.5 at all six energy points; tau+tau- background at 3.671 GeV simulated with KKMC.")
  .note(:no_kinematic_fit, "No 4C kinematic fit is performed - the Kalman 1C pi0 mass constraint and the K_S0 secondary-vertex fit are the only fits. The offline pi0 angular and K_S0 daughter / decay-length cuts act in place of a kinematic fit.")
  .note(:no_pid, "No PID is applied to the K_S0 daughter tracks; every charged track passing the track selection is used as pi+/pi- for the K_S0 vertex.")
  .note(:pi0_angular_cut, "Offline pi0 acceptance cut on the lab polar angle of the pi0 photons: |cos(theta_gamma)| < 0.8 for p_pi0 < 0.3 GeV/c and < 0.95 for p_pi0 > 0.3 GeV/c.")
  .note(:ks0_decay_length, "K_S0 candidates are required to have decay length > 2 sigma, applied offline after the secondary-vertex fit.")
  .note(:ks0_daughter_cuts, "K_S0 daughter-track cuts |cos(theta)| < 0.93, |Vz| < 30 cm, |Vxy| < 10 cm are applied to the pi+/pi- tracks.")

# Full selection chain: tracks -> photons -> pi0 (1C Kalman) -> K_S0 (secondary vertex).
event_selection = Selection.new
event_selection.select_track {                     # charged-track selection
                    cos_theta        0.93      # |cos(theta)| < 0.93
                    Vz               30.0      # |Vz| < 30 cm (K_S0 daughter requirement)
                    Vr               10.0      # |Vxy| < 10 cm (K_S0 daughter requirement)
                    nChrp            ">=2"     # at least two positively charged tracks
                    nChrn            ">=2"     # at least two negatively charged tracks
                  }
               .select_photon {                     # photon selection
                    tdc_emc_start    0         # EMC time window start
                    tdc_emc_end      14        # EMC time window end (700 ns)
                    angle_to_track   10.0      # angle to the nearest charged track > 10 degrees
                    energyThreshold_b 0.025    # > 25 MeV in the EMC barrel
                    energyThreshold_e 0.050    # > 50 MeV in the EMC endcap
                  }
               .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma, 1C mass constraint
                    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                    chi2_cut         200       # chi2 < 200
                    npi0             ">=1"     # at least one pi0 candidate
                  }
               .assign({:chrgp => :pip, :chrgn => :pim})   # no PID: all remaining tracks treated as pi+/pi-
               .secondary_vertex_fit([:pip, :pim]) {       # K_S0 -> pi+ pi- from a common secondary vertex
                    build_virtual_particle(:K_S0).by_minimizing_mass_difference
                    remove_used_particle_from_candidate_list
                  }
# No 4C kinematic fit is applied, so the selection stops after the K_S0 secondary-vertex fit.

# Generate the complete algorithm for the inclusive topology in the decay card.
my_algorithm.with_decay_card(decay_card_inclusive).apply(event_selection)

# Execute on the six real-data points, their inclusive MC, and the background MC.
root_files = my_algorithm.execute_on(all_data + all_incMC + exMC_qed + [exMC_tautau])