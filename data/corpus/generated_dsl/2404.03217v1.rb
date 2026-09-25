# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC sample

# Decay card for psi(2S) -> pi0 h_c, h_c -> K_S0 K+ pi-, K_S0 -> pi+ pi-, pi0 -> gamma gamma
# (charge-conjugate channel implied by h_c decay)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 K_S0 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC sample: one million events for the signal decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "psip_to_pi0_hc_to_kskpim"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
# Save the MC configuration for later use (optional)
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "PsippPi0Hc"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # CMS energy 3.686 GeV

event_selection = Selection.new
event_selection
    .select_track {                     # charged track selection
        cos_theta 0.93                  # |cos(theta)| < 0.93
        Vz        10.0                  # |Vz| < 10 cm
        Vr        1.0                   # Vr < 1 cm
        nChrp     ">=2"                 # at least two positive tracks
        nChrn     ">=2"                 # at least two negative tracks
    }
    .select_photon {                    # photon selection
        tdc_emc_start     0             # EMC time window [0, 700] ns
        tdc_emc_end       14
        energyThreshold_b 0.025         # > 25 MeV in the EMC barrel
        energyThreshold_e 0.050         # > 50 MeV in the EMC endcap
        angle_to_track    10.0          # > 10 deg from any charged track
        nGam              ">=2"         # at least two photons
    }
    .pid(method: :probability) {        # particle identification (probability method)
        prob_cut 0.001                  # PID probability > 0.001
        identify :kaon, against: [:pion, :proton]   # K+/K- separated from pi and p
        identify :pion, against: [:kaon, :proton]   # pi+/pi- separated from K and p
        nkp   ">=1"                     # at least one K+
        npim  ">=2"                     # at least two pi-
    }
    .secondary_vertex_fit([:pip, :pim]) {   # K_S0 -> pi+ pi- secondary vertex
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list   # remove the used pi+ / pi-
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma, 1C mass constraint
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25                     # chi2 < 25
        npi0 ">=1"                      # at least one pi0 candidate
    }
    .kinematic_fit([:K_S0, :kp, :pim, :pi0]) {  # 5C fit: 4C energy-momentum + 1C pi0 mass
        nominal                         # nominal fit; corrected four-momenta are saved
        constrain_four_momentum         # 4C energy-momentum conservation (CMS)
        invariant_mass_of(:pi0).constrain_to_nominal_mass_of(:pi0)  # 1C pi0 mass constraint
        chi2_cut 200                    # loose cut here; optimal tight cut applied in ROOT
    }

# BOSS-side procedure without a dedicated DSL primitive: K_S0 vertex-quality selection.
my_algorithm.note(:k_s0_selection, "K_S0 candidates are kept only if the pi+pi- invariant mass lies in (0.487, 0.511) GeV/c^2 and the K_S0 flight-length significance L/dL > 2; applied on the virtual K_S0 track after the secondary vertex fit and before the final kinematic fit")

# Attach the decay card and render the full selection into the BOSS algorithm
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC, and the signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])