### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")        # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")     # Matching inclusive MC

# Decay card for the double-tag signal J/psi -> Xi0 anti-Xi0 with the full subsequent decays
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000  Xi0  anti-Xi0        PHSP;
    Enddecay

    Decay Xi0
    1.0000  gamma  Sigma0        PHSP;
    Enddecay

    Decay Sigma0
    1.0000  gamma  Lambda0       PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-             HypWK;
    Enddecay

    Decay anti-Xi0
    1.0000  anti-Lambda0  pi0    PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-  pi+        HypWK;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma        PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the full signal chain (1M events, single mode)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_Xi0Xibar0"
  config.related_dataset = jpsi_data        # Associated real dataset
  config.events          = 1_000_000        # 1M events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "Xi0Xibar0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # J/psi CMS energy (GeV)

event_selection = Selection.new
  .select_track {                 # Charged track selection
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        100.0               # Loose |Vz| < 100 cm
    Vr        10.0                # Vr < 10 mm
    nChrp     ">=2"               # At least 2 positive tracks
    nChrn     ">=2"               # At least 2 negative tracks
    nNet      "==0"               # Net charge zero
  }
  .select_photon {                # Photon selection
    tdc_emc_start     0           # TDC start
    tdc_emc_end       14          # TDC end
    angle_to_track    10.0        # >10 deg from nearest charged track
    energyThreshold_b 0.025       # Barrel E > 25 MeV
    energyThreshold_e 0.050       # Endcap E > 50 MeV
    nGam ">=4"                    # At least 4 photons
  }
  .pid(method: :probability) {    # PID with the probability method
    prob_cut 0.001                # PID probability > 0.001
    identify :proton, against: [:kaon, :pion]    # p / p-bar vs K, pi
    identify :pion,   against: [:kaon, :proton]  # pi / pi-bar vs K, p
    nprp ">=1"                    # At least one p
    nprm ">=1"                    # At least one p-bar
    npip ">=1"                    # At least one pi+
    npim ">=1"                    # At least one pi-
  }
  .secondary_vertex_fit([:prp, :pim]) {          # Lambda -> p pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {          # anti-Lambda -> p-bar pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {      # pi0 -> gamma gamma (1-C mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 30                                  # chi2 < 30
    npi0 ">=1"                                   # At least one pi0
  }
  .secondary_vertex_fit([:Lambda_bar, :pi0]) {   # Tag anti-Xi0 -> anti-Lambda pi0
    build_virtual_particle(:Xi_bar0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Xi_bar0, :Lambda, :gamma, :gamma]) {   # Final kinematic fit
    nominal                                      # Nominal fit (corrected 4-momenta used)
    constrain_four_momentum                      # 4C constraint to CMS energy-momentum
    invariant_mass_of(:gamma, :Lambda).constrain_to_nominal_mass_of(:Sigma0)  # gamma Lambda -> Sigma0
    chi2_cut 200                                 # Loose chi2 < 200 (tight cut in ROOT)
  }

# Background evaluation cannot be expressed in the event-selection DSL
my_algorithm.note(:background_evaluation,
  "Dominant Xi0 -> Lambda pi0 background with a misidentified photon (pi0 faking a photon) " \
  "is evaluated via MC and sidebands; alpha_gamma is subsequently extracted from the decay angular distribution")

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])