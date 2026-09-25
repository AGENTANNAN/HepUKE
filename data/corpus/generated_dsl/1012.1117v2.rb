# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi(3097) real data (3.097 GeV)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC sample

# Decay card for the signal process: J/psi -> gamma eta', eta' -> eta pi+ pi-, eta -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'                PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi-               PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma               PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 500k events of the single mode J/psi -> gamma eta'
# (eta' -> eta pi+ pi- three-body / Dalitz-type decay, eta -> gamma gamma)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gamma_etap_eta_pipi"
  config.related_dataset = jpsi_data      # Associated real dataset for matching conditions
  config.events          = 500000         # Number of events to generate
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiGammaEtaP"
my_algorithm = Algorithm.new(alg_name)                                   # Create algorithm object
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])               # Set the header file
            .set_constant({"ECMS" => [:double, 3.097]})                  # ECMS = 3.097 GeV (J/psi)
            .set_alias({"std::vector<double>" => "Vdouble"})             # Alias for long type names

# Build the event selection chain
event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta 0.93        # |cos(theta)| < 0.93
                  Vz        10.0        # |Vz| < 10 cm along the beam axis
                  Vr        1.0         # Vr < 1 cm in the transverse plane
                  nChrp     "==1"       # Exactly one positively charged track
                  nChrn     "==1"       # Exactly one negatively charged track
                  nNet      "==0"       # Net charge zero
                }
               .select_photon {         # Photon selection
                  tdc_emc_start     0       # EMC TDC start (0)
                  tdc_emc_end       14      # EMC TDC end (14)
                  angle_to_track    20.0    # Isolation: at least 20 deg from nearest charged track
                  energyThreshold_b 0.025   # Min energy in EMC barrel region (25 MeV)
                  energyThreshold_e 0.050   # Min energy in EMC endcap region (50 MeV)
                  nGam              ">=3"   # At least three photons (gamma from J/psi + gamma gamma from eta)
                }
               .pid(method: :probability) {          # PID via the probability method
                  prob_cut 0.001                     # PID probability > 0.001
                  identify :pion, against: [:kaon]   # Identify pion (pi+ and pi-) against kaon
                  npip     ">=1"                     # Require at least one pi+
                }
                # 4C kinematic fit to gamma gamma gamma pi+ pi-
                # The candidate combination (photon assignment) with the minimum chi2 is chosen automatically.
               .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {
                  nominal                 # Nominal fit: its corrected four-momenta are saved
                  constrain_four_momentum # 4C energy-momentum constraint to the CMS system
                  chi2_cut 200            # Loose chi2 cut in BOSS (optimal tight cut applied in ROOT)
                }
                # NOTE: the eta is taken from the gamma-gamma pair whose invariant mass is closest
                # to nominal m(eta), and the eta signal window 0.518-0.578 GeV/c^2 is applied at the
                # ROOT level (post-fit) -- therefore it is not encoded here.

# Attach the decay card and render the selection into the BOSS C++ algorithm
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute the identical selection on real data, inclusive MC and exclusive signal MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])