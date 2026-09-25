# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
data_3097 = DatasetManager.real_data.find("708_3097")       # J/psi peak real data at 3.097 GeV
incMC_3097 = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC sample

# Decay card for the signal process J/psi -> gamma A0, A0 -> mu+ mu- (phase-space)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma A0 PHSP;
    Enddecay

    Decay A0
    1.0000 mu+ mu- PHSP;
    Enddecay

    End
DECAYCARD

# Create the exclusive MC sample for the signal process (120k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gammaA0_mumu"
  config.related_dataset = data_3097   # Associate with the J/psi real data sample
  config.events          = 120000      # 120k signal events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiGammaA0MuMu"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])   # header file of the algorithm
            .set_constant({"ECMS" => [:double, 3.097]})      # CMS energy = 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"}) # alias for long type names

# Build the event selection chain (applied identically to data, inclusive MC, exclusive MC)
event_selection = Selection.new
event_selection.select_track {                 # Charged track selection
                  cos_theta   0.93             # |cos(theta)| < 0.93
                  Vz          10.0             # |Vz| < 10 cm
                  Vr          1.0              # Vr < 1 cm in transverse plane
                  nChrp       ">=1"            # At least one positive track
                  nChrn       ">=1"            # At least one negative track
                }
               .select_photon {                # Photon selection
                  tdc_emc_start     0          # TDC start time
                  tdc_emc_end       14         # TDC end time
                  angle_to_track    10.0       # Min angle to nearest charged track (degrees)
                  energyThreshold_b 0.025      # Min energy for EMC barrel region (GeV)
                  energyThreshold_e 0.050      # Min energy for EMC endcap region (GeV)
                  nGam              ">=1"      # At least one photon
                }
               .pid(method: :probability) {    # Particle identification via the probability method
                  # Tracks with p > 0.5 GeV/c are treated as leptons;
                  # lepton with EMC energy > 0.6 GeV -> electron, otherwise muon.
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.5,
                                                 treat_as_electron_if_energy_above: 0.6
                  nlp   ">=1"                  # At least one lepton with positive charge
                  nlm   ">=1"                  # At least one lepton with negative charge
                }
               .kinematic_fit([:gamma, :lp, :lm]) {  # 4C kinematic fit to the gamma l+ l- system
                  nominal                     # Nominal fit — corrected four-momentum is used
                  constrain_four_momentum     # 4C energy-momentum conservation
                  chi2_cut 200                # Loose chi2 < 200 (tight cut applied later in ROOT)
                }

# Generate the algorithm for the process defined in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on the real data, inclusive MC, and exclusive MC samples
root_files = my_algorithm.execute_on([data_3097, incMC_3097, exMC_signal])