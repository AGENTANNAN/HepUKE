# ============================================================================
# J/psi -> omega eta  (omega -> pi+ pi- pi0 [Dalitz], eta -> gamma gamma,
#                      pi0 -> gamma gamma)  @ sqrt(s) = 3.097 GeV
# BOSS part: dataset preparation + event selection up to (and including) the
# 6C kinematic fit. The Lambda' veto and the omega signal-region cut are
# ROOT-level (post-fit) selections and are therefore not expressed here.
# ============================================================================

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi(3097) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

# Decay card for the signal process (EvtGen syntax / EvtGen particle names)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 omega eta        PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0      OMEGA_DALITZ;
    Enddecay

    Decay eta
    1.0000 gamma gamma      PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma      PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for J/psi -> omega eta with the full decay chain (24 million events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_omegaeta"
  config.related_dataset = jpsi_data          # associated real dataset (for better simulation)
  config.events          = 24_000_000         # 24 million signal MC events
  config.decay_card      = decay_card_signal  # decay card defining the process
  config.cross_section   = :default           # default cross-section
end

### Event selection (BOSS) ###
alg_name = "OmegaEta"                                   # J/psi -> omega eta
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]}) # ECMS = 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {                 # charged track selection
      cos_theta   0.93            # |cos(theta)| < 0.93
      Vz          10.0            # |Vz| < 10 cm
      Vr          1.0             # Vr < 1 cm (transverse plane)
      nChrp       ">=1"           # at least one positively charged track
      nChrn       ">=1"           # at least one negatively charged track
      nNet        "==0"           # net charge zero
  }
  .select_photon {                # photon selection
      tdc_emc_start     0         # TDC start time
      tdc_emc_end       14        # TDC end time
      angle_to_track    10.0      # min angle to nearest charged track (degrees)
      energyThreshold_b 0.025     # E > 25 MeV in the EMC barrel
      energyThreshold_e 0.050     # E > 50 MeV in the EMC endcap
      nGam              ">=4"     # at least four photons (pi0 -> gg and eta -> gg)
  }
  .pid(method: :probability) {    # particle identification, probability method
      prob_cut 0.001              # PID probability > 0.001
      identify :pion, against: [:kaon, :proton]  # pi+ and pi- (charge-conjugation shorthand)
      npip ">=1"                  # at least one pi+
      npim ">=1"                  # at least one pi-
  }
  # 6C kinematic fit over pi+ pi- gamma gamma gamma gamma:
  #   4C (four-momentum conservation) + pi0 mass constraint + eta mass constraint
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma]) {
      nominal                     # nominal fit: its corrected four-momenta are saved
      constrain_four_momentum     # 4C energy-momentum conservation
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 200                # loose chi2 cut (tight cut applied later in ROOT)
  }

# Generate the algorithm for the process defined in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute the algorithm on real data, inclusive MC and signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])