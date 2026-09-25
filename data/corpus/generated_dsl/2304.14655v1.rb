# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/ψ real data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Matching inclusive MC at 3.097 GeV

# Decay card for the signal process: J/ψ → Σ+ Σ̄−, Σ+ → n π+, Σ̄− → n̄ π− (phase space)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Sigma+ anti-Sigma- PHSP;
    Enddecay

    Decay Sigma+
    1.000 n0 pi+ PHSP;
    Enddecay

    Decay anti-Sigma-
    1.000 anti-n0 pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 500k events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_JpsiSigmaSigma"  # Signal MC name
  config.related_dataset = jpsi_data                  # Associated real dataset
  config.events         = 500000                      # 500k signal events
  config.decay_card     = decay_card_signal           # Decay card defining the process
  config.cross_section  = :default                    # Default cross section
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "JpsiSigmaSigma"                 # Algorithm name
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])          # Header file
            .set_constant({"ECMS" => [:double, 3.097]})            # ECMS = 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})       # Type alias

# Build the event selection chain
event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta 0.93        # |cosθ| < 0.93
                  Vz        30.0        # |Vz| < 30 cm
                  Vr        10.0        # Vr < 10 cm in the transverse plane
                  nChrp     "==1"       # exactly one positively charged track
                  nChrn     "==1"       # exactly one negatively charged track
                  nNet      "==0"       # zero net charge
                }
               .pid(method: :probability) {   # PID: separate pions from kaons/protons
                  prob_cut 0.001              # PID probability > 0.001
                  identify :pion, against: [:kaon, :proton]  # π+ and π− at once
                  npip "==1"                  # one π+ candidate
                  npim "==1"                  # one π− candidate
                }
               .kinematic_fit([:pip, :pim]) { # 4C kinematic fit on the π+π− pair
                  nominal                     # nominal fit (its four-momenta are saved)
                  constrain_four_momentum     # four-momentum conservation
                  chi2_cut 200                # χ² < 200 (loose; tight cut applied in ROOT)
                }

# Neutron / anti-neutron reconstruction from EMC showers is a BOSS-side procedure with
# no corresponding DSL primitive — preserve it as a note. The second Σ-mass-constrained
# fit (free neutron momenta) and the Σ+ → p π0 / Σ̄− → p̄ π0 ratio channels are
# ROOT-level steps and are intentionally outside the BOSS spec.
my_Algorithm.note(:neutron_reconstruction, "Neutrons and anti-neutrons reconstructed from EMC showers with E > 600 MeV, lateral moment > 20, angle to any charged track > 10 degrees, and EMC time 0-700 ns; the most energetic shower candidate is taken as the (anti-)neutron.")

# Generate the complete algorithm for the process in the decay card
my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute the algorithm on the specified datasets, producing ROOT files
root_files = my_Algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])