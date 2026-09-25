# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi real data at 3.097 GeV (1.31e9 events)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Matching inclusive MC sample

# Decay card for the baryon-number-violating signal J/psi -> Lambda_c+ e-
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda_c+ e- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 p+ K- pi+ PHSP;
    Enddecay

    End
DECAYCARD

# 100k exclusive-MC events for the signal process (same decay card)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_Lambdac_e"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "LambdacE"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {              # Charged track selection
                  cos_theta 0.93            # |cos(theta)| < 0.93
                  Vz        10.0            # |Vz| < 10 cm
                  Vr        1.0             # Vr < 1 cm
                  nChrp     "==2"           # exactly two positive tracks
                  nChrn     "==2"           # exactly two negative tracks
                  nNet      "==0"           # net charge zero
                }
               .pid(method: :probability) {                   # Per-track probability PID (highest CL)
                  prob_cut 0.001                              # confidence level > 0.001
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.1,
                                                 treat_as_electron_if_energy_above: 0.6
                  identify :pion,   against: [:kaon, :proton] # pi+/pi- vs K, p
                  identify :kaon,   against: [:pion, :proton] # K+/K- vs pi, p
                  identify :proton, against: [:kaon, :pion]   # p/pbar vs K, pi
                  nlm   "==1"       # exactly one e-
                  nprp  "==1"       # exactly one proton
                  nkm   "==1"       # exactly one K-
                  npip  "==1"       # exactly one pi+
                }
               .kinematic_fit([:prp, :km, :pip, :em]) {   # nominal 4C fit, p K- pi+ e- mass assignment
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])