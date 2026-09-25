# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# 2.396 GeV point of the 2024 R-scan (one of the 5 points 2.396 / 2.500 / 2.646 / 2.900 / 3.080 GeV)
data_2396  = DatasetManager.real_data.find("713_2396")      # real data, BOSS 713
incMC_2396 = DatasetManager.inclusive_mc.find("713_2396")   # matching inclusive MC

# Decay card for the continuum signal e+e- -> Lambda anti-Lambda.
# ConExc models ISR (up to 2nd order) and the vacuum-polarisation-corrected Born
# cross section sigma_0(m). The DSL auto-detects the `ConExc` token, switches to the
# no-KKMC simulation template and injects `Particle vpho <ECMS> 0.0`, so `Particle vpho`
# is omitted in the card. `ConExc 8` selects the Lambda anti-Lambda mode of the model.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 8;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for continuum Lambda anti-Lambda production (ISR + vacuum polarisation).
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_2396_lambdalambdabar_conexc"
  config.related_dataset = data_2396
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "LambdaFF"                       # time-like Lambda electromagnetic form factor
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 2.396]})   # sqrt(s) = 2.396 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# Single-tag reconstruction variant (BOSS-side) that cannot be expressed in this chain.
my_Algorithm.note(:single_tag_selection,
  "single-tag (only one Lambda reconstructed) selections use a tighter +/-4.7 MeV/c^2 " \
  "p-pi invariant-mass window, a Lambda momentum within +/-3 sigma of expectation, and " \
  "require a combined production+decay vertex chi^2 < 8 for the anti-Lambda.")

event_selection = Selection.new
    .select_track {                                    # charged-track selection
        cos_theta 0.93                                 # |cos(theta)| < 0.93
        Vz        30.0                                 # |Vz| < 30 cm
        Vr        10.0                                 # Vr < 10 cm
        nChrp     ">=2"                                # at least two positive tracks
        nChrn     ">=2"                                # at least two negative tracks
    }
    .pid(method: :probability) {                       # per-track probability PID
        prob_cut 0.001                                 # PID probability > 0.001
        identify :proton, against: [:kaon, :pion]      # p and p-bar
        identify :pion,   against: [:kaon, :proton]    # pi+ and pi-
        nprp ">=1"                                     # at least one p
        nprm ">=1"                                     # at least one p-bar
        npip ">=1"                                     # at least one pi+
        npim ">=1"                                     # at least one pi-
    }
    .secondary_vertex_fit([:prp, :pim]) {              # Lambda -> p pi-
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:prm, :pip]) {              # anti-Lambda -> p-bar pi+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :Lambda_bar]) {           # 4C fit to e+e- -> p pi- p-bar pi+
        nominal
        constrain_four_momentum
        chi2_cut 130
    }

my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the ConExc signal MC
root_files = my_Algorithm.execute_on([data_2396, incMC_2396, exMC_signal])