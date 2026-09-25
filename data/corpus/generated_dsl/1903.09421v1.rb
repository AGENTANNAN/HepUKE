### Dataset description ###
data_2396  = DatasetManager.real_data.find("713_Rscan_2396")       # 2.396 GeV real data (~66.9 pb^-1)
incMC_2396 = DatasetManager.inclusive_mc.find("713_Rscan_2396")    # Matching inclusive MC at 2.396 GeV

# Decay card for the signal process e+e- -> Lambda anti-Lambda (continuum production)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Create the exclusive MC sample for e+e- -> Lambda anti-Lambda
exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_2396_LambdaLambdaBar"
    config.related_dataset = data_2396          # Associated real dataset
    config.events = 500000                       # 500k signal events
    config.decay_card = decay_card_signal
    config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "LambdaLambdaBar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 2.396]})           # sqrt(s) = 2.396 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# Build the event selection chain (no photon is required for this final state)
event_selection = Selection.new
event_selection.select_track {                 # Charged track selection
        cos_theta 0.93                         # |cos(theta)| < 0.93
        Vz        10.0                         # |Vz| < 10 cm
        Vr        1.0                          # Vr < 1 cm
        nChrp     ">=2"                        # At least 2 positive tracks
        nChrn     ">=2"                        # At least 2 negative tracks
        nNet      "==0"                        # Net charge zero
    }
    .pid(method: :probability) {               # PID: probability method
        prob_cut 0.001                         # PID probability > 0.001
        identify :proton, against: [:pion]     # Identify p+ / p- against pions
        nprp ">=1"                             # At least one proton
        nprm ">=1"                             # At least one anti-proton
    }
    # Low-momentum (p < 0.2 GeV/c) tracks are treated as pions, not protons
    .remove(:prp) { condition "three_momentum_of(:prp) < 0.2" }
    .remove(:prm) { condition "three_momentum_of(:prm) < 0.2" }
    .remove([:prp <= :chrgp])                  # Remove protons from positive charged list
    .remove([:prm <= :chrgn])                  # Remove anti-protons from negative charged list
    .assign({:chrgp => :pip, :chrgn => :pim})  # Remaining tracks assigned as pi+ / pi-
    .secondary_vertex_fit([:prp, :pim]) {      # Lambda -> p pi-
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:prm, :pip]) {      # anti-Lambda -> anti-p pi+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    # Nominal 4C kinematic fit to the Lambda anti-Lambda system
    .kinematic_fit([:Lambda, :Lambda_bar]) {
        nominal                                # Nominal fit
        constrain_four_momentum                # 4C energy-momentum constraint
        chi2_cut 50                            # chi2 < 50
    }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([data_2396, incMC_2396, exMC_signal])