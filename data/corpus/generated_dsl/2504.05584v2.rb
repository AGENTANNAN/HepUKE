### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")     # psi(3770) real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC sample

# Decay card for the signal process e+e- -> Lambda anti-Lambda -> p anti-p pi+ pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 Lambda0 anti-Lambda0    PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                  HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+             HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 500k events for psi(3770) -> Lambda anti-Lambda
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_LambdaLambdabar"
  config.related_dataset = data_3773
  config.events         = 500000
  config.decay_card     = decay_card_signal
  config.cross_section  = :default
end

### Event selection (BOSS) ###
alg_name = "LambdaLambdabar"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})           # ECMS = 3.773 GeV
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
    .select_track {                                      # Charged track selection
        cos_theta  0.93                                  # |cos(theta)| < 0.93
        Vz         10.0                                  # |Vz| < 10 cm
        Vr         1.0                                   # Vr < 1 cm
        nChrp      ">=2"                                 # at least 2 positive tracks
        nChrn      ">=2"                                 # at least 2 negative tracks
        nNet       "==0"                                 # net charge zero
    }
    .pid(method: :probability) {                         # PID by probability method
        prob_cut   0.001                                 # probability > 0.001
        identify   :proton, against: [:kaon, :pion]      # identify proton and anti-proton vs K and pi
        nprp       ">=1"                                 # at least one proton
        nprm       ">=1"                                 # at least one anti-proton
    }
    .remove([:prp <= :chrgp, :prm <= :chrgn])            # remove identified protons/anti-protons from charged lists
    .assign({:chrgp => :pip, :chrgn => :pim})            # remaining positive/negative tracks assigned as pi+/pi-
    .secondary_vertex_fit([:prp, :pim]) {                # Reconstruct Lambda -> p pi-
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:prm, :pip]) {                # Reconstruct anti-Lambda -> anti-p pi+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:Lambda, :Lambda_bar]) {             # Nominal 4C kinematic fit to Lambda anti-Lambda
        nominal
        constrain_four_momentum                          # constrain total four-momentum to CMS energy
        chi2_cut 200                                     # loose chi2 cut in BOSS; tighter chi2 < 100 applied later in ROOT
    }

# Generate algorithm from the decay card and apply the event selection
alg.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC, and exclusive signal MC
root_files = alg.execute_on([data_3773, incMC_3773, exMC_signal])