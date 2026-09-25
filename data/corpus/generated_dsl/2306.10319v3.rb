# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi(3097) real data sample
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive MC sample

# Decay card for the signal process (EvtGen format, phase space).
# Representative channel of the four charge-conjugate channels measured:
#   J/psi -> anti-Lambda0 pi+ Sigma- ,  anti-Lambda0 -> anti-p- pi+ ,  Sigma- -> n pi-
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 anti-Lambda0 pi+ Sigma-     PHSP;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+                 PHSP;
    Enddecay

    Decay Sigma-
    1.000 n pi-                       PHSP;
    Enddecay

    End
DECAYCARD

# 500,000 exclusive MC events for J/psi -> anti-Lambda0 pi+ Sigma-
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_Lambdabar_pi_Sigma"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "JpsiLambdaBarPiSigma"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})       # CMS energy 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# Charge-symmetric selection shared by all four channels.
# The Sigma is NOT reconstructed; it is inferred from the recoil mass of the (Lambda/anti-Lambda) pi system.
event_selection = Selection.new
event_selection.select_track {                                 # Charged-track selection
                  cos_theta 0.93                               # |cos(theta)| < 0.93
                  nTot      ">=3"                              # at least three charged tracks
                }
               .pid(method: :probability) {                    # PID, probability method
                  prob_cut 0.001                               # PID probability > 0.001
                  identify :proton, against: [:pion, :kaon]    # p / p-bar
                  identify :pion,   against: [:kaon, :proton]  # pi+/pi-
                }
               .secondary_vertex_fit([:prp, :pim]) {           # Lambda -> p pi-
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               .secondary_vertex_fit([:prm, :pip]) {           # anti-Lambda -> p-bar pi+
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               .partial_miss([3]) {                            # Sigma (recID 3) left undetected, inferred from recoil
                  require_recoil_mass 1.14, 1.24               # recoil mass of (Lambda/anti-Lambda)pi in [1.14, 1.24] GeV/c^2
                }

# BOSS-side procedures that have no dedicated DSL construct
my_algorithm
  .note(:prompt_pion_dca,
        "the prompt pion (the pion not consumed by the Lambda/anti-Lambda secondary vertex) " \
        "is required to have DCA < 10 cm in z and < 1 cm in xy with respect to the interaction " \
        "point; this separates prompt pions from the displaced Lambda/anti-Lambda daughter tracks")
  .note(:secondary_vertex_selection,
        "Lambda/anti-Lambda candidates are required to satisfy |M(p pi) - M(Lambda)| < 5 MeV/c^2 " \
        "and a decay-length significance > 2; among the surviving combinations the candidate " \
        "minimising [M(p pi) - M(Lambda)]^2 is retained per event")
  .note(:pwa_efficiency_reconstruction,
        "for the PWA efficiency a full Sigma -> n pi reconstruction is used, with a 1C kinematic " \
        "fit to the missing neutron (chi^2 < 30) together with a Lambda vertex fit (chi^2 < 30); " \
        "the resulting efficiency is ~31-35% after data/MC corrections")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])