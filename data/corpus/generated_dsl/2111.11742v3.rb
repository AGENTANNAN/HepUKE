### Dataset description ###
data_3773 = DatasetManager.real_data.find("712_3773")          # ψ(3770) real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")      # Corresponding inclusive MC sample

# Decay card for the signal process e+e- -> Lambda anti-Lambda (EvtGen format / EvtGen particle names)
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
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

# Exclusive MC sample for Lambda anti-Lambda (200M events, phase space + hyperon weak decays)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3773_LambdaLambdabar"
  config.related_dataset = data_3773          # associate with the ψ(3770) real dataset
  config.events = 200_000_000                 # 200M events
  config.decay_card = decay_card_signal
  config.cross_section = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "LambdaLambdabar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.773]})   # double ECMS = 3.773 GeV

# Full selection chain (applied identically to data and MC)
event_selection = Selection.new
event_selection.select_track {          # Charged track selection
    cos_theta   0.93      # |cos(theta)| < 0.93
    Vr          1.0       # Vr < 1 cm in the transverse plane
    Vz          10.0      # |Vz| < 10 cm along the beam direction
    nChrp       "==2"     # exactly two positive tracks
    nChrn       "==2"     # exactly two negative tracks
    nNet        "==0"     # net charge zero
  }
  .pid(method: :probability) {          # Particle identification
    prob_cut   0.001                    # PID probability > 0.001
    identify :proton, against: [:kaon, :pion]  # identify p+ and anti-p- (charge-conjugation shorthand)
    nprp       ">=1"                    # at least one proton
    nprm       ">=1"                    # at least one anti-proton
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])   # remove identified (anti-)protons from charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})   # remaining positive/negative tracks treated as pions
  .secondary_vertex_fit([:prp, :pim]) {       # Reconstruct Lambda -> p pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {       # Reconstruct anti-Lambda -> anti-p pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Nominal 4-momentum-constrained kinematic fit to Lambda anti-Lambda.
  # The +/-5 MeV/c^2 Lambda/anti-Lambda mass windows are applied as pre-selection on the fit candidates.
  .kinematic_fit([:Lambda, :Lambda_bar]) {
    nominal
    invariant_mass_of(:Lambda).within(1.1107, 1.1207)       # M(p pi-) within 5 MeV of the Lambda mass
    invariant_mass_of(:Lambda_bar).within(1.1107, 1.1207)   # M(anti-p pi+) within 5 MeV of the Lambda mass
    constrain_four_momentum                                 # 4C energy-momentum constraint
    chi2_cut 200                                            # chi2 < 200
  }

# Generate the complete algorithm for the process defined in the decay card
my_algorithm
  .note(:background_veto, "Require a positive Lambda (anti-Lambda) decay length in the
    secondary vertex fit to suppress non-Lambda background; not expressible in the DSL.")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC, and signal exclusive MC
root_files = my_algorithm.execute_on([data_3773, incMC_3773, exMC_signal])