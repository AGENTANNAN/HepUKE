# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
data_3686 = DatasetManager.real_data.find("709_3686")        # psi(3686) real data set
incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")    # matching inclusive MC

# Decay card for the signal process psi(3686) -> Lambda anti-Lambda (phase space),
# with hyperon weak decays for both Lambda and anti-Lambda.
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
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

# Exclusive signal MC: 5,000,000 events of psi(3686) -> Lambda anti-Lambda
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_to_LambdaLambdabar"
  config.related_dataset = data_3686
  config.events = 5_000_000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "LambdaLambdabarPol"
my_algorithm = Algorithm.new(alg_name)
# Inexpressible BOSS-side procedures captured as notes
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .set_alias({"std::vector<double>" => "Vdouble"})
            .note(:pid_correction_method, "Momentum-dependent PID: tracks with p < 0.6 GeV/c are
              taken as pions and tracks with p > 0.8 GeV/c as protons; the 0.6-0.8 GeV/c overlap
              region is excluded from both hypotheses.")
            .note(:decay_length_significance, "Lambda / anti-Lambda candidates required to have
              decay-length significance L/sigma_L > 2 from the secondary vertex fit.")
            .note(:best_pair_selection, "Best Lambda anti-Lambda pair selected by minimizing
              (M(p pi-) - m_Lambda)^2 + (M(pbar pi+) - m_Lambda)^2, within Lambda mass windows
              [1.108, 1.123] GeV/c^2.")

event_selection = Selection.new
event_selection.select_track {          # Charged track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz 10.0                             # |Vz| < 10 cm
    Vr 1.0                              # Vr < 1 cm in the transverse plane
    nTot ">=4"                          # At least four charged tracks
  }
  .pid(method: :probability) {          # PID by the probability method
    prob_cut 0.001                      # PID probability > 0.001
    identify :proton, against: [:pion]  # protons (p+ and pbar) separated from pions
    nprp ">=1"                          # At least one proton
    nprm ">=1"                          # At least one anti-proton
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])   # Remove identified protons/anti-protons from charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})   # Remaining positive/negative tracks -> pi+ / pi-
  .secondary_vertex_fit([:prp, :pim]) {       # Fit p pi- to a common vertex to form Lambda
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {       # Fit pbar pi+ to a common vertex to form anti-Lambda
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:prp, :pim, :prm, :pip]) {  # Nominal 4C kinematic fit to p pi- pbar pi+
    nominal
    constrain_four_momentum                   # 4C energy-momentum constraint
    chi2_cut 100                              # chi^2 < 100
  }

# Generate the algorithm from the decay card and the selection chain
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC and exclusive signal MC
root_files = my_algorithm.execute_on([data_3686, incMC_3686, exMC_signal])