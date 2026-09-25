# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# Eight real-data c.m. energy points: 4.682, 4.699, 4.740, 4.750, 4.781, 4.843, 4.918, 4.951 GeV
data_points = [
  DatasetManager.real_data.find("706_4680"),  # 4.682 GeV
  DatasetManager.real_data.find("706_4700"),  # 4.699 GeV
  DatasetManager.real_data.find("707_4740"),  # 4.740 GeV
  DatasetManager.real_data.find("707_4750"),  # 4.750 GeV
  DatasetManager.real_data.find("707_4780"),  # 4.781 GeV
  DatasetManager.real_data.find("707_4840"),  # 4.843 GeV
  DatasetManager.real_data.find("707_4914"),  # 4.918 GeV
  DatasetManager.real_data.find("707_4946"),  # 4.951 GeV
]

# Corresponding inclusive MC at each energy point
incMC_samples = [
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
  DatasetManager.inclusive_mc.find("707_4740"),
  DatasetManager.inclusive_mc.find("707_4750"),
  DatasetManager.inclusive_mc.find("707_4780"),
  DatasetManager.inclusive_mc.find("707_4840"),
  DatasetManager.inclusive_mc.find("707_4914"),
  DatasetManager.inclusive_mc.find("707_4946"),
]

# Decay card: e+e- -> K_S0 K_S0 psi(3686); psi(3686) -> J/psi pi+ pi-; J/psi -> e+ e-
decay_card_ee = <<~DECAYCARD
    Decay psi(4260)
    1.000 K_S0 K_S0 psi(3686) PHSP;
    Enddecay

    Decay psi(3686)
    1.000 J/psi pi+ pi- PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: same process, but J/psi -> mu+ mu-
decay_card_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.000 K_S0 K_S0 psi(3686) PHSP;
    Enddecay

    Decay psi(3686)
    1.000 J/psi pi+ pi- PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive signal MC for each J/psi decay mode, generated at every energy point
exMCs_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_kskspsi_ee"
  config.events = 100000
  config.decay_card = decay_card_ee
  config.cross_section = :default
end

exMCs_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_kskspsi_mumu"
  config.events = 100000
  config.decay_card = decay_card_mumu
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "KSKSpsi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.682]}) # ECMS is overridden per energy point at execution
            .set_alias({"std::vector<double>" => "Vdouble"})

# The e+e- and mu+mu- channels share one selection chain, differing only in the lepton flavour
event_selection = Selection.new
event_selection.select_track {                 # Charged track selection
                  cos_theta 0.93               # |cos(theta)| < 0.93
                  Vz        10.0               # |Vz| < 10 cm
                  Vr        1.0                # Vr < 1 cm
                  nChrp     ">=1"              # At least one positive track
                  nChrn     ">=1"              # At least one negative track
                  nTot      ">=6"              # At least six tracks in total
                }
               .pid(method: :probability) {    # Lepton identification (high-momentum e/mu)
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.2,   # p > 1.2 GeV/c -> lepton
                                                 treat_as_electron_if_energy_above: 0.6   # EMC energy > 0.6 GeV -> electron, else muon
                  nlp "==1"                    # Exactly one positive lepton
                  nlm "==1"                    # Exactly one negative lepton
                }
               .remove([:lp <= :chrgp, :lm <= :chrgn])  # Remove the identified leptons from the charged lists
               .assign({:chrgp => :pip, :chrgn => :pim}) # Remaining tracks assigned as a pi+ pi- pair
               .secondary_vertex_fit([:pip, :pim]) {     # Secondary vertex fit -> first K_S0
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               .secondary_vertex_fit([:pip, :pim]) {     # Secondary vertex fit -> second K_S0 (all K_S0 K_S0 combinations kept)
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               .kinematic_fit([:K_S0, :K_S0, :lp, :lm]) { # Kinematic fit to K_S0 K_S0 l+ l-
                  nominal                                 # Nominal fit: corrected four-momenta are saved
                  invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi) # Constrain M(l+ l-) to the J/psi mass
                  constrain_four_momentum                 # Four-momentum conservation (missing X left unconstrained)
                  chi2_cut 200                            # Loose chi2 < 200 (tight cut applied later in ROOT)
                }

# Generate the complete algorithm for the signal process in the decay card
my_algorithm.with_decay_card(decay_card_ee).apply(event_selection)

# Execute the algorithm on all datasets, producing ROOT files
root_files = my_algorithm.execute_on(data_points + incMC_samples + exMCs_ee + exMCs_mumu)