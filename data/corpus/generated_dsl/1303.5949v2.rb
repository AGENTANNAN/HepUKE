# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
data_4260 = DatasetManager.real_data.find("703_4260")        # Real dataset at 4.260 GeV
incMC_sample = DatasetManager.inclusive_mc.find("703_4260")  # Corresponding inclusive MC sample

# Decay card for e+e- -> pi+ pi- J/psi (J/psi -> e+ e-)
decay_card_signal_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Decay card for e+e- -> pi+ pi- J/psi (J/psi -> mu+ mu-)
decay_card_signal_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Create 100k-event exclusive MC samples for both J/psi decay modes
exMC_signal_ee = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_4260_pipijpsi_ee"
    config.related_dataset = data_4260
    config.events          = 100000
    config.decay_card      = decay_card_signal_ee
    config.cross_section   = :default
end

exMC_signal_mumu = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_4260_pipijpsi_mumu"
    config.related_dataset = data_4260
    config.events          = 100000
    config.decay_card      = decay_card_signal_mumu
    config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "pipiJpsi"
my_algorithm = Algorithm.new(alg_name) # Create a new algorithm object
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"]) # Set the header file of the algorithm
            .set_constant({"ECMS" => [:double, 4.260]})     # Define a constant double variable ECMS = 4.260 GeV
            .set_alias({"std::vector<double>" => "Vdouble"}) # Alias for long type names

# Common event selection chain for both the J/psi -> e+e- and J/psi -> mu+mu- channels
event_selection = Selection.new
event_selection.select_track {           # Charged track selection
                  cos_theta   0.93       # |cos(theta)| < 0.93
                  Vz          10.0       # |Vz| < 10 cm along the beam direction
                  Vr          1.0        # Vr < 1 cm in the transverse plane
                  nChrp       "==2"      # Exactly 2 positively charged tracks
                  nChrn       "==2"      # Exactly 2 negatively charged tracks
                  nNet        "==0"      # Net charge summed over all tracks to be zero
                }
               .pid(method: :probability) {  # PID with hybrid pion + high-momentum lepton identification
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6  # track with p>1.0 GeV/c -> lepton; lepton with EMC eraw>0.6 GeV -> electron, else muon
                  identify :pion, against: [:kaon]  # pi+ and pi- (charge-conjugation shorthand); separate pion from kaon
                  npip   "==1"   # Require exactly one identified pi+
                  npim   "==1"   # Require exactly one identified pi-
                  nlp    "==1"   # Require exactly one identified positive lepton
                  nlm    "==1"   # Require exactly one identified negative lepton
               }
               # 4C kinematic fit of the pi+ pi- l+ l- system to the initial four-momentum
               .kinematic_fit([:pip, :pim, :lp, :lm]) do
                  nominal                  # Mark this fit as the nominal one; its corrected four-momenta are used
                  constrain_four_momentum  # 4C energy-momentum conservation constraint
                  chi2_cut 60              # Require chi^2 < 60
               }

# Both J/psi decay modes share the same selection chain; attach the e+e- decay card
my_algorithm.with_decay_card(decay_card_signal_ee).apply(event_selection)

# Execute the algorithm on the real data, inclusive MC, and both exclusive MC samples
root_files = my_algorithm.execute_on([data_4260, incMC_sample, exMC_signal_ee, exMC_signal_mumu])