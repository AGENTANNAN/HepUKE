# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# Real data at the three energy points (BOSS 703)
data_3810 = DatasetManager.real_data.find("703_3810")   # ≈3807.7 MeV
data_3872 = DatasetManager.real_data.find("703_3872")   # ≈3867.4 MeV (absorbs the paper's 3871.3 MeV point)
data_3900 = DatasetManager.real_data.find("703_3900")   # ≈3896.2 MeV
# Corresponding inclusive MC samples
incMC_3810 = DatasetManager.inclusive_mc.find("703_3810")
incMC_3872 = DatasetManager.inclusive_mc.find("703_3872")
incMC_3900 = DatasetManager.inclusive_mc.find("703_3900")

# Decay card for the signal process: psi(4260) -> pi+ pi- J/psi, J/psi -> e+ e-
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1  pi+  pi-  J/psi    PHSP;
    Enddecay

    Decay J/psi
    1  e+  e-    PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the e+ e- channel only: 200k events, default cross section
exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_703_pipijpsi_ee"
    config.related_dataset = data_3872
    config.events          = 200000
    config.decay_card      = decay_card_signal
    config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "pipiJpsi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.872]})

event_selection = Selection.new
event_selection.select_track {               # Charged track selection
        cos_theta 0.93                        # |cos(theta)| < 0.93
        Vz        10.0                        # |Vz| < 10 cm
        Vr        1.0                         # Vr < 1 cm
        nChrp "==2"                           # exactly 2 positive tracks
        nChrn "==2"                           # exactly 2 negative tracks
        nNet  "==0"                           # net charge summed over all tracks = 0
    }
    .pid(method: :probability) {              # High-momentum lepton identification
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,  # p > 1.0 GeV -> lepton
                                       treat_as_electron_if_energy_above: 1.1   # EMC E > 1.1 GeV -> electron, else muon
        nlp "==1"                             # exactly one positive lepton
        nlm "==1"                             # exactly one negative lepton
    }
    .remove([:lp <= :chrgp, :lm <= :chrgn])   # remove identified leptons from the charged-track lists
    .pid(method: :probability) {              # Remaining tracks identified as pions
        prob_cut 0.001                        # PID probability > 0.001
        identify :pion, against: [:kaon, :proton]
        npip "==1"                            # exactly one pi+
        npim "==1"                            # exactly one pi-
    }
    .kinematic_fit([:pip, :pim, :lp, :lm]) {  # 4C + J/psi mass constraint
        nominal                               # nominal fit: corrected four-momenta are saved
        constrain_four_momentum               # total four-momentum = CMS energy
        invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)  # m(l+l-) = m(J/psi)
        chi2_cut 60                           # chi^2 < 60
    }

# The J/psi -> mu+mu- channel needs its own separate exclusive MC sample (not generated here)
my_algorithm
    .note(:mumu_channel_sample, "The J/psi -> mu+ mu- channel requires a separate exclusive MC sample (psi(4260) -> pi+ pi- J/psi, J/psi -> mu+ mu-); only the e+ e- channel is generated in this spec.")
    .with_decay_card(decay_card_signal)
    .apply(event_selection)

# Run the selection on the three real-data samples, their inclusive MC, and the signal MC
root_files = my_algorithm.execute_on([data_3810, data_3872, data_3900,
                                      incMC_3810, incMC_3872, incMC_3900,
                                      exMC_signal])