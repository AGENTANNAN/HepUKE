### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # 3.097 GeV J/psi real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

# Decay card for J/psi -> Lambda Sigma0 + c.c., Sigma0 -> gamma Lambda, Lambda -> p pi-
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    0.5 Lambda0 anti-Sigma0 PHSP;
    0.5 anti-Lambda0 Sigma0 PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0 gamma anti-Lambda0 PHSP;
    Enddecay

    Decay Sigma0
    1.0 gamma Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# 100k exclusive MC events for this decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "jpsi_lambda_sigma0_exclusive_mc"
    config.related_dataset = jpsi_data
    config.events          = 100000
    config.decay_card      = decay_card_signal
    config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "JpsiLambdaSigma0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # sqrt(s) = 3.097 GeV

event_selection = Selection.new
event_selection
    .select_track {                       # charged track selection
        cos_theta 0.93                    # |cos(theta)| < 0.93
        Vz        100.0                   # |Vz| < 100 cm
        Vr        10.0                    # Vr < 10 cm
        nChrp     ">=2"                   # >= 2 positive tracks
        nChrn     ">=2"                   # >= 2 negative tracks
        nNet      "==0"                   # net charge zero
    }
    .select_photon {                      # photon selection
        tdc_emc_start     0               # TDC timing window 0-14
        tdc_emc_end       14
        angle_to_track    10.0            # angle to nearest charged track > 10 deg
        energyThreshold_b 0.025           # >= 25 MeV (barrel)
        energyThreshold_e 0.050           # >= 50 MeV (endcap)
        nGam              ">=1"           # at least one photon
    }
    # No PID: positive tracks serve as proton (Lambda) and pi+ (Lambda_bar) candidates;
    # negative tracks serve as pi- (Lambda) and anti-proton (Lambda_bar) candidates.
    .assign({:chrgp => :prp, :chrgn => :pim})
    .assign({:chrgp => :pip, :chrgn => :prm})
    .secondary_vertex_fit([:prp, :pim]) {          # Lambda -> p pi-
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:prm, :pip]) {          # anti-Lambda -> anti-p pi+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    # 4C kinematic fit to gamma Lambda Lambda_bar
    .kinematic_fit([:gamma, :Lambda, :Lambda_bar]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:Lambda).within(1.1107, 1.1207)       # |M(p pi-) - M_Lambda| < 5 MeV
        invariant_mass_of(:Lambda_bar).within(1.1107, 1.1207)   # same window for anti-Lambda
        invariant_mass_of(:gamma, :Lambda).larger_than(1.135)   # M(gamma Lambda) > 1.135 GeV
        chi2_cut 200                                            # loose chi2; tight chi2_4C < 30 applied in ROOT
    }

my_algorithm.note(:candidate_selection, "when several gamma Lambda Lambda_bar combinations pass the selection, the combination minimising sqrt((M(p pi-) - M_Lambda)^2 + (M(anti-p pi+) - M_Lambda)^2) is chosen; the BOSS kinematic fit selects the smallest chi2_4C combination by default")
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])