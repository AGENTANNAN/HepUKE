### Dataset description ###
# J/psi (3.097 GeV) real data and its inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal process: J/psi -> Lambda anti-Lambda,
# Lambda -> n pi0 (neutron undetected), anti-Lambda -> anti-p pi+, pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 n0 pi0 PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive signal MC with exactly this decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_jpsi_lambda_lambdabar"
    config.related_dataset = jpsi_data
    config.events = 100000
    config.decay_card = decay_card_signal
    config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiLambdaLambdabar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # J/psi peak energy
            .set_alias({"std::vector<double>" => "Vdouble"})

# Charged-track, PID, anti-Lambda vertex fit, photon, pi0 Kalman fit,
# then partial reconstruction treating the neutron as missing.
event_selection = Selection.new
event_selection.select_track {
        cos_theta 0.93        # |cos(theta)| < 0.93
        Vz 10.0               # |Vz| < 10 cm
        Vr 1.0                # Vr < 1 cm
        nChrp ">=1"           # at least one positive track
        nChrn ">=1"           # at least one negative track
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]   # identify (anti-)proton vs K and pi
        nprm ">=1"                                  # at least one (anti-)proton for anti-Lambda
    }
    .remove([:prp <= :chrgp, :prm <= :chrgn])       # drop identified (anti-)protons from charged lists
    .assign({:chrgp => :pip, :chrgn => :pim})       # remaining positives -> pi+, negatives -> pi-
    .secondary_vertex_fit([:prm, :pip]) {           # build anti-Lambda from anti-p pi+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        energyThreshold_b 0.025    # 25 MeV (barrel)
        energyThreshold_e 0.050    # 50 MeV (endcap)
        nGam ">=2"                 # at least two photons
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {       # reconstruct pi0 from gamma gamma
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 100
        npi0 ">=1"
    }
    # Neutron undetected: partial reconstruction treats recID 3 (n0, from Lambda0 -> n0 pi0) as missing;
    # all other recIDs are automatically tagged. Recoil is P4_cms - p_Lambda0(reconstructed) - p_antiLambda0.
    .partial_miss([3]) {
        best_combination_by_mass :Lambda_bar, 1.115683   # best anti-Lambda combination near nominal mass
        require_recoil_mass 0.926, 0.957                 # recoil (neutron) mass window
    }

# BOSS-side step with no direct DSL expression: proton-multiplicity vs negative-track-multiplicity veto.
my_algorithm.note(:background_veto,
    "events whose proton multiplicity is <= the negative-track multiplicity are rejected at BOSS level to suppress neutron-related backgrounds")

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and exclusive signal MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])