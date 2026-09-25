# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # ψ(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # Inclusive ψ(2S) MC
cont_data  = DatasetManager.real_data.find("709_3650")       # 3.65 GeV continuum data for the QED background estimate

# Decay card for the signal process ψ(2S) → p p̄ η, η → γγ (phase space, EvtGen format)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 p+ anti-p- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC sample for the signal channel
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_ppbar_eta_gg"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsipToPpbarEta"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

# Selection chain: tracks → photons → PID → isolated photons → 4C kinematic fit
event_selection = Selection.new
    .select_track {
        cos_theta 0.93     # |cosθ| < 0.93
        Vz        20.0     # |Vz| < 20 cm in beam direction
        Vr        2.0      # Vr < 2 cm in the transverse plane
        nChrp     "==1"    # exactly one positive track
        nChrn     "==1"    # exactly one negative track
        nNet      "==0"    # net charge zero
    }
    .select_photon {
        tdc_emc_start     0      # EMC timing window start
        tdc_emc_end       14     # EMC timing window end
        energyThreshold_b 0.025  # > 25 MeV in the barrel (|cosθ| < 0.80)
        energyThreshold_e 0.050  # > 50 MeV in the endcap (0.86 < |cosθ| < 0.92)
        nGam              ">=2"  # two photons from η → γγ
    }
    .pid(method: :probability) {
        prob_cut  0.001                              # PID probability > 0.001
        identify :proton, against: [:pion, :kaon]    # p+ and p̄ vs pions/kaons
        nprp      "==1"                              # exactly one proton
        nprm      "==1"                              # exactly one antiproton
    }
    .select_isolated_photon {
        angle_to_prp_track 10.0  # photon ≥ 10° from the proton track
        angle_to_prm_track 30.0  # photon ≥ 30° from the antiproton track
        nGam               ">=2" # at least two isolated photons remain
    }
    # 4C kinematic fit to the p p̄ γγ hypothesis; among photon permutations the
    # combination with the smallest χ² is chosen automatically.
    .kinematic_fit([:prp, :prm, :gamma, :gamma]) {
        nominal                 # nominal fit — corrected four-momenta are saved
        constrain_four_momentum # constrain total four-momentum to the CMS energy (4C)
        chi2_cut 20             # χ² < 20
    }

# Generate the complete algorithm for the signal process in the decay card.
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC, continuum data and signal exclusive MC.
# NOTE: post-kinematic-fit cuts (p(p)/p(p̄) > 300 MeV/c, M(p p̄) < 3.067 GeV/c²,
# M(p p̄) < 3.1 − 0.75·M(γγ) GeV/c², |M(γγ) − M_η| < 21 MeV/c², η sidebands at
# 0.43/0.65 GeV) are applied in the ROOT analysis, not here.
root_files = my_algorithm.execute_on([psip_data, psip_incMC, cont_data, exMC_signal])