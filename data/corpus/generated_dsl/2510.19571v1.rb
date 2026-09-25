### Dataset description ###
# ψ(3686) real data and the corresponding inclusive MC sample (2.712e9 ψ(3686) events)
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for ψ(3686) → Ξ0 Ξ̄0, Ξ0 → π0 Λ, Λ → pπ−, and the charge-conjugate chain (π0 → γγ)
decay_card_xi0_xibar0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Xi0 anti-Xi0 PHSP;
    Enddecay

    Decay Xi0
    1.0000 pi0 Lambda0 PHSP;
    Enddecay

    Decay anti-Xi0
    1.0000 pi0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive signal MC for the full Ξ0 Ξ̄0 decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exmc_psip_xi0_xibar0"
    config.related_dataset = psip_data
    config.events          = 100000
    config.decay_card      = decay_card_xi0_xibar0
    config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "Xi0Xibar0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # √s = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {          # Charged track selection
        cos_theta 0.93                  # |cosθ| < 0.93
        Vz        10.0                  # |Vz| < 10 cm
        Vr        10.0                  # Vr < 1 cm (10 mm) in the transverse plane
        nChrp     ">=2"                 # At least two positive tracks (p+, π+)
        nChrn     ">=2"                 # At least two negative tracks (π−, p̄)
        nNet      "==0"                 # Net charge zero
    }
    .select_photon {                    # Photon selection
        tdc_emc_start     0             # TDC start
        tdc_emc_end       14            # TDC end
        angle_to_track    10.0          # Min angle to nearest charged track (degrees)
        energyThreshold_b 0.025         # 25 MeV threshold in the barrel
        energyThreshold_e 0.050         # 50 MeV threshold in the endcap
        nGam              ">=4"         # At least four photons (two π0 → γγ)
    }
    .pid(method: :probability) {        # Probability-method PID
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]   # identify p+ and p̄ at once
        nprp ">=1"                      # At least one proton
        nprm ">=1"                      # At least one anti-proton
    }
    .remove([:prp <= :chrgp])           # Remove protons from the positive charged list
    .remove([:prm <= :chrgn])           # Remove anti-protons from the negative charged list
    .assign({:chrgp => :pip, :chrgn => :pim})   # Remaining tracks treated as pions
    .secondary_vertex_fit([:prp, :pim]) {       # Λ → p π−
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:prm, :pip]) {       # Λ̄ → p̄ π+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    # 4C kinematic fit on Λ Λ̄ γγγγ, with the photon pairs constrained to the π0 mass
    .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma, :gamma, :gamma]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 200
    }

my_algorithm
    .note(:pid_correction_method, "the paper's momentum-dependent proton/pion assignment " \
          "(p > 0.5 GeV/c for the (anti)proton track) is approximated here by the standard " \
          "BESIII probability PID (prob_cut 0.001, identify :proton against [:kaon, :pion]); " \
          "the residual efficiency difference is to be evaluated as a systematic uncertainty")
    .with_decay_card(decay_card_xi0_xibar0).apply(event_selection)

# Execute on real data, inclusive MC and signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])