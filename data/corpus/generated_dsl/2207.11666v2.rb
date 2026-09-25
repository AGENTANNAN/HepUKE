### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data @ 3.686 GeV (BOSS 7.0.9)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # ψ(3686) inclusive MC

# Decay card: ψ(3686) → Λ Λ̄ ω, Λ → p π−, Λ̄ → p̄ π+, ω → π+ π− π0, π0 → γ γ
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Lambda0 anti-Lambda0 omega   PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-   HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+   HypWK;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0   OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma   PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 500k events for ψ(3686) → Λ Λ̄ ω (same sub-decays)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3686_LambdaLambdabarOmega"
    config.related_dataset = psip_data
    config.events          = 500000
    config.decay_card      = decay_card_signal
    config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "LambdaLambdabarOmega"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # 3.686 GeV centre-of-mass energy
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                    cos_theta 0.93    # |cosθ| < 0.93
                    Vz        10.0    # |Vz| < 10 cm
                    Vr        1.0     # Vr < 1 cm
                    nChrp     ">=3"   # at least six charged tracks in total
                    nChrn     ">=3"
                    nNet      "==0"   # net charge zero
                }
               .select_photon {
                    tdc_emc_start     0      # EMC time window 0–700 ns
                    tdc_emc_end       14
                    energyThreshold_b 0.025  # E > 25 MeV (barrel)
                    energyThreshold_e 0.050  # E > 50 MeV (endcap)
                    angle_to_track    10.0   # > 10° from any charged track
                    nGam              ">=2"  # two photons needed for π0 → γγ
                }
               .pid(method: :probability) {                  # probability-method PID
                    prob_cut 0.001                           # PID probability > 0.001
                    identify :prp,  against: [:pion, :kaon]    # protons vs π/K
                    identify :prm,  against: [:pion, :kaon]    # anti-protons vs π/K
                    identify :pion, against: [:kaon, :proton]  # π+ and π− vs K/p
                    nprp ">=1"
                    nprm ">=1"
                    npip ">=2"    # one π+ for Λ̄, one for ω
                    npim ">=2"    # one π− for Λ, one for ω
                }
               .secondary_vertex_fit([:prp, :pim]) {         # Λ → p π−
                    build_virtual_particle(:Lambda).by_minimizing_mass_difference
                    remove_used_particle_from_candidate_list
                }
               .secondary_vertex_fit([:prm, :pip]) {         # Λ̄ → p̄ π+
                    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                    remove_used_particle_from_candidate_list
                }
                # 5C kinematic fit: 4-momentum conservation + m(γγ) = m(π0)
               .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim, :gamma, :gamma]) {
                    nominal                                                               # nominal fit
                    constrain_four_momentum                                               # 4C
                    invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)                # π0 window 115–150 MeV
                    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # 5th constraint
                    chi2_cut 200                                                          # loose χ² cut (tight cut in ROOT)
                }

# Capture the inexpressible |M−M_Λ| < 5 MeV window on both Λ/Λ̄ reconstructions
my_algorithm
  .note(:lambda_mass_window,
        "Λ and Λ̄ candidates are required to satisfy |M(pπ) − m_Λ| < 5 MeV (m_Λ = 1.1157 GeV). " \
        "by_minimizing_mass_difference selects the combination closest to the nominal Λ mass; " \
        "the explicit ±5 MeV window is imposed as a cut on the reconstructed invariant mass.")
  .with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])