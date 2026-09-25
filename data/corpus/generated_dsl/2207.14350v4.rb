### Dataset preparation ###
# ψ(3686) real data and inclusive MC from the 709_3686 set
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: ψ(3686) → Λ Λ̄ π0, Λ → pπ⁻, Λ̄ → p̄π⁺, π0 → γγ
decay_card_pi0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Lambda0 anti-Lambda0 pi0 PHSP;
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

# Decay card: ψ(3686) → Λ Λ̄ η, η → γγ (same final state p p̄ π⁺π⁻γγ)
decay_card_eta = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Lambda0 anti-Lambda0 eta PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for each signal mode (π0 and η)
exMC_pi0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_Lambdabar_pi0"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_pi0
    config.cross_section   = :default
end

exMC_eta = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_Lambdabar_eta"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_eta
    config.cross_section   = :default
end

### Event selection (BOSS) ###
# Both signal modes (Λ Λ̄ π0 and Λ Λ̄ η) share the identical final state
# p p̄ π⁺π⁻γγ and the identical selection criteria, so a single Algorithm
# instance covers both (shared-final-state rule).
alg_name = "LambdabarGG"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track {               # charged track selection
                    cos_theta 0.93       # |cosθ| < 0.93
                    Vz        10.0       # |Vz| < 10 cm
                    Vr        1.0        # Vr < 1 cm
                    nChrp     "==2"      # two positive tracks (p, π⁺)
                    nChrn     "==2"      # two negative tracks (p̄, π⁻)
                    nNet      "==0"      # net charge zero
                }
               .select_photon {              # photon selection
                    tdc_emc_start     0     # EMC time window [0, 700] ns
                    tdc_emc_end       14
                    angle_to_track    10.0  # ≥ 10° from any charged track
                    energyThreshold_b 0.025 # E > 25 MeV in barrel
                    energyThreshold_e 0.050 # E > 50 MeV in endcap
                    nGam              ">=2" # at least two photons
                }
               .pid(method: :probability) {  # probability PID, prob > 0.001
                    prob_cut 0.001
                    identify :proton, against: [:pion, :kaon]    # p / p̄ vs π, K
                    identify :pion,   against: [:kaon, :proton]  # π⁺ / π⁻ vs K, p
                    nprp "==1"      # p
                    nprm "==1"      # p̄
                    npip "==1"      # π⁺
                    npim "==1"      # π⁻
                }
               .secondary_vertex_fit([:prp, :pim]) {   # Λ → p π⁻
                    build_virtual_particle(:Lambda).by_minimizing_mass_difference
                    remove_used_particle_from_candidate_list
                }
               .secondary_vertex_fit([:prm, :pip]) {   # Λ̄ → p̄ π⁺
                    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                    remove_used_particle_from_candidate_list
                }
               .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) {  # 4C fit
                    nominal                    # nominal fit
                    constrain_four_momentum    # energy-momentum conservation (Λ Λ̄ γγ)
                    chi2_cut 200               # χ² < 200
                    invariant_mass_of(:Lambda).within(1.111, 1.121)      # Λ mass window
                    invariant_mass_of(:Lambda_bar).within(1.111, 1.121)  # Λ̄ mass window
                }

my_Algorithm.with_decay_card(decay_card_pi0).apply(event_selection)
root_files = my_Algorithm.execute_on([psip_data, psip_incMC, exMC_pi0, exMC_eta])