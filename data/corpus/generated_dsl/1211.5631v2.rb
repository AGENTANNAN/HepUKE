# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

### Decay cards (EvtGen format) ###
# ψ' → pbar K+ Σ0 , Σ0 → γ Λ , Λ → p π-
decay_card_sigma0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 anti-p- K+ Sigma0    PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0    PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-    HypWK;
    Enddecay

    End
DECAYCARD

# ψ' → γ χ_c0 , χ_c0 → pbar K+ Λ , Λ → p π-
decay_card_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0    PHSP;
    Enddecay

    Decay chi_c0
    1.0000 anti-p- K+ Lambda0    PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-    HypWK;
    Enddecay

    End
DECAYCARD

# ψ' → γ χ_c1 , χ_c1 → pbar K+ Λ , Λ → p π-
decay_card_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1    PHSP;
    Enddecay

    Decay chi_c1
    1.0000 anti-p- K+ Lambda0    PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-    HypWK;
    Enddecay

    End
DECAYCARD

# ψ' → γ χ_c2 , χ_c2 → pbar K+ Λ , Λ → p π-
decay_card_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2    PHSP;
    Enddecay

    Decay chi_c2
    1.0000 anti-p- K+ Lambda0    PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-    HypWK;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples (200k events each) ###
exMC_sigma0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_pbarK_Sigma0"
    config.related_dataset = psip_data
    config.events          = 200_000
    config.decay_card      = decay_card_sigma0
    config.cross_section   = :default
end

exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_gamma_chic0_pbarKLambda"
    config.related_dataset = psip_data
    config.events          = 200_000
    config.decay_card      = decay_card_chic0
    config.cross_section   = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_gamma_chic1_pbarKLambda"
    config.related_dataset = psip_data
    config.events          = 200_000
    config.decay_card      = decay_card_chic1
    config.cross_section   = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_gamma_chic2_pbarKLambda"
    config.related_dataset = psip_data
    config.events          = 200_000
    config.decay_card      = decay_card_chic2
    config.cross_section   = :default
end

### Event selection — Σ0 mode (ψ' → pbar K+ Σ0) ###
alg_name_sigma0 = "PsipToPbarKSigma0"
alg_sigma0 = Algorithm.new(alg_name_sigma0)
alg_sigma0.set_header(["#{alg_name_sigma0}Alg/#{alg_name_sigma0}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .note(:tight_vertex_track_cuts,
                "the pbar and K+ tracks entering the pbar/K+ common-vertex fit are " \
                "required to satisfy |Vz| < 10 cm and Vr < 1 cm, tighter than the base " \
                "event-level |Vz| < 30 cm, |Vr| < 15 cm; a per-particle track-quality " \
                "cut the current DSL cannot express")

# Build the full Σ0 event-selection chain
sel_sigma0 = Selection.new
    .select_track {                     # charged track selection
        cos_theta 0.93                  # |cosθ| < 0.93
        Vz        30.0                  # |Vz| < 30 cm
        Vr        15.0                  # |Vr| < 15 cm
        nChrp     ">=2"                 # at least two positive tracks
        nChrn     ">=2"                 # at least two negative tracks
    }
    .select_photon {                    # photon selection
        tdc_emc_start     0             # EMC TDC in [0, 14]
        tdc_emc_end       14
        energyThreshold_b 0.025         # 25 MeV in barrel
        energyThreshold_e 0.050         # 50 MeV in endcap
        nGam              ">=1"         # at least one photon
    }
    .pid(method: :probability) {        # probability-method PID
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]     # p / pbar vs K, π
        identify :kaon,   against: [:pion, :proton]   # K+ / K- vs π, p
        identify :pion,   against: [:kaon, :proton]   # π+ / π- vs K, p
        nprp "==1"                      # exactly one p
        nprm "==1"                      # exactly one pbar
        nkp  "==1"                      # exactly one K+
        npim "==1"                      # exactly one π-
    }
    .secondary_vertex_fit([:prp, :pim]) {             # Λ → p π- secondary vertex
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:gamma, :prm, :kp, :Lambda]) {    # 4C fit on γ pbar K+ Λ
        nominal                                        # nominal fit
        vertex_fit([1, 2])                             # pbar (idx 1), K+ (idx 2) common vertex
        constrain_four_momentum                        # 4C energy-momentum constraint
        chi2_cut 200                                   # χ² < 200
        invariant_mass_of(:prp, :pim).within(1.1087, 1.1227)      # |M(pπ-) − M_Λ| < 7 MeV
        invariant_mass_of(:gamma, :Lambda).within(1.1776, 1.2076) # |M(γΛ) − M_Σ0| < 15 MeV
    }

alg_sigma0.with_decay_card(decay_card_sigma0).apply(sel_sigma0)
alg_sigma0.execute_on([psip_data, psip_incMC, exMC_sigma0])

### Event selection — χ_cJ modes (ψ' → γ χ_cJ → γ pbar K+ Λ), shared by J = 0, 1, 2 ###
alg_name_chicJ = "PsipToGammaChicJPbarKLambda"
alg_chicJ = Algorithm.new(alg_name_chicJ)
alg_chicJ.set_header(["#{alg_name_chicJ}Alg/#{alg_name_chicJ}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .note(:tight_vertex_track_cuts,
               "the pbar and K+ tracks entering the pbar/K+ common-vertex fit are " \
               "required to satisfy |Vz| < 10 cm and Vr < 1 cm, tighter than the base " \
               "event-level |Vz| < 30 cm, |Vr| < 15 cm; a per-particle track-quality " \
               "cut the current DSL cannot express")

# Identical final state (γ pbar K+ p π-) and identical selection for χ_c0, χ_c1, χ_c2
sel_chicJ = Selection.new
    .select_track {
        cos_theta 0.93
        Vz        30.0
        Vr        15.0
        nChrp     ">=2"
        nChrn     ">=2"
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=1"
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :kaon,   against: [:pion, :proton]
        identify :pion,   against: [:kaon, :proton]
        nprp "==1"
        nprm "==1"
        nkp  "==1"
        npim "==1"
    }
    .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:gamma, :prm, :kp, :Lambda]) {
        nominal
        vertex_fit([1, 2])                             # pbar (idx 1), K+ (idx 2) common vertex
        constrain_four_momentum
        chi2_cut 200
        invariant_mass_of(:prp, :pim).within(1.1087, 1.1227)      # same Λ mass window
        invariant_mass_of(:gamma, :Lambda).out_of(1.1776, 1.2076) # veto |M(γΛ) − M_Σ0| < 15 MeV
    }

# One algorithm for all three χ_cJ decays (shared final state and selection)
alg_chicJ.with_decay_card(decay_card_chic0).apply(sel_chicJ)
alg_chicJ.execute_on([psip_data, psip_incMC, exMC_chic0, exMC_chic1, exMC_chic2])