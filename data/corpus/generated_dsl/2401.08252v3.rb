# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset preparation ###
data_3686  = DatasetManager.real_data.find("709_3686")        # ψ(3686) real data at 3.686 GeV
incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")     # corresponding inclusive MC
# Off-resonance continuum data at 3.773 GeV: processed with the SAME selection and
# later scaled / subtracted in ROOT (luminosity & cross-section scaling, n = 1, 1/s).
data_3773  = DatasetManager.real_data.find("712_3773")

# Decay card for the signal ψ(3686) → Ω− K+ anti-Ξ0. The charge-conjugate channel
# (+c.c.) is covered by the charge-symmetric selection below. anti-Ξ0 (and its
# decay anti-Ξ0 → anti-Λ π0) appears in the card but is NOT reconstructed.
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Omega- K+ anti-Xi0    PHSP;
    Enddecay

    Decay Omega-
    1.0000 Lambda0 K-    PHSP;
    Enddecay

    Decay anti-Xi0
    1.0000 anti-Lambda0 pi0    PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-    HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+    HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 1,000,000 events of ψ(3686) → Ω− K+ anti-Ξ0
exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3686_OmegaKXi"
    config.related_dataset = data_3686
    config.events          = 1000000
    config.decay_card      = decay_card_signal
    config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name  = "OmegaKXiRecoil"
omega_alg = Algorithm.new(alg_name)
omega_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         # Inexpressible BOSS-side selection details captured for the systematic layer
         .note(:kaon_selection,
               "the K+ originating from the interaction point is required to have |Vz| < 10 cm
                and |Vxy| < 1 cm; if several K+ candidates survive, the one with the highest
                kaon PID confidence is retained")
         .note(:lambda_mass_window,
               "Λ candidates are required to satisfy the mass window M(p π−) ∈ [1.111, 1.121] GeV
                after the p π− secondary vertex fit")
         .note(:omega_vertex_chi2,
               "the Ω− (Λ K−) secondary vertex fit is required to have vertex χ² < 200")
         .note(:omega_reduced_mass,
               "Ω− candidates are required to satisfy the reduced mass
                M(Λ K−) − M(p π−) + M_Λ_PDG ∈ [1.663, 1.681] GeV")

event_selection = Selection.new
    .select_track {                 # charged-track quality selection
        cos_theta 0.93              # |cosθ| < 0.93
        Vz        10.0              # |Vz| < 10 cm
        Vr        1.0               # Vr < 1 cm
        nChrp     ">=2"             # at least two positive tracks
        nChrn     ">=2"             # at least two negative tracks
    }
    .pid(method: :probability) {    # probability-method PID
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]   # CL_p > CL_K, CL_p > CL_π, CL_p > 0.001
        nprp     ">=1"                              # at least one proton
        identify :kaon,   against: [:pion]          # kaons separated from pions
        identify :pion,   against: [:kaon, :proton] # remaining tracks taken as pions
    }
    .secondary_vertex_fit([:prp, :pim]) {           # Λ → p π−
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:Lambda, :km]) {         # Ω− → Λ K− (common vertex)
        build_virtual_particle(:"Omega-").by_minimizing_mass_difference  # closest to nominal Ω− mass
        remove_used_particle_from_candidate_list
    }
    # Partial reconstruction of Ω− K+; anti-Ξ0 is left unreconstructed and only its
    # recoil mass is required. No global kinematic fit / photon selection is applied.
    .partial_rec([1, 2]) {
        require_recoil_mass 1.282, 1.352            # recoil mass against Ω− K+ = M(anti-Ξ0)
    }

omega_alg.with_decay_card(decay_card_signal).apply(event_selection)

# Real data, inclusive MC, off-resonance continuum data (for background subtraction)
# and the signal exclusive MC.
root_files = omega_alg.execute_on([data_3686, incMC_3686, data_3773, exMC_signal])