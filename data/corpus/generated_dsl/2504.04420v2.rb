# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding ψ(3686) inclusive MC sample

# Decay card for the signal ψ(3686) → Ξ− K_S^0 Ω+ (charge conjugate is analysed as well).
# Note the decay-card index order used by the partial reconstruction further below:
#   0 psi(2S), 1 Xi-, 2 K_S0, 3 anti-Omega-, 4 Lambda0, 5 pi-(Xi), 6 p+, 7 pi-(Lambda), 8 pi+, 9 pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Xi- K_S0 anti-Omega- PHSP;
    Enddecay

    Decay Xi-
    1.0000 Lambda0 pi- PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC sample for the signal process (2,000,000 events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_XimKS0Omegap"
  config.related_dataset = psip_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "XimKS0Omegap"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .set_alias({"std::vector<double>" => "Vdouble"})
            # BOSS-side procedures that have no dedicated DSL construct are preserved as notes
            .note(:background_veto, "ψ(3686) -> pi+ pi- J/psi background is vetoed by requiring the recoil mass against the pi+pi- pair assigned to the K_S^0 candidate, M_recoil(pi+pi-) = sqrt((P_cms - p_pip - p_pim)^2), to lie outside [3.09, 3.105] GeV")
            .note(:xi_vertex_fit, "the Xi- candidate is built by a common vertex fit of Lambda and pi- with a positive decay-length requirement; the best Xi- candidate is chosen by minimising sqrt(((M(Lambda pi-)-M_Xi-)/sigma_Xi-)^2 + ((M(p pi-)-M_Lambda)/sigma_Lambda)^2), approximated here by the mass-difference minimisation of the secondary vertex fit")
            .note(:charge_conjugation, "the charge-conjugate channel psi(3686) -> Xi+ K_S^0 Omega- is reconstructed together with the signal mode using the same event selection")

# Event selection chain: tracks -> PID -> Λ / K_S^0 -> Ξ− -> partial reconstruction (missing Ω+)
event_selection = Selection.new
event_selection.select_track {                       # charged track selection
                  cos_theta 0.93                     # |cos(theta)| < 0.93
                  Vz        20.0                     # |Vz| < 20 cm
                  Vr        1.0                      # Vr < 1 cm
                  nChrp     ">=2"                    # at least two positively charged tracks
                  nChrn     ">=3"                    # at least three negatively charged tracks
                }
               .pid(method: :probability) {          # probability PID
                  prob_cut 0.001                     # PID probability > 0.001
                  identify :proton, against: [:kaon, :pion]   # identify p+ (and p-bar)
                  nprp   ">=1"                       # at least one proton
                }
               .remove([:prp <= :chrgp])             # remove proton candidates from the positive track list
               .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks treated as pions
               .secondary_vertex_fit([:prp, :pim]) { # Λ → p π− : secondary vertex fit, min. mass difference
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               .invariant_mass_of(:Lambda).between(1.111, 1.120)   # Λ mass window [1.111, 1.120] GeV
               .secondary_vertex_fit([:pip, :pim]) {  # K_S^0 → π+ π− : secondary vertex fit
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               .invariant_mass_of(:K_S0).between(0.489, 0.506)     # K_S^0 mass window [0.489, 0.506] GeV
               .secondary_vertex_fit([:Lambda, :pim]) {   # Ξ− → Λ π− : vertex fit of Λ and π−
                  build_virtual_particle(:Xi_m).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               .invariant_mass_of(:Xi_m).between(1.313, 1.330)     # Ξ− mass window [1.313, 1.330] GeV
               # Partial reconstruction (recoil-mass technique): tag Ξ− and K_S^0, leave the Ω+
               # (decay-card recID 3) undetected and infer it from the recoil four-momentum.
               # partial_miss replaces the kinematic fit — no kinematic fit is applied.
               .partial_miss([3]) {
                  best_combination_by_mass :Xi_m,   1.32171    # M(Λπ−) → M_Ξ−
                  best_combination_by_mass :Lambda, 1.115683   # M(pπ−) → M_Λ
                  best_combination_by_mass :K_S0,   0.497611   # M(π+π−) → M_K_S^0
                }

# Generate the complete algorithm for the process defined by the decay card
my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute the algorithm on real data, inclusive MC and the signal exclusive MC
root_files = my_Algorithm.execute_on([psip_data, psip_incMC, exMC_signal])