# ===================== Dataset preparation =====================

# ψ(3686) real data and matching inclusive MC sample (sample name = BOSS version _ CMS energy)
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: inclusive K_S0 production in phase space, K_S0 → π+π−
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 K_S0 X PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 1,000,000 events, ψ(4260) → K_S0 X (phase space), K_S0 → π+π−
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_ks0_inclusive"
  config.related_dataset = psip_data          # generated at the ψ(3686) energy point
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# ===================== Event selection (BOSS) =====================

alg_name = "InclusiveKS0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection
  .select_track {                    # Charged-track quality selection
    cos_theta 0.93                   # |cosθ| < 0.93
    Vz        20.0                   # |Vz| < 20 cm
    Vr        10.0                   # Vr < 10 cm
    nTot      ">=3"                  # at least three good tracks
  }
  .pid(method: :probability) {       # PID by the probability method
    prob_cut 0.001                   # PID probability > 0.001
    identify :pion, against: [:kaon, :proton]  # π+ and π- vs K and p
    npip "==1"                       # exactly one positive pion assigned
    npim "==1"                       # exactly one negative pion assigned
  }
  .remove(:pip) { condition "ep_ratio_of(:pip) > 0.9" }  # reject electron-like π+ (E/p < 0.9)
  .remove(:pim) { condition "ep_ratio_of(:pim) > 0.9" }  # reject electron-like π- (E/p < 0.9)
  .secondary_vertex_fit([:pip, :pim]) {   # Reconstruct K_S0 from the π+π- pair
    build_virtual_particle(:K_S0).by_minimizing_mass_difference  # candidate with mass closest to K_S0
    remove_used_particle_from_candidate_list                      # do not reuse the fitted pions
  }
  .kinematic_fit([:K_S0]) {          # 4C kinematic fit on the K_S0
    nominal                          # nominal fit → corrected four-momenta are used
    constrain_four_momentum          # constrain total four-momentum to the CMS energy
    chi2_cut 200                     # loose χ² cut (tight cut applied in ROOT)
  }

# Procedures on the BOSS side that cannot be expressed in the current DSL
my_algorithm
  .note(:ks_decay_length, "K_S0 candidates are required to have a decay length L > 0.4 cm; the
    reconstruction picks the candidate minimizing the mass difference that also has the longest
    decay length among the π+π- secondary-vertex candidates. No dedicated DSL primitive exists
    for the flight-length cut/ranking, so this is applied in the generated code / cross-checked
    against the secondary-vertex fit output.")
  .note(:additional_track, "require at least one additional good charged track (not used by the
    K_S0) with R_xy < 1 cm and R_z < 10 cm, i.e. a recoil-side track opposite the reconstructed
    K_S0; applied on the leftover charged-track list after the secondary-vertex fit.")
  .note(:ks_mass_window, "the K_S0 signal region |M(π+π-) - M_K_S0| < 11 MeV/c² and the sideband
    regions 21-42 MeV/c² on either side are defined for background estimation and applied as cuts
    / fits in the ROOT analysis stage, not in the BOSS selection.")
  .note(:dataset_scan, "the lineshape uses 22 energy points spanning 3.640-3.701 GeV around the
    ψ(3686); only the 709_3686 data/inclusive-MC sample is configured here, the remaining scan
    points are added through their own dataset entries before the lineshape fit.")

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the exclusive signal MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])