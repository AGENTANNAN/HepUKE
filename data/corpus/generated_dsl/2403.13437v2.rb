### Dataset preparation ###
# ψ(3686) data and inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# ---- Signal mode I: ψ(3686) → Ω− Ω+, Ω− → Σ0 π− (ΔS=2), Σ0 treated inclusively;
#      tag side Ω+ → Λ K+ with Λ → p π− (charge conjugate)
decay_card_sigma0 = <<~DECAYCARD
  Decay psi(2S)
  1.000  Omega-  Omega+   PHSP;
  Enddecay

  Decay Omega-
  1.000  Sigma0  pi-      PHSP;
  Enddecay

  Decay Sigma0
  1.000  Lambda0  gamma   PHSP;
  Enddecay

  Decay Lambda0
  1.000  p+  pi-         HypWK;
  Enddecay

  Decay Omega+
  1.000  Lambda0  K+      PHSP;
  Enddecay

  End
DECAYCARD

# ---- Signal mode II: ψ(3686) → Ω− Ω+, Ω− → n K− (ΔS=2), neutron not reconstructed;
#      tag side Ω+ → Λ K+ with Λ → p π− (charge conjugate)
decay_card_nK = <<~DECAYCARD
  Decay psi(2S)
  1.000  Omega-  Omega+   PHSP;
  Enddecay

  Decay Omega-
  1.000  n0  K-          PHSP;
  Enddecay

  Decay Omega+
  1.000  Lambda0  K+      PHSP;
  Enddecay

  Decay Lambda0
  1.000  p+  pi-         HypWK;
  Enddecay

  End
DECAYCARD

# ---- Exclusive MC samples: 1.27M events for each signal mode ----
exMC_sigma0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_OmegaToSigma0Pi"
  config.related_dataset = psip_data
  config.events          = 1_270_000
  config.decay_card      = decay_card_sigma0
  config.cross_section   = :default
end

exMC_nK = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_OmegaToNK"
  config.related_dataset = psip_data
  config.events          = 1_270_000
  config.decay_card      = decay_card_nK
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ================= Signal mode I: Ω− → Σ0 π− =================
alg_name_sigma0 = "OmegaToSigma0Pi"
alg_sigma0 = Algorithm.new(alg_name_sigma0)
alg_sigma0.set_header(["#{alg_name_sigma0}Alg/#{alg_name_sigma0}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_sigma0 = Selection.new
  .select_track {                               # charged track quality cuts
    cos_theta 0.93                              # |cosθ| < 0.93
    Vz        10.0                              # |Vz| < 10 cm
    Vr        1.0                               # Vr < 1 cm
    nChrp     ">=2"                             # at least two positive tracks
    nChrn     ">=2"                             # at least two negative tracks
    nNet      "==0"                             # net charge zero
  }
  .pid(method: :probability) {                  # PID, probability method
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # tag Λ → p
    identify :kp,     against: [:pion, :proton] # tag K+
    nprp ">=1"
    nkp  ">=1"
  }
  .remove([:prp <= :chrgp])                     # strip identified particles from charged lists
  .remove([:kp  <= :chrgp])
  .assign({:chrgp => :pip, :chrgn => :pim})     # remaining tracks treated as pions
  .secondary_vertex_fit([:prp, :pim]) {         # tag-side Λ → p π− via secondary vertex fit
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Kp, :pim]) {        # 4C fit on Λ K+ π−
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_sigma0.with_decay_card(decay_card_sigma0).apply(sel_sigma0)

# ================= Signal mode II: Ω− → n K− =================
alg_name_nK = "OmegaToNK"
alg_nK = Algorithm.new(alg_name_nK)
alg_nK.set_header(["#{alg_name_nK}Alg/#{alg_name_nK}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})

sel_nK = Selection.new
  .select_track {                               # charged track quality cuts
    cos_theta 0.93                              # |cosθ| < 0.93
    Vz        10.0                              # |Vz| < 10 cm
    Vr        1.0                               # Vr < 1 cm
    nChrp     "==2"                             # exactly two positive tracks
    nChrn     "==2"                             # exactly two negative tracks
    nNet      "==0"                             # net charge zero
  }
  .pid(method: :probability) {                  # PID, probability method
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # tag Λ → p
    identify :kaon,   against: [:pion, :proton] # tag K+ and signal K−
    nprp ">=1"
    nkp  ">=1"
    nkm  ">=1"
  }
  .remove([:prp <= :chrgp])                     # strip identified particles from charged lists
  .remove([:kp  <= :chrgp])
  .remove([:km  <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})     # remaining tracks treated as pions
  .secondary_vertex_fit([:prp, :pim]) {         # tag-side Λ → p π− via secondary vertex fit
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Kp, :Km]) {         # 4C fit on Λ K+ K−
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_nK.with_decay_card(decay_card_nK).apply(sel_nK)

# Note: Ω+ tag mass M(Λ K+) ∈ [1.664, 1.680] GeV/c² and recoil-mass window
# ∈ [1.652, 1.695] GeV/c² are applied at the ROOT level (post-kinematic-fit), hence
# deliberately not expressed here.

### Execute on datasets ###
root_files_sigma0 = alg_sigma0.execute_on([psip_data, psip_incMC, exMC_sigma0])
root_files_nK     = alg_nK.execute_on([psip_data, psip_incMC, exMC_nK])