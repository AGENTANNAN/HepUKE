# ============================================================================
# Single-baryon-tagged baryon-pair production on the J/ψ and ψ(3686) resonances
#
#   J/ψ(3.097 GeV)  -> Σ(1385)0 Σ̄(1385)0
#   J/ψ(3.097 GeV)  -> Ξ0 Ξ̄0
#   ψ(3686)(3.686)  -> Σ(1385)0 Σ̄(1385)0
#   ψ(3686)(3.686)  -> Ξ0 Ξ̄0
#
#   with  Σ(1385)0 / Ξ0 -> Λ π0 ,  Λ -> p π− ,  π0 -> γ γ
#
# The anti-baryon is NOT reconstructed; it is inferred from the recoil mass
# of the reconstructed π0 Λ system (partial reconstruction).  Four independent
# final states -> four Algorithm objects, each with its own decay card.
# ============================================================================

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/ψ real data, 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # J/ψ inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")     # ψ(3686) real data, 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # ψ(3686) inclusive MC

### Decay cards (EvtGen format) ###
# J/ψ -> Σ(1385)0 Σ̄(1385)0
decay_card_jpsi_Sigma = <<~DECAYCARD
  Decay J/psi
  1.0000 Sigma*0 anti-Sigma*0 PHSP;
  Enddecay

  Decay Sigma*0
  1.0000 Lambda0 pi0 PHSP;
  Enddecay

  Decay anti-Sigma*0
  1.0000 anti-Lambda0 pi0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# J/ψ -> Ξ0 Ξ̄0
decay_card_jpsi_Xi = <<~DECAYCARD
  Decay J/psi
  1.0000 Xi0 anti-Xi0 PHSP;
  Enddecay

  Decay Xi0
  1.0000 Lambda0 pi0 PHSP;
  Enddecay

  Decay anti-Xi0
  1.0000 anti-Lambda0 pi0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ψ(3686) -> Σ(1385)0 Σ̄(1385)0
decay_card_psip_Sigma = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Sigma*0 anti-Sigma*0 PHSP;
  Enddecay

  Decay Sigma*0
  1.0000 Lambda0 pi0 PHSP;
  Enddecay

  Decay anti-Sigma*0
  1.0000 anti-Lambda0 pi0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ψ(3686) -> Ξ0 Ξ̄0
decay_card_psip_Xi = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Xi0 anti-Xi0 PHSP;
  Enddecay

  Decay Xi0
  1.0000 Lambda0 pi0 PHSP;
  Enddecay

  Decay anti-Xi0
  1.0000 anti-Lambda0 pi0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC samples (1,000,000 events per mode) ###
exMC_jpsi_Sigma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_Sigma1385_Sigma1385bar"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_jpsi_Sigma
  config.cross_section   = :default
end

exMC_jpsi_Xi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_Xi0_Xi0bar"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_jpsi_Xi
  config.cross_section   = :default
end

exMC_psip_Sigma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_Sigma1385_Sigma1385bar"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_psip_Sigma
  config.cross_section   = :default
end

exMC_psip_Xi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_Xi0_Xi0bar"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_psip_Xi
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Common selection shared by all four channels:
#   tracks -> photons -> PID (proton + pion) -> π0 (1C Kalman) -> Λ (secondary vertex)
base_selection = Selection.new
    .select_track {                       # charged track quality cuts
        cos_theta 0.93                    # |cosθ| < 0.93
        Vz        20.0                    # |Vz| < 20 cm
        Vr        20.0                    # Vr < 20 cm
        nChrp     ">=1"                   # at least one positive track
        nChrn     ">=1"                   # at least one negative track
    }
    .select_photon {                      # photon selection
        tdc_emc_start     0               # EMC TDC window 0–14
        tdc_emc_end       14
        energyThreshold_b 0.025           # > 25 MeV in the barrel
        energyThreshold_e 0.050           # > 50 MeV in the endcap
        nGam              ">=2"           # at least two photons
    }
    .pid(method: :probability) {          # PID: probability method
        prob_cut 0.001                    # probability threshold 0.001
        identify :proton, against: [:kaon, :pion]   # p / p̄ vs K, π
        identify :pion,   against: [:kaon, :proton] # π+ / π− vs K, p
        nprp ">=1"                        # at least one proton
        npim ">=1"                        # at least one π−
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # π0 -> γγ (1C mass constraint)
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 20                       # χ² < 20
        npi0 ">=1"                        # at least one π0
    }
    .secondary_vertex_fit([:prp, :pim]) {       # Λ -> p π− secondary vertex
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }

# J/ψ -> Σ(1385)0 Σ̄(1385)0 : reconstruct the Σ(1385)0 (π0Λ), infer Σ̄(1385)0 from recoil
sel_jpsi_Sigma = base_selection.dup
    .partial_rec([1]) do                  # recID 1 = Σ(1385)0 (top-level daughter)
        best_combination_by_mass :Sigma_star0, 1.3846   # closest π0Λ combination to Σ(1385)0
        require_recoil_mass 1.3046, 1.4646              # anti-Σ(1385)0 recoil mass ±80 MeV/c²
    end

# J/ψ -> Ξ0 Ξ̄0 : reconstruct the Ξ0 (π0Λ), infer Ξ̄0 from recoil
sel_jpsi_Xi = base_selection.dup
    .partial_rec([1]) do                  # recID 1 = Ξ0 (top-level daughter)
        best_combination_by_mass :Xi0, 1.3149           # closest π0Λ combination to Ξ0
        require_recoil_mass 1.2649, 1.3649              # anti-Ξ0 recoil mass ±50 MeV/c²
    end

# ψ(3686) -> Σ(1385)0 Σ̄(1385)0
sel_psip_Sigma = base_selection.dup
    .partial_rec([1]) do                  # recID 1 = Σ(1385)0
        best_combination_by_mass :Sigma_star0, 1.3846
    end

# ψ(3686) -> Ξ0 Ξ̄0
sel_psip_Xi = base_selection.dup
    .partial_rec([1]) do                  # recID 1 = Ξ0
        best_combination_by_mass :Xi0, 1.3149
    end

### Algorithms ###
# --- J/ψ -> Σ(1385)0 Σ̄(1385)0 ---
alg_jpsi_Sigma = Algorithm.new("JpsiSigma1385Sigma1385bar")
alg_jpsi_Sigma
    .set_header(["JpsiSigma1385Sigma1385barAlg/JpsiSigma1385Sigma1385bar.h"])
    .set_constant({"ECMS" => [:double, 3.097]})
    .note(:lambda_mass_window,
          "J/ψ mode: require |M(pπ−) − M(Λ)| < 5 MeV/c² on the secondary-vertex Λ candidate")
    .note(:signal_window,
          "final window |M(π0Λ) − M(Σ(1385)0)| < 34 MeV/c² applied in ROOT after the partial reconstruction")
alg_jpsi_Sigma.with_decay_card(decay_card_jpsi_Sigma).apply(sel_jpsi_Sigma)
alg_jpsi_Sigma.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_Sigma])

# --- J/ψ -> Ξ0 Ξ̄0 ---
alg_jpsi_Xi = Algorithm.new("JpsiXi0Xi0bar")
alg_jpsi_Xi
    .set_header(["JpsiXi0Xi0barAlg/JpsiXi0Xi0bar.h"])
    .set_constant({"ECMS" => [:double, 3.097]})
    .note(:lambda_mass_window,
          "J/ψ mode: require |M(pπ−) − M(Λ)| < 5 MeV/c² on the secondary-vertex Λ candidate")
    .note(:signal_window,
          "final window |M(π0Λ) − M(Ξ0)| < 10 MeV/c² applied in ROOT after the partial reconstruction")
alg_jpsi_Xi.with_decay_card(decay_card_jpsi_Xi).apply(sel_jpsi_Xi)
alg_jpsi_Xi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_Xi])

# --- ψ(3686) -> Σ(1385)0 Σ̄(1385)0 ---
alg_psip_Sigma = Algorithm.new("PsipSigma1385Sigma1385bar")
alg_psip_Sigma
    .set_header(["PsipSigma1385Sigma1385barAlg/PsipSigma1385Sigma1385bar.h"])
    .set_constant({"ECMS" => [:double, 3.686]})
    .note(:signal_window,
          "final window |M(π0Λ) − M(Σ(1385)0)| < 35 MeV/c² applied in ROOT after the partial reconstruction")
    .note(:background_veto,
          "veto J/ψ transitions: require |M_recoil(π+π−) − M(J/ψ)| > 5 MeV/c² and |M_recoil(π0π0) − M(J/ψ)| > 15 MeV/c²")
alg_psip_Sigma.with_decay_card(decay_card_psip_Sigma).apply(sel_psip_Sigma)
alg_psip_Sigma.execute_on([psip_data, psip_incMC, exMC_psip_Sigma])

# --- ψ(3686) -> Ξ0 Ξ̄0 ---
alg_psip_Xi = Algorithm.new("PsipXi0Xi0bar")
alg_psip_Xi
    .set_header(["PsipXi0Xi0barAlg/PsipXi0Xi0bar.h"])
    .set_constant({"ECMS" => [:double, 3.686]})
    .note(:signal_window,
          "final window |M(π0Λ) − M(Ξ0)| < 11 MeV/c² applied in ROOT after the partial reconstruction")
    .note(:background_veto,
          "veto J/ψ transitions: require |M_recoil(π+π−) − M(J/ψ)| > 5 MeV/c² and |M_recoil(π0π0) − M(J/ψ)| > 15 MeV/c²")
alg_psip_Xi.with_decay_card(decay_card_psip_Xi).apply(sel_psip_Xi)
alg_psip_Xi.execute_on([psip_data, psip_incMC, exMC_psip_Xi])