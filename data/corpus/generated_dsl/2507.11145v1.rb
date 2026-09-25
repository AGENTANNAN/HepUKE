### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")        # J/ψ real data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")     # Corresponding inclusive MC

# ------------------------------------------------------------------
# Mode I: J/ψ → γ Λ Λ̄ , via Λ(1520) → γ Λ ;  Λ → p π⁻ , Λ̄ → p̄ π⁺
# ------------------------------------------------------------------
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda(1520)0 anti-Lambda0          PHSP;
    Enddecay

    Decay Lambda(1520)0
    1.0000 gamma Lambda0                       PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                              HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                         HypWK;
    Enddecay

    End
DECAYCARD

# ------------------------------------------------------------------
# Mode II: J/ψ → γ γ Λ Λ̄ , via Λ(1520) → γ Σ0 , Σ0 → γ Λ
# ------------------------------------------------------------------
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda(1520)0 anti-Lambda0          PHSP;
    Enddecay

    Decay Lambda(1520)0
    1.0000 gamma Sigma0                        PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0                       PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                              HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                         HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 100k events for each mode
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gLambdaLambdabar_modeI"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_ggLambdaLambdabar_modeII"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# ==================================================================
# Mode I event selection
# ==================================================================
alg_modeI = Algorithm.new("JpsiGammaLambdaLambdabar")
alg_modeI.set_header(["JpsiGammaLambdaLambdabarAlg/JpsiGammaLambdaLambdabar.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:lambda_mass_window, "Λ and Λ̄ candidates must satisfy |M(pπ) - M_Λ| < 5 MeV/c²; the by_minimizing_mass_difference step picks the combination closest to the nominal Λ mass and the residual 5 MeV/c² acceptance window is imposed as an event-level cut")

sel_modeI = Selection.new
  .select_track {
    cos_theta  0.93        # |cosθ| < 0.93
    Vz         20.0        # |Vz| < 20 cm
    Vr         10.0        # Vr < 10 cm
    nChrp      "==2"       # exactly two positive tracks
    nChrn      "==2"       # exactly two negative tracks
    nNet       "==0"       # net charge zero
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025    # 25 MeV (barrel)
    energyThreshold_e 0.050    # 50 MeV (endcap)
    nGam              ">=1"    # Mode I: at least one photon
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p+ / p̄ vs K, π
    identify :pion,   against: [:kaon, :proton] # π+ / π- vs K, p
    nprp ">=1"
    nprm ">=1"
    npip ">=1"
    npim ">=1"
  }
  .secondary_vertex_fit([:prp, :pim]) {          # Λ → p π⁻
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {          # Λ̄ → p̄ π⁺
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :Lambda, :Lambda_bar]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
    # veto photons that form a Σ0 with Λ or Λ̄ within ±20 MeV/c²
    invariant_mass_of(:gamma, :Lambda).out_of(1.1726, 1.2126)
    invariant_mass_of(:gamma, :Lambda_bar).out_of(1.1726, 1.2126)
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

# ==================================================================
# Mode II event selection
# ==================================================================
alg_modeII = Algorithm.new("JpsiGGammaLambdaLambdabar")
alg_modeII.set_header(["JpsiGGammaLambdaLambdabarAlg/JpsiGGammaLambdaLambdabar.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:lambda_mass_window, "Λ and Λ̄ candidates must satisfy |M(pπ) - M_Λ| < 5 MeV/c²; the by_minimizing_mass_difference step picks the combination closest to the nominal Λ mass and the residual 5 MeV/c² acceptance window is imposed as an event-level cut")
          .note(:sigma0_selection, "Mode II requires the LOWER-energy photon to form |M(γ_low Λ) - M_Σ0| < 10 MeV/c²; the DSL cannot single out the lower-energy photon, so the Σ0 mass window is applied over the γΛ combination and the lower-energy assignment is enforced at the analysis level")

sel_modeII = Selection.new
  .select_track {
    cos_theta  0.93
    Vz         20.0
    Vr         10.0
    nChrp      "==2"
    nChrn      "==2"
    nNet       "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025    # 25 MeV (barrel)
    energyThreshold_e 0.050    # 50 MeV (endcap)
    nGam              ">=2"    # Mode II: at least two photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion,   against: [:kaon, :proton]
    nprp ">=1"
    nprm ">=1"
    npip ">=1"
    npim ">=1"
  }
  .secondary_vertex_fit([:prp, :pim]) {          # Λ → p π⁻
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {          # Λ̄ → p̄ π⁺
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :gamma, :Lambda, :Lambda_bar]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
    # require the γΛ system to be consistent with a Σ0 (±10 MeV/c²)
    invariant_mass_of(:gamma, :Lambda).within(1.1826, 1.2026)
    # veto γγ combinations consistent with a π0 (±30 MeV/c²)
    invariant_mass_of(:gamma, :gamma).out_of(0.105, 0.165)
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

# ==================================================================
# Execute both algorithms
# ==================================================================
root_files_modeI  = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])
root_files_modeII = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])