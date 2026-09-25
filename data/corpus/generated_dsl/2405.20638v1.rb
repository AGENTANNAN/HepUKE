### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data, BOSS 709 @ 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # Matching inclusive MC sample

# Decay cards — ψ(3686) → γ χ_cJ, χ_cJ → Λ Λ̄ φ, one card per χ_cJ state.
# All three χ_cJ share the identical final state γ Λ Λ̄ K⁺K⁻ (through Λ→pπ⁻, Λ̄→p̄π⁺, φ→K⁺K⁻).
decay_card_chi_c0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.0000 Lambda0 anti-Lambda0 phi PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

decay_card_chi_c1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 Lambda0 anti-Lambda0 phi PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

decay_card_chi_c2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2 PHSP;
    Enddecay

    Decay chi_c2
    1.0000 Lambda0 anti-Lambda0 phi PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Exclusive MC, 200k events for each χ_cJ mode
exMC_chi_c0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chi_c0_LLbarphi"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_chi_c0
  config.cross_section   = :default
end

exMC_chi_c1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chi_c1_LLbarphi"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_chi_c1
  config.cross_section   = :default
end

exMC_chi_c2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chi_c2_LLbarphi"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_chi_c2
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# The three χ_cJ states share identical final states and identical selection criteria,
# so a single Algorithm instance serves all of them.
alg_name   = "ChiCJGammaLLbarPhi"
chi_cJ_alg = Algorithm.new(alg_name)
chi_cJ_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})       # √s = 3.686 GeV
          .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {
      cos_theta 0.93     # |cosθ| < 0.93
      Vz        10.0     # |Vz| < 10 cm
      Vr        1.0      # Vr < 1 cm (transverse plane)
      nChrp     ">=3"    # at least three positive tracks
      nChrn     ">=3"    # at least three negative tracks
      nNet      "==0"    # net charge zero
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      energyThreshold_b 0.025   # 25 MeV barrel threshold
      energyThreshold_e 0.050   # 50 MeV endcap threshold
      angle_to_track    10.0    # angle to any charged track > 10°
      nGam              ">=1"   # at least one photon
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]    # p+ / anti-p- vs K and π
      identify :kaon,   against: [:pion, :proton]  # K+ / K- vs π and p
      nprp ">=1"    # at least one p+
      nprm ">=1"    # at least one anti-p-
      nkp  ">=1"    # at least one K+
      nkm  ">=1"    # at least one K-
  }
  # Remove identified (anti-)protons and kaons; remaining tracks are taken as pions
  .remove([:prp <= :chrgp, :kp <= :chrgp, :prm <= :chrgn, :km <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  # Reconstruct Λ → pπ⁻ from a displaced secondary vertex
  .secondary_vertex_fit([:prp, :pim]) {
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  # Reconstruct anti-Λ → anti-p- π+ from a displaced secondary vertex
  .secondary_vertex_fit([:prm, :pip]) {
      build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  # Final 4C kinematic fit to γ Λ Λ̄ K⁺K⁻
  .kinematic_fit([:gamma, :Lambda, :Lambda_bar, :kp, :km]) {
      nominal                   # nominal fit — corrected four-momenta are used
      constrain_four_momentum   # 4-momentum conservation against the CMS
      chi2_cut 60               # χ² < 60
  }

# Single shared algorithm for all three χ_cJ; χ_c1 card used to generate the
# kinematic variables (identical final-state particles for χ_c0, χ_c1, χ_c2).
chi_cJ_alg.with_decay_card(decay_card_chi_c1).apply(event_selection)
root_files = chi_cJ_alg.execute_on([psip_data, psip_incMC, exMC_chi_c0, exMC_chi_c1, exMC_chi_c2])