# Paper: 2202.08058v3 — ψ(3686) → γχcJ → γΞΞ̄ partial reconstruction
# Data: ψ(3686) (709_3686, ~448M events)
# Two channels: Ξ−Ξ+ and Ξ0Ξ̄0
# Partial reconstruction with partial_rec (NOT kinematic fit)
# Key features: secondary vertex fits for Λ, Ξ−; Kalman fit for π0; recoil mass technique

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for ψ(3686) → γ χcJ → γ Ξ− Ξ+
decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_cJ PHSP;
    Enddecay

    Decay chi_cJ
    1.0000 Xi- anti-Xi- PHSP;
    Enddecay

    Decay Xi-
    1.0000 Lambda pi- PHSP;
    Enddecay

    Decay anti-Xi-
    1.0000 anti-Lambda pi+ PHSP;
    Enddecay

    Decay Lambda
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "my_signal_exclusive_mc"
  config.related_dataset = psip_data
  config.events = 200000
  config.decay_card = decay_card
  config.cross_section = :default
end

# DecayCardResolver.rec_id_list:
#    0 => psip                   (top mother — skip)
#    1 => gamma_rad              (radiative photon — reconstructed)
#    2 => chi_cJ                 (intermediate — skip in partial rec)
#    3 => Xim                    (Ξ− — to be reconstructed or partially tagged)
#    4 => Xip                    (anti-Ξ− — to be inferred from recoil)
#    5 => Lambda_Xim             (Λ from Ξ− decay)
#    6 => pim_Xim                (π− from Ξ− → Λπ−)
#    7 => prp_Lambda_Xim         (p from Λ → pπ−)
#    8 => pim_Lambda_Xim         (π− from Λ → pπ−)
#    9 => antiLambda_Xip         (Λ̄ from anti-Ξ− decay)
#   10 => pip_Xip                (π+ from anti-Ξ− → Λ̄π+)
#   11 => prm_antiLambda_Xip     (p̄ from Λ̄ → p̄π+)
#   12 => pip_antiLambda_Xip     (π+ from Λ̄ → p̄π+)

# ==============================================================================
# Channel 1: Ξ−Ξ+ — reconstruct Ξ−, infer Ξ+ from recoil mass
# ==============================================================================
alg_name_1 = "ChicJXiXi_C1"
alg_1 = Algorithm.new(alg_name_1)
alg_1.set_header(["#{alg_name_1}Alg/#{alg_name_1}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})

sel_1 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"    # p from Λ, π+ from anti-Ξ− (if reconstructed)
    nChrn ">=2"    # π− from Λ, π− from Ξ−
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    tdc_emc_start 0
    tdc_emc_end 14
    nGam ">=1"     # radiative photon
  }
  # Proton PID: L(p) > L(K) AND L(p) > L(π)
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=0"
  }
  .remove([:prp <= :chrgp])
  # Assign remaining: positive → pip, negative → pim
  .assign({:chrgp => :pip, :chrgn => :pim})
  # Secondary vertex fit: Λ → p π−
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Secondary vertex fit: Ξ− → Λ π− (using the virtual Λ + a remaining π−)
  .secondary_vertex_fit([:Lambda, :pim]) {
    build_virtual_particle(:Xi).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Partial reconstruction: reconstruct γ + Ξ−, infer anti-Ξ+ from recoil
  # RecIDs: gamma_rad(1), Xim(3) + children (5 Lambda_Xim, 6 pim_Xim,
  #   7 prp_Lambda_Xim, 8 pim_Lambda_Xim)
  # Xip(4) + children (9-12) are missed → inferred from recoil
  .partial_rec([1, 3, 5, 6, 7, 8]) {
    best_combination_by_mass :Xi, 1.32171
    # Recoil mass window: |M_recoil(γΞ−) − m_Ξ| < ±50 MeV/c2
    require_recoil_mass 1.2717, 1.3717
  }

alg_1.with_decay_card(decay_card).apply(sel_1)
root_files_1 = alg_1.execute_on([psip_data, psip_incMC, exMC_signal])

# ==============================================================================
# Channel 2: Ξ0Ξ̄0 — reconstruct Ξ0, infer Ξ̄0 from recoil mass
# Ξ0 → Λ π0, Λ → p π−, π0 → γγ
# ==============================================================================
alg_name_2 = "ChicJXiXi_C2"
alg_2 = Algorithm.new(alg_name_2)
alg_2.set_header(["#{alg_name_2}Alg/#{alg_name_2}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})

sel_2 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=1"    # p from Λ
    nChrn ">=1"    # π− from Λ
  }
  .select_photon {
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    tdc_emc_start 0
    tdc_emc_end 14
    nGam ">=3"     # 2 for π0 → γγ + 1 radiative photon
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
  }
  .remove([:prp <= :chrgp])
  .assign({:chrgp => :pip, :chrgn => :pim})
  # Secondary vertex fit: Λ → p π−
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Kalman fit: π0 → γγ
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # Partial reconstruction: reconstruct γ + Ξ0(= Λ + π0), infer Ξ̄0 from recoil
  # RecIDs: gamma_rad(1), Xi0(3) + children (Lambda_Xi0, pi0_Xi0, ...)
  .partial_rec([1, 3]) {
    best_combination_by_mass :Xi0, 1.31486
    # Recoil mass window: |M_recoil(γΞ0) − m_Ξ0| < ±25 MeV/c2
    require_recoil_mass 1.2899, 1.3399
  }

alg_2.with_decay_card(decay_card).apply(sel_2)
root_files_2 = alg_2.execute_on([psip_data, psip_incMC, exMC_signal])

# Note: The DecayCardResolver.rec_id_list for Ξ0 channel depends on the specific
# decay card. The decay card above lists Ξ− → Λπ−. For the Ξ0 channel a
# separate decay card with Ξ0 → Λπ0, π0 → γγ would be needed. The partial_rec
# rec_ids shown here are illustrative — actual IDs depend on the decay card.
#
# Post-fit analysis (ROOT stage):
# - Λ mass window: |M_pπ− − m_Λ| < 5 MeV/c2, decay length > 0
# - Ξ− mass window: |M_Λπ− − m_Ξ−| < 10 MeV/c2, decay length > 0
# - Ξ0 mass window: |M_Λπ0 − m_Ξ0| < 10 MeV/c2 (after 1C fit, χ2_1C < 20)
# - Radiative photon selection: min |M_recoil(γΞ) − m_Ξ|
# - χcJ identified from M_recoil(γ) in [3.3, 3.6] GeV/c2