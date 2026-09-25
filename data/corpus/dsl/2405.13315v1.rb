# BESIII DSL: 2405.13315v1 — Observation of χ_cJ → ΛΛ̄ω
# ψ(3686) data, (27.12±0.14)×10^8 events, BOSS 709
# ψ(3686) → γ χ_cJ, χ_cJ → Λ Λ̄ ω, Λ→pπ⁻, Λ̄→p̄π⁺, ω→π⁺π⁻π⁰, π⁰→γγ
# Three χ_cJ states (J=0,1,2) share identical event selection

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# =============================================
# Decay cards for χ_c0, χ_c1, χ_c2 signal MC
# All share the same decay chain; only χ_cJ mass differs
# =============================================

# χ_c0 → ΛΛ̄ω
decay_card_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.000 Lambda0 anti-Lambda0 omega PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# χ_c1 → ΛΛ̄ω
decay_card_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.000 Lambda0 anti-Lambda0 omega PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# χ_c2 → ΛΛ̄ω
decay_card_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 PHSP;
    Enddecay

    Decay chi_c2
    1.000 Lambda0 anti-Lambda0 omega PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC for each χ_cJ state
exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "chic0_LambdaLambdabar_omega"
  config.related_dataset = psip_data
  config.events = 200000
  config.decay_card = decay_card_chic0
  config.cross_section = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "chic1_LambdaLambdabar_omega"
  config.related_dataset = psip_data
  config.events = 200000
  config.decay_card = decay_card_chic1
  config.cross_section = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "chic2_LambdaLambdabar_omega"
  config.related_dataset = psip_data
  config.events = 200000
  config.decay_card = decay_card_chic2
  config.cross_section = :default
end

# =============================================
# Algorithm — single algorithm for all three χ_cJ states
# (identical event selection; signal yields extracted via ROOT-level simultaneous fit)
# =============================================
alg = Algorithm.new("ChicJToLambdaLambdabarOmega")
alg.set_header(["ChicJToLambdaLambdabarOmegaAlg/ChicJToLambdaLambdabarOmega.h"])
   .set_constant({ "ECMS" => [:double, 3.686] })
   .note(:pid_method, "Paper uses highest-CL PID among pion/proton/kaon hypotheses. DSL approximates with probability PID (prob_cut 0.001) identifying protons against pion/kaon; remaining tracks assigned as pions.")
   .note(:lambda_decay_length, "Decay length of Λ/Λ̄ required > 0; handled implicitly by secondary_vertex_fit which reconstructs displaced vertices.")
   .note(:best_lambda_combination, "Best ΛΛ̄ pair selected by minimizing ΔM = sqrt((M(pπ⁻)-m(Λ))² + (M(p̄π⁺)-m(Λ̄))²). DSL approximates with by_minimizing_mass_difference in each secondary vertex fit.")
   .note(:competing_jpsi_veto, "χ²_signal < χ²_bkg veto against ψ(3686)→π⁺π⁻J/ψ requires reinterpreting tracks under different PID hypotheses; not expressible in DSL — applied at ROOT level.")
   .note(:mass_vetoes, "Mass vetoes from Table I (Σ*, Ξ, J/ψ, Σ⁰ windows) are applied at ROOT level after kinematic fit selects best combination.")
   .note(:omega_sideband, "ω sideband regions [0.693,0.747] and [0.819,0.873] GeV/c² used for simultaneous fit; ROOT-level.")
   .note(:lambda_mass_window, "Λ/Λ̄ signal mass window ±8 MeV/c² around nominal Λ mass; ω signal window [0.756, 0.810] GeV/c²; ROOT-level.")
   .note(:body3_model, "BODY3 generator used for χ_cJ→ΛΛ̄ω to model intermediate structures; data-driven Dalitz-plot reweighting applied at ROOT level.")
   .note(:five_c_chi2_optimization, "χ²_5C optimized to < 30 via FOM (S/√(S+B) maximization); paper-specific tight cut replaces default chi2_cut 200.")

# =============================================
# Event selection chain
# =============================================
sel = Selection.new
sel.select_track {
       cos_theta 0.93
       Vz 10.0
       Vr 1.0
       nChrp ">=3"
       nChrn ">=3"
       nNet "==0"
     }
     .select_photon {
       tdc_emc_start 0
       tdc_emc_end 14
       angle_to_track 10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam ">=3"
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       nprp ">=1"
       nprm ">=1"
     }
     .remove([:prp <= :chrgp, :prm <= :chrgn])
     # Remaining charged tracks are pions (including those from Λ/Λ̄ decay and ω decay)
     .assign({:chrgp => :pip, :chrgn => :pim})
     # Reconstruct Λ → p π⁻ with secondary vertex fit
     .secondary_vertex_fit([:prp, :pim]) {
       build_virtual_particle(:Lambda).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
     # Reconstruct Λ̄ → p̄ π⁺ with secondary vertex fit
     .secondary_vertex_fit([:prm, :pip]) {
       build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
     # Reconstruct π⁰ → γγ (one π⁰ from ω decay; remaining γ is the radiative photon from ψ(3686)→γ χ_cJ)
     .kalman_kinematic_fit([:gamma, :gamma]) {
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
       chi2_cut 25
       npi0 ">=1"
     }
     # Nominal 5C kinematic fit: 4C + π⁰ mass constraint
     # Particles: Λ, Λ̄, π⁺(from ω), π⁻(from ω), γ(radiative), π⁰(from ω)
     .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim, :gamma, :pi0]) {
       nominal
       constrain_four_momentum
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
       chi2_cut 30
     }
     # Competing hypothesis: ψ(3686) → ΛΛ̄ω (no radiative photon, 4C + π⁰ mass)
     # χ² stored for ROOT-level veto: χ²_signal < χ²_psip_to_LambdaLambdabar_omega
     .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim, :pi0]) {
       constrain_four_momentum
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
     }

alg.with_decay_card(decay_card_chic0).apply(sel)
alg.execute_on([psip_data, psip_incMC, exMC_chic0, exMC_chic1, exMC_chic2])