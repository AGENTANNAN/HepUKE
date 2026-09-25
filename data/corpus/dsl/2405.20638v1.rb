# BESIII DSL: 2405.20638v1 — Study of χ_cJ → ΛΛ̄φ
# ψ(3686) data, (2712.4±14.3)×10^6 events, BOSS 709
# ψ(3686) → γ χ_cJ, χ_cJ → Λ Λ̄ φ, Λ→pπ⁻, Λ̄→p̄π⁺, φ→K⁺K⁻
# Three χ_cJ states (J=0,1,2) share identical event selection

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# =============================================
# Decay cards for χ_c0, χ_c1, χ_c2 signal MC
# ψ(3686) → γ χ_cJ, χ_cJ → Λ Λ̄ φ
# BODY3 generator used for χ_cJ→ΛΛ̄φ (data-driven Dalitz-plot model)
# =============================================

decay_card_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.000 Lambda0 anti-Lambda0 phi PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay
    End
DECAYCARD

decay_card_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.000 Lambda0 anti-Lambda0 phi PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay
    End
DECAYCARD

decay_card_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 PHSP;
    Enddecay

    Decay chi_c2
    1.000 Lambda0 anti-Lambda0 phi PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay
    End
DECAYCARD

# Exclusive MC for each χ_cJ state
exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "chic0_LambdaLambdabar_phi"
  config.related_dataset = psip_data
  config.events = 200000
  config.decay_card = decay_card_chic0
  config.cross_section = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "chic1_LambdaLambdabar_phi"
  config.related_dataset = psip_data
  config.events = 200000
  config.decay_card = decay_card_chic1
  config.cross_section = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "chic2_LambdaLambdabar_phi"
  config.related_dataset = psip_data
  config.events = 200000
  config.decay_card = decay_card_chic2
  config.cross_section = :default
end

# =============================================
# Algorithm — single algorithm for all three χ_cJ states
# (identical event selection; signal yields extracted via ROOT-level simultaneous fit)
# =============================================
alg = Algorithm.new("ChicJToLambdaLambdabarPhi")
alg.set_header(["ChicJToLambdaLambdabarPhiAlg/ChicJToLambdaLambdabarPhi.h"])
   .set_constant({ "ECMS" => [:double, 3.686] })
   .note(:pid_method, "Paper uses highest-CL PID among pion/kaon/proton hypotheses. DSL approximates with probability PID (prob_cut 0.001) identifying protons and kaons against other hadrons. Some tracks may be misidentified; this is a BOSS-level approximation.")
   .note(:lambda_decay_length, "Decay length of Λ/Λ̄ required > 0; handled implicitly by secondary_vertex_fit which reconstructs displaced vertices.")
   .note(:best_lambda_combination, "Best ΛΛ̄ pair selected by minimizing ΔM = sqrt((M(pπ⁻)-m(Λ))² + (M(p̄π⁺)-m(Λ̄))²). DSL approximates with by_minimizing_mass_difference in each secondary vertex fit.")
   .note(:phi_mass_window, "φ signal region: |M(K⁺K⁻) - m_φ| < 18 MeV/c²; sideband: 1055–1127 MeV/c². Applied at ROOT level.")
   .note(:omega_veto, "Ω⁻/Ω̄⁺ veto: |M(ΛK⁻) - m_Ω⁻| > 12 MeV/c² and |M(Λ̄K⁺) - m_Ω̄⁺| > 12 MeV/c². Applied at ROOT level.")
   .note(:lambda_mass_window, "Λ/Λ̄ signal mass window ±6 MeV/c² around nominal Λ mass; 2D sidebands used for non-ΛΛ̄ background estimation. ROOT-level.")
   .note(:phi_sideband, "φ sideband region [1055, 1127] MeV/c² used for simultaneous fit with normalization factor f_φ=0.71. ROOT-level.")
   .note(:body3_model, "BODY3 generator used for χ_cJ→ΛΛ̄φ to model intermediate structures; data-driven Dalitz-plot reweighting applied at ROOT level.")
   .note(:four_c_chi2_optimization, "χ²_4C optimized to < 60 via FOM (ε/(a/2+√N_bkg) with a=3). Paper-specific tight cut replaces default chi2_cut 200.")

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
       nGam ">=1"
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       identify :kaon, against: [:pion, :proton]
       nprp ">=1"
       nprm ">=1"
       nkp ">=1"
       nkm ">=1"
     }
     .remove([:prp <= :chrgp, :prm <= :chrgn])
     .remove([:kp <= :chrgp, :km <= :chrgn])
     # Remaining charged tracks are pions (Λ decay daughters)
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
     # 4C kinematic fit: ψ(3686) → γ Λ Λ̄ K⁺ K⁻
     # Four-momentum conservation; no mass constraints needed
     .kinematic_fit([:gamma, :Lambda, :Lambda_bar, :kp, :km]) {
       nominal
       constrain_four_momentum
       chi2_cut 60
     }

alg.with_decay_card(decay_card_chic0).apply(sel)
alg.execute_on([psip_data, psip_incMC, exMC_chic0, exMC_chic1, exMC_chic2])