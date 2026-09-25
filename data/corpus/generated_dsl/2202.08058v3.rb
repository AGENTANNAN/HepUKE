# ============================================================================
# ψ(3686) → γ χ_cJ → γ Ξ Ξ̄  via partial reconstruction (π⁺π⁻-independent)
#   channel 1 : Ξ⁻ Ξ⁺        rec: γ + Ξ⁻ (→ Λπ⁻, Λ → pπ⁻); Ξ⁺ inferred from recoil
#   channel 2 : Ξ⁰ Ξ̄⁰        rec: γ + Ξ⁰ (→ Λπ⁰, Λ → pπ⁻, π⁰ → γγ); Ξ̄⁰ from recoil
#   Partial reconstruction replaces the kinematic fit entirely.
# ============================================================================

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data (~448M events)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # ψ(3686) inclusive MC

# --- Decay card, charged channel : χ_cJ → Ξ⁻ Ξ̄⁺ (= Ξ⁻Ξ⁺) ---
decay_card_charged = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 Xi- anti-Xi+ PHSP;
    Enddecay

    Decay Xi-
    1.000 Lambda0 pi- HypWK;
    Enddecay

    Decay anti-Xi+
    1.000 anti-Lambda0 pi+ HypWK;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# --- Decay card, neutral channel : χ_cJ → Ξ⁰ Ξ̄⁰ ---
decay_card_neutral = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 Xi0 anti-Xi0 PHSP;
    Enddecay

    Decay Xi0
    1.000 Lambda0 pi0 HypWK;
    Enddecay

    Decay anti-Xi0
    1.000 anti-Lambda0 pi0 HypWK;
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

# --- Exclusive signal MC (200k events per channel) ---
exMC_charged = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chic_xim_xip"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_charged
  config.cross_section   = :default
end

exMC_neutral = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chic_xi0_xi0bar"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_neutral
  config.cross_section   = :default
end

### Event selection — channel 1 : Ξ⁻Ξ⁺ ###
alg_name_chrg = "GamChicXimXip"
alg_chrg = Algorithm.new(alg_name_chrg)
alg_chrg.set_header(["#{alg_name_chrg}Alg/#{alg_name_chrg}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .set_alias({"std::vector<double>" => "Vdouble"})

# RecID map (decay-card traversal): 0 psi(2S), 1 gamma, 2 chi_c1, 3 Xi-,
# 4 anti-Xi+, 5 Lambda0, 6 pi-(Xi-), 7 p+(Lambda), 8 pi-(Lambda), 9 anti-Lambda0,
# 10 pi+(anti-Xi+), 11 anti-p-, 12 pi+(anti-Lambda0)
sel_chrg = Selection.new
sel_chrg.select_track {
            cos_theta 0.93      # |cosθ| < 0.93
            Vz        10.0      # |Vz| < 10 cm
            Vr        1.0       # Vr < 1 cm
            nChrp     ">=2"     # at least two positive tracks
            nChrn     ">=2"     # at least two negative tracks
          }
        .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    10.0    # ≥10° separation from charged tracks
            energyThreshold_b 0.025   # 25 MeV (barrel)
            energyThreshold_e 0.050   # 50 MeV (endcap)
            nGam              ">=1"   # at least one radiative photon
          }
        .pid(method: :probability) {
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]   # p vs K/π
            nprp ">=1"                                   # ≥1 proton
          }
        .remove([:prp <= :chrgp])                        # remove proton from positive-track list
        .assign({:chrgp => :pip, :chrgn => :pim})        # remaining tracks → π⁺ / π⁻
        .secondary_vertex_fit([:prp, :pim]) {            # Λ → pπ⁻
            build_virtual_particle(:Lambda).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
        .secondary_vertex_fit([:Lambda, :pim]) {         # Ξ⁻ → Λπ⁻
            build_virtual_particle(:Xim).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
        .partial_rec([1, 3, 5, 6, 7, 8]) {               # reconstruct γ + Ξ⁻ (Λ pre-built)
            best_combination_by_mass :Xim, 1.32171       # best Ξ⁻ mass
            require_recoil_mass 1.2717, 1.3717           # m_Ξ ± 50 MeV (Ξ⁺ from recoil)
          }

alg_chrg.with_decay_card(decay_card_charged).apply(sel_chrg)
alg_chrg.execute_on([psip_data, psip_incMC, exMC_charged])

### Event selection — channel 2 : Ξ⁰Ξ̄⁰ ###
alg_name_neut = "GamChicXi0Xi0bar"
alg_neut = Algorithm.new(alg_name_neut)
alg_neut.set_header(["#{alg_name_neut}Alg/#{alg_name_neut}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .set_alias({"std::vector<double>" => "Vdouble"})

# RecID map: 0 psi(2S), 1 gamma, 2 chi_c1, 3 Xi0, 4 anti-Xi0, 5 Lambda0, 6 pi0(Xi0),
# 7 p+(Lambda), 8 pi-(Lambda), 9 gamma(pi0), 10 gamma(pi0), 11 anti-Lambda0,
# 12 pi0(anti-Xi0), 13 anti-p-, 14 pi+, 15 gamma, 16 gamma
sel_neut = Selection.new
sel_neut.select_track {
            cos_theta 0.93
            Vz        10.0
            Vr        1.0
            nChrp     ">=1"     # at least one positive track
            nChrn     ">=1"     # at least one negative track
          }
        .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam              ">=3"   # at least three photons
          }
        .pid(method: :probability) {
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]   # same proton PID
            nprp ">=1"
          }
        .remove([:prp <= :chrgp])
        .assign({:chrgp => :pip, :chrgn => :pim})
        .secondary_vertex_fit([:prp, :pim]) {            # same Λ → pπ⁻ vertex fit
            build_virtual_particle(:Lambda).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
        .kalman_kinematic_fit([:gamma, :gamma]) {        # π⁰ → γγ (1-C mass constraint)
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"                                   # ≥1 π⁰
          }
        .partial_rec([1, 3, 5, 6, 7, 8, 9, 10]) {        # reconstruct γ + Ξ⁰ (Λ, π⁰ pre-built)
            best_combination_by_mass :Xi0, 1.31486       # best Ξ⁰ mass
            require_recoil_mass 1.2899, 1.3399           # m_Ξ⁰ ± 25 MeV (Ξ̄⁰ from recoil)
          }

alg_neut.with_decay_card(decay_card_neutral).apply(sel_neut)
alg_neut.execute_on([psip_data, psip_incMC, exMC_neutral])