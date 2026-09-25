# ============================================================================
# ψ(2S)(3.686 GeV) → Ω⁻ Ω̄⁺ : absolute branching fractions of Ω⁻ decays
# Single tag: Ω̄⁺ → Λ̄ K⁺  (Λ̄ → p̄ π⁺)
# Signal modes: A (Ω⁻ → Ξ⁰π⁻,  Ξ⁰ → Λπ⁰)
#               B (Ω⁻ → Ξ⁻π⁰,  Ξ⁻ → Λπ⁻)
#               C (Ω⁻ → ΛK⁻)
# ============================================================================

### ---------------------------- Datasets ---------------------------- ###
psip_data  = DatasetManager.real_data.find("709_3686")      # 3.686 GeV real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # matching inclusive MC

### ---------------------------- Decay cards ---------------------------- ###
# Mode A : ψ(2S) → Ω⁻ Ω̄⁺ ; Ω⁻ → Ξ⁰ π⁻ ; Ξ⁰ → Λ π⁰ ; tag Ω̄⁺ → Λ̄ K⁺ ; Λ(Λ̄) → pπ⁻
decay_card_modeA = <<~DECAYCARD
  Decay psi(2S)
  1.000 Omega- anti-Omega+ PHSP;
  Enddecay

  Decay anti-Omega+
  1.000 anti-Lambda0 K+ PHSP;
  Enddecay

  Decay anti-Lambda0
  1.000 anti-p- pi+ HypWK;
  Enddecay

  Decay Omega-
  1.000 Xi0 pi- PHSP;
  Enddecay

  Decay Xi0
  1.000 Lambda0 pi0 PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- HypWK;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode B : ψ(2S) → Ω⁻ Ω̄⁺ ; Ω⁻ → Ξ⁻ π⁰ ; Ξ⁻ → Λ π⁻ ; tag Ω̄⁺ → Λ̄ K⁺ ; Λ(Λ̄) → pπ⁻
decay_card_modeB = <<~DECAYCARD
  Decay psi(2S)
  1.000 Omega- anti-Omega+ PHSP;
  Enddecay

  Decay anti-Omega+
  1.000 anti-Lambda0 K+ PHSP;
  Enddecay

  Decay anti-Lambda0
  1.000 anti-p- pi+ HypWK;
  Enddecay

  Decay Omega-
  1.000 Xi- pi0 PHSP;
  Enddecay

  Decay Xi-
  1.000 Lambda0 pi- PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- HypWK;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode C : ψ(2S) → Ω⁻ Ω̄⁺ ; Ω⁻ → Λ K⁻ ; tag Ω̄⁺ → Λ̄ K⁺ ; Λ(Λ̄) → pπ⁻
decay_card_modeC = <<~DECAYCARD
  Decay psi(2S)
  1.000 Omega- anti-Omega+ PHSP;
  Enddecay

  Decay anti-Omega+
  1.000 anti-Lambda0 K+ PHSP;
  Enddecay

  Decay anti-Lambda0
  1.000 anti-p- pi+ HypWK;
  Enddecay

  Decay Omega-
  1.000 Lambda0 K- PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- HypWK;
  Enddecay

  End
DECAYCARD

### ------------------- Exclusive MC (1M events each) ------------------- ###
exMC_modeA = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_OmegaOmegaBar_modeA"
  config.related_dataset = psip_data
  config.events          = 1000000
  config.decay_card      = decay_card_modeA
  config.cross_section   = :default
end

exMC_modeB = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_OmegaOmegaBar_modeB"
  config.related_dataset = psip_data
  config.events          = 1000000
  config.decay_card      = decay_card_modeB
  config.cross_section   = :default
end

exMC_modeC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_OmegaOmegaBar_modeC"
  config.related_dataset = psip_data
  config.events          = 1000000
  config.decay_card      = decay_card_modeC
  config.cross_section   = :default
end

# ============================================================================
# Mode A :  Ω⁻ → Ξ⁰ π⁻,  Ξ⁰ → Λ π⁰
# ============================================================================
alg_name_A = "OmegaModeA"
alg_A = Algorithm.new(alg_name_A)
alg_A.set_header(["#{alg_name_A}Alg/#{alg_name_A}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})
     .note(:lambda_mass_window,
           "Λ / Λ̄ candidates required to satisfy |M(pπ⁻) − m_Λ| < 11 MeV after the secondary-vertex fit")

sel_A = Selection.new
sel_A.select_track {                       # charged track selection
        cos_theta 0.93                     # |cosθ| < 0.93
        Vz        10.0                     # |Vz| < 10 cm
        Vr        1.0                      # Vr < 1 cm
        nChrp     ">=2"                    # at least 2 positive tracks
        nChrn     ">=3"                    # at least 3 negative tracks
      }
     .select_photon {                      # photon selection (modes A, B)
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025            # barrel E > 25 MeV
        energyThreshold_e 0.050            # endcap  E > 50 MeV
        angle_to_track    10.0             # angle to any charged track > 10 deg
        nGam              ">=2"            # at least two photons
      }
     .pid(method: :probability) {          # PID (probability method)
        prob_cut 0.001
        identify :proton, against: [:pion, :kaon]   # p vs π/K
        identify :kaon,   against: [:pion, :proton] # K vs π/p
        nkp ">=1"                          # at least one K⁺ (tag)
      }
     .remove([:prp <= :chrgp, :prm <= :chrgn, :kp <= :chrgp, :km <= :chrgn])
     .assign({:chrgp => :pip, :chrgn => :pim})     # remaining tracks treated as π⁺/π⁻
     .kalman_kinematic_fit([:gamma, :gamma]) {     # 1C fit → π⁰ from γγ
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)   # M(γγ) window
        chi2_cut 200
        npi0 ">=1"
      }
     .secondary_vertex_fit([:prp, :pim]) {         # signal Λ : p π⁻
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
     .secondary_vertex_fit([:prm, :pip]) {         # tag Λ̄ : p̄ π⁺
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
     .for_each(:pim) {                             # signal π⁻ = highest-energy π⁻
        best { maximize { energy } }
        store :sig_pim_energy
      }
     .kinematic_fit([:Lambda_bar, :kp, :Lambda, :pi0, :pim]) {   # final 4C fit
        nominal
        constrain_four_momentum
        chi2_cut 200
        invariant_mass_of(:Lambda_bar, :kp).within(1.664, 1.680) # Ω tag mass window
      }

alg_A.with_decay_card(decay_card_modeA).apply(sel_A)
root_files_A = alg_A.execute_on([psip_data, psip_incMC, exMC_modeA])

# ============================================================================
# Mode B :  Ω⁻ → Ξ⁻ π⁰,  Ξ⁻ → Λ π⁻
# ============================================================================
alg_name_B = "OmegaModeB"
alg_B = Algorithm.new(alg_name_B)
alg_B.set_header(["#{alg_name_B}Alg/#{alg_name_B}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})
     .note(:lambda_mass_window,
           "Λ / Λ̄ candidates required to satisfy |M(pπ⁻) − m_Λ| < 11 MeV after the secondary-vertex fit")

sel_B = Selection.new
sel_B.select_track {                       # charged track selection
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     ">=2"
        nChrn     ">=3"
      }
     .select_photon {                      # photon selection (modes A, B)
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        angle_to_track    10.0
        nGam              ">=2"
      }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:pion, :kaon]
        identify :kaon,   against: [:pion, :proton]
        nkp ">=1"                          # at least one K⁺ (tag)
      }
     .remove([:prp <= :chrgp, :prm <= :chrgn, :kp <= :chrgp, :km <= :chrgn])
     .assign({:chrgp => :pip, :chrgn => :pim})
     .kalman_kinematic_fit([:gamma, :gamma]) {     # 1C fit → π⁰ from γγ
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)
        chi2_cut 200
        npi0 ">=1"
      }
     .secondary_vertex_fit([:prp, :pim]) {         # signal Λ : p π⁻
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
     .secondary_vertex_fit([:prm, :pip]) {         # tag Λ̄ : p̄ π⁺
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
     .for_each(:pim) {                             # signal π⁻ = highest-energy π⁻
        best { maximize { energy } }
        store :sig_pim_energy
      }
     .kinematic_fit([:Lambda_bar, :kp, :Lambda, :pi0, :pim]) {   # final 4C fit
        nominal
        constrain_four_momentum
        chi2_cut 200
        invariant_mass_of(:Lambda_bar, :kp).within(1.664, 1.680) # Ω tag mass window
      }

alg_B.with_decay_card(decay_card_modeB).apply(sel_B)
root_files_B = alg_B.execute_on([psip_data, psip_incMC, exMC_modeB])

# ============================================================================
# Mode C :  Ω⁻ → Λ K⁻
# ============================================================================
alg_name_C = "OmegaModeC"
alg_C = Algorithm.new(alg_name_C)
alg_C.set_header(["#{alg_name_C}Alg/#{alg_name_C}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})
     .note(:lambda_mass_window,
           "Λ / Λ̄ candidates required to satisfy |M(pπ⁻) − m_Λ| < 11 MeV after the secondary-vertex fit")

sel_C = Selection.new
sel_C.select_track {                       # charged track selection (no photons in mode C)
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     ">=2"
        nChrn     ">=3"
      }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:pion, :kaon]
        identify :kaon,   against: [:pion, :proton]
        nkp ">=1"                          # at least one K⁺ (tag)
        nkm ">=1"                          # at least one K⁻ (signal)
      }
     .remove([:prp <= :chrgp, :prm <= :chrgn, :kp <= :chrgp, :km <= :chrgn])
     .assign({:chrgp => :pip, :chrgn => :pim})
     .secondary_vertex_fit([:prp, :pim]) {         # signal Λ : p π⁻
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
     .secondary_vertex_fit([:prm, :pip]) {         # tag Λ̄ : p̄ π⁺
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
     .kinematic_fit([:Lambda_bar, :kp, :Lambda, :km]) {          # final 4C fit
        nominal
        constrain_four_momentum
        chi2_cut 200
        invariant_mass_of(:Lambda_bar, :kp).within(1.664, 1.680) # Ω tag mass window
      }

alg_C.with_decay_card(decay_card_modeC).apply(sel_C)
root_files_C = alg_C.execute_on([psip_data, psip_incMC, exMC_modeC])