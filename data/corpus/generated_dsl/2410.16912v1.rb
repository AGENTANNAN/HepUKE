# =====================================================================
# Λ_c⁺ → Λ K_S⁰ K⁺ / Λ K_S⁰ π⁺ / Λ K*⁺ (K*⁺ → K_S⁰ π⁺)
# with Λ → p π⁻  and  K_S⁰ → π⁺ π⁻
# Real data + inclusive MC at the seven energy points 4.600 – 4.700 GeV
# (703-4600, 706-4610, 706-4620, 706-4640, 706-4660, 706-4680, 706-4700; ~4.5 fb⁻¹)
# =====================================================================

### ------------------------- Datasets ------------------------- ###
data_4600 = DatasetManager.real_data.find("703_4600")   # 4.600 GeV
data_4610 = DatasetManager.real_data.find("706_4610")   # 4.611 GeV
data_4620 = DatasetManager.real_data.find("706_4620")   # 4.628 GeV
data_4640 = DatasetManager.real_data.find("706_4640")   # 4.641 GeV
data_4660 = DatasetManager.real_data.find("706_4660")   # 4.661 GeV
data_4680 = DatasetManager.real_data.find("706_4680")   # 4.682 GeV
data_4700 = DatasetManager.real_data.find("706_4700")   # 4.699 GeV
data_points = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_points = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

### ------------------------ Decay cards ----------------------- ###
# Mode 1: Λ_c⁺ → Λ K_S⁰ K⁺
decay_card_LcToLambdaKsK = <<~DECAYCARD
    Decay Lambda_c+
    1.000 Lambda0 K_S0 K+      PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi-               PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi-              PHSP;
    Enddecay

    End
DECAYCARD

# Mode 2: Λ_c⁺ → Λ K_S⁰ π⁺ (non-resonant)
decay_card_LcToLambdaKsPi = <<~DECAYCARD
    Decay Lambda_c+
    1.000 Lambda0 K_S0 pi+     PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi-               PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi-              PHSP;
    Enddecay

    End
DECAYCARD

# Mode 3: Λ_c⁺ → Λ K*⁺, K*⁺ → K_S⁰ π⁺  (resonance described in the decay card)
decay_card_LcToLambdaKstar = <<~DECAYCARD
    Decay Lambda_c+
    1.000 Lambda0 K*+          PHSP;
    Enddecay

    Decay K*+
    1.000 K_S0 pi+             PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi-               PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi-              PHSP;
    Enddecay

    End
DECAYCARD

### -------------------- Exclusive MC samples ------------------- ###
# 100k events per decay mode, generated for every energy point (energy scan)
exMC_K = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "exmc_LcToLambdaKsK"
    config.events        = 100000
    config.decay_card    = decay_card_LcToLambdaKsK
    config.cross_section = :default
end

exMC_pi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "exmc_LcToLambdaKsPi"
    config.events        = 100000
    config.decay_card    = decay_card_LcToLambdaKsPi
    config.cross_section = :default
end

exMC_Kstar = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "exmc_LcToLambdaKstar"
    config.events        = 100000
    config.decay_card    = decay_card_LcToLambdaKstar
    config.cross_section = :default
end

### ================================================================ ###
###  Algorithm 1 : Λ_c⁺ → Λ K_S⁰ K⁺                                 ###
### ================================================================ ###
alg_name_K = "LcToLambdaKsK"
alg_K = Algorithm.new(alg_name_K)
alg_K.set_header(["#{alg_name_K}Alg/#{alg_name_K}.h"])
     .set_constant({ "ECMS" => [:double, 4.65] })   # representative CMS energy of the scan
     .note(:secondary_vertex_selection,
           "Λ (from p π⁻) and K_S⁰ (from π⁺ π⁻) are reconstructed by secondary-vertex fits " \
           "constrained to the nominal Λ and K_S⁰ masses; require vertex-fit χ² < 100 and at " \
           "least one candidate for each; mass windows M(Λ) ∈ [1.090, 1.140] GeV/c² and " \
           "M(K_S⁰) ∈ [0.450, 0.540] GeV/c²; flight significance > 2σ.")
     .note(:delta_e_selection,
           "After the 4C kinematic fit (χ² < 200), if several candidate combinations survive, " \
           "keep the one with the smallest |ΔE| and require ΔE ∈ [−0.02, 0.02] GeV.")

sel_K = Selection.new
sel_K.select_track {                       # charged track selection
        cos_theta 0.93                     # |cosθ| < 0.93
        Vz 10.0                            # |Vz| < 10 cm
        Vr 1.0                             # Vr < 1 cm
        nChrp ">=3"                        # at least 3 positive tracks
        nChrn ">=2"                        # at least 2 negative tracks
        nNet  "==1"                        # net charge +1
     }
     .select_photon {                      # photon selection
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0                # > 10° to nearest charged track
        energyThreshold_b 0.025            # 25 MeV (barrel)
        energyThreshold_e 0.050            # 50 MeV (endcap)
        nGam "<4"                          # fewer than four photons
     }
     .pid(method: :probability) {          # PID : probability method
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]   # p vs K, π
        identify :kaon,   against: [:pion]          # K vs π  (K⁺ mode)
        nprp ">=1"                                  # at least one proton
        nkp  ">=1"                                  # at least one kaon
     }
     .remove([:prp <= :chrgp, :kp <= :chrgp])       # remove p and K from positive-track list
     .assign({:chrgp => :pip, :chrgn => :pim})      # remaining + tracks → π⁺, − tracks → π⁻
     .secondary_vertex_fit([:prp, :pim]) {          # Λ → p π⁻
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
     }
     .secondary_vertex_fit([:pip, :pim]) {          # K_S⁰ → π⁺ π⁻
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
     }
     .kinematic_fit([:Lambda, :K_S0, :kp]) {        # 4C fit to Λ K_S⁰ K⁺
        nominal
        constrain_four_momentum
        chi2_cut 200
     }

alg_K.with_decay_card(decay_card_LcToLambdaKsK).apply(sel_K)

### ================================================================ ###
###  Algorithm 2 : Λ_c⁺ → Λ K_S⁰ π⁺                                 ###
### ================================================================ ###
alg_name_pi = "LcToLambdaKsPi"
alg_pi = Algorithm.new(alg_name_pi)
alg_pi.set_header(["#{alg_name_pi}Alg/#{alg_name_pi}.h"])
      .set_constant({ "ECMS" => [:double, 4.65] })
      .note(:secondary_vertex_selection,
            "Λ (from p π⁻) and K_S⁰ (from π⁺ π⁻) are reconstructed by secondary-vertex fits " \
            "constrained to the nominal Λ and K_S⁰ masses; require vertex-fit χ² < 100 and at " \
            "least one candidate for each; mass windows M(Λ) ∈ [1.090, 1.140] GeV/c² and " \
            "M(K_S⁰) ∈ [0.450, 0.540] GeV/c²; flight significance > 2σ.")
      .note(:delta_e_selection,
            "After the 4C kinematic fit (χ² < 200), if several candidate combinations survive, " \
            "keep the one with the smallest |ΔE| and require ΔE ∈ [−0.02, 0.02] GeV.")

sel_pi = Selection.new
sel_pi.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=3"
        nChrn ">=2"
        nNet  "==1"
     }
     .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam "<4"
     }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]   # p vs K, π
        identify :pion,   against: [:kaon]          # π vs K  (π⁺ mode)
        nprp ">=1"                                  # at least one proton
     }
     .remove([:prp <= :chrgp])                      # remove protons from positive-track list
     .assign({:chrgp => :pip, :chrgn => :pim})      # remaining + tracks → π⁺, − tracks → π⁻
     .secondary_vertex_fit([:prp, :pim]) {          # Λ → p π⁻
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
     }
     .secondary_vertex_fit([:pip, :pim]) {          # K_S⁰ → π⁺ π⁻
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
     }
     .kinematic_fit([:Lambda, :K_S0, :pip]) {       # 4C fit to Λ K_S⁰ π⁺
        nominal
        constrain_four_momentum
        chi2_cut 200
     }

alg_pi.with_decay_card(decay_card_LcToLambdaKsPi).apply(sel_pi)

### ================================================================ ###
###  Algorithm 3 : Λ_c⁺ → Λ K*⁺, K*⁺ → K_S⁰ π⁺                     ###
###  Same event-selection chain as the π⁺ mode; the K*⁺ resonance    ###
###  is described only in the decay card.                           ###
### ================================================================ ###
alg_name_Kstar = "LcToLambdaKstar"
alg_Kstar = Algorithm.new(alg_name_Kstar)
alg_Kstar.set_header(["#{alg_name_Kstar}Alg/#{alg_name_Kstar}.h"])
        .set_constant({ "ECMS" => [:double, 4.65] })
        .note(:secondary_vertex_selection,
              "Λ (from p π⁻) and K_S⁰ (from π⁺ π⁻) are reconstructed by secondary-vertex fits " \
              "constrained to the nominal Λ and K_S⁰ masses; require vertex-fit χ² < 100 and at " \
              "least one candidate for each; mass windows M(Λ) ∈ [1.090, 1.140] GeV/c² and " \
              "M(K_S⁰) ∈ [0.450, 0.540] GeV/c²; flight significance > 2σ.")
        .note(:delta_e_selection,
              "After the 4C kinematic fit (χ² < 200), if several candidate combinations survive, " \
              "keep the one with the smallest |ΔE| and require ΔE ∈ [−0.02, 0.02] GeV.")

sel_Kstar = sel_pi.dup    # identical selection chain to the Λ K_S⁰ π⁺ mode

alg_Kstar.with_decay_card(decay_card_LcToLambdaKstar).apply(sel_Kstar)

### --------------------------- Execution --------------------------- ###
root_files_K     = alg_K.execute_on(data_points + incMC_points + exMC_K)
root_files_pi    = alg_pi.execute_on(data_points + incMC_points + exMC_pi)
root_files_Kstar = alg_Kstar.execute_on(data_points + incMC_points + exMC_Kstar)