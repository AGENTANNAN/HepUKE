# ============================================================
# BOSS DSL — J/psi -> eta' K Kbar pi  (search for the h1(1380))
# Two signal modes sharing the same photon and PID criteria.
# ============================================================

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi 3.097 GeV real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding inclusive MC

# ---- Decay card: Mode I   J/psi -> eta' K+ K- pi0 ----
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.000 eta' K+ K- pi0 PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---- Decay card: Mode II  J/psi -> eta' K+ K- K_S0 ----
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.000 eta' K+ K- K_S0 PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# ---- 100k-event exclusive MC for each mode ----
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_etaprime_KpKm_pi0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_etaprime_KS_KpKm"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# ============================================================
# Mode I : J/psi -> eta' K+ K- pi0
# ============================================================
algI_name = "EtaprimeKKpiModeI"
algI = Algorithm.new(algI_name)
algI.set_header(["#{algI_name}Alg/#{algI_name}.h"])
    .set_constant({"ECMS" => [:double, 3.097]})

selI = Selection.new
selI.select_track {                       # exactly 4 tracks, >=2 positive, >=2 negative
      cos_theta 0.93                      # |cos(theta)| < 0.93
      Vz        10.0                      # |Vz| < 10 cm
      Vr        1.0                       # Vr < 1 cm
      nChrp     ">=2"
      nChrn     ">=2"
      nNet      "==0"
      nTot      "==4"
    }
    .select_photon {                      # >=4 photons (eta -> gamma gamma, pi0 -> gamma gamma)
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=4"
    }
    .pid(method: :probability) {          # pi/K separation
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      identify :pion, against: [:kaon, :proton]
      nkp  ">=1"
      nkm  ">=1"
      npip ">=1"
      npim ">=1"
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {      # pi0 -> gamma gamma  (1-C)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      invariant_mass_of(:gamma, :gamma).within(0.115, 0.155)   # pi0 window 0.02 GeV
      chi2_cut 25
      npi0 ">=1"
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {      # eta -> gamma gamma  (1-C)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      invariant_mass_of(:gamma, :gamma).within(0.518, 0.578)   # eta window 0.03 GeV
      chi2_cut 25
      neta ">=1"
    }
    .kinematic_fit([:pip, :pim, :kp, :km, :pi0, :eta]) {  # 4C to pi+ pi- K+ K- pi0 eta
      nominal
      constrain_four_momentum
      invariant_mass_of(:pip, :pim, :eta).within(0.928, 0.988)   # eta' window 0.03 GeV
      chi2_cut 100
    }

algI.with_decay_card(decay_card_modeI).apply(selI)

# ============================================================
# Mode II : J/psi -> eta' K_S0 K+ K-
# ============================================================
algII_name = "EtaprimeKSKKModeII"
algII = Algorithm.new(algII_name)
algII.set_header(["#{algII_name}Alg/#{algII_name}.h"])
     .set_constant({"ECMS" => [:double, 3.097]})
     .note(:k_s0_vertex_cut,
           "K_S0 -> pi+ pi- secondary-vertex fit required vertex chi2 < 100 and
            |M(pi+pi-) - m_K_S0| < 0.01 GeV; combination closest to the nominal
            K_S0 mass is kept by minimizing the mass difference")

selII = Selection.new
selII.select_track {                      # exactly 6 tracks, >=3 positive, >=3 negative
       cos_theta 0.93                     # |cos(theta)| < 0.93
       Vz        10.0                     # |Vz| < 10 cm
       Vr        1.0                      # Vr < 1 cm
       nChrp     ">=3"
       nChrn     ">=3"
       nNet      "==0"
       nTot      "==6"
     }
     .select_photon {                     # >=2 photons (eta -> gamma gamma)
       tdc_emc_start     0
       tdc_emc_end       14
       angle_to_track    10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam              ">=2"
     }
     .pid(method: :probability) {         # same pi/K PID as Mode I
       prob_cut 0.001
       identify :kaon, against: [:pion, :proton]
       identify :pion, against: [:kaon, :proton]
       nkp  ">=1"
       nkm  ">=1"
       npip ">=1"
       npim ">=1"
     }
     .kalman_kinematic_fit([:gamma, :gamma]) {     # eta -> gamma gamma  (1-C)
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
       invariant_mass_of(:gamma, :gamma).within(0.518, 0.578)   # eta window 0.03 GeV
       chi2_cut 25
       neta ">=1"
     }
     .secondary_vertex_fit([:pip, :pim]) {         # K_S0 -> pi+ pi-
       build_virtual_particle(:K_S0).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
     .kinematic_fit([:pip, :pim, :kp, :km, :K_S0, :eta]) {  # 4C to pi+ pi- K+ K- K_S0 eta
       nominal
       constrain_four_momentum
       invariant_mass_of(:pip, :pim, :eta).within(0.928, 0.988)   # eta' window 0.03 GeV
       chi2_cut 100
     }

algII.with_decay_card(decay_card_modeII).apply(selII)

# ============================================================
# Execute both algorithms
# ============================================================
root_files_I  = algI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])
root_files_II = algII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])