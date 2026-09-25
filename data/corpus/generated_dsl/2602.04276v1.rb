# =====================================================================
# J/psi(3.097 GeV) -> Lambda Lambdabar ; Lambda -> p pi-
# The Lambdabar annihilates at rest on a proton of the beam-pipe
# cooling oil.  Six Lambdabar-p annihilation channels are reconstructed.
# BOSS scope only: datasets, decay cards, exclusive MC and the event
# selection chain up to (and including) the kinematic fit.
# =====================================================================

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # 3.097 GeV J/psi real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding inclusive MC

### Decay cards (EvtGen) — one per signal channel ###
# channel 1 : Lambdabar p -> K_S0 pi+
decay_card_ch1 = <<~DECAYCARD
  Decay J/psi
  1.0 Lambda K_S0 pi+ PHSP;
  Enddecay

  Decay Lambda
  1.0 p+ pi- HypWK;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# channel 2 : Lambdabar p -> K+ pi+ pi-
decay_card_ch2 = <<~DECAYCARD
  Decay J/psi
  1.0 Lambda K+ pi+ pi- PHSP;
  Enddecay

  Decay Lambda
  1.0 p+ pi- HypWK;
  Enddecay

  End
DECAYCARD

# channel 3 : Lambdabar p -> K_S0 pi+ pi+ pi-
decay_card_ch3 = <<~DECAYCARD
  Decay J/psi
  1.0 Lambda K_S0 pi+ pi+ pi- PHSP;
  Enddecay

  Decay Lambda
  1.0 p+ pi- HypWK;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# channel 4 : Lambdabar p -> K+ pi+ pi+ pi- pi-
decay_card_ch4 = <<~DECAYCARD
  Decay J/psi
  1.0 Lambda K+ pi+ pi+ pi- pi- PHSP;
  Enddecay

  Decay Lambda
  1.0 p+ pi- HypWK;
  Enddecay

  End
DECAYCARD

# channel 5 : Lambdabar p -> K+ pi+ pi+ pi+ pi- pi- pi-
decay_card_ch5 = <<~DECAYCARD
  Decay J/psi
  1.0 Lambda K+ pi+ pi+ pi+ pi- pi- pi- PHSP;
  Enddecay

  Decay Lambda
  1.0 p+ pi- HypWK;
  Enddecay

  End
DECAYCARD

# channel 6 : Lambdabar p -> K+ K+ K-
decay_card_ch6 = <<~DECAYCARD
  Decay J/psi
  1.0 Lambda K+ K+ K- PHSP;
  Enddecay

  Decay Lambda
  1.0 p+ pi- HypWK;
  Enddecay

  End
DECAYCARD

### Exclusive MC samples (200k events per channel) ###
exMC_ch1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_LambdabarP_KsPi"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_ch1
  config.cross_section   = :default
end

exMC_ch2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_LambdabarP_KPiPi"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_ch2
  config.cross_section   = :default
end

exMC_ch3 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_LambdabarP_Ks2PiPi"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_ch3
  config.cross_section   = :default
end

exMC_ch4 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_LambdabarP_K2Pi2Pi"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_ch4
  config.cross_section   = :default
end

exMC_ch5 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_LambdabarP_K3Pi3Pi"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_ch5
  config.cross_section   = :default
end

exMC_ch6 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_LambdabarP_2KK"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_ch6
  config.cross_section   = :default
end

# Nominal masses used for the mass windows (GeV/c^2)
# m_Lambda ~ 1.1157, m_K_S0 ~ 0.4976

### ================= Event selection (BOSS) ================= ###

# ---------------------------------------------------------------------
# Channel 1 : Lambdabar p -> K_S0 pi+     (Lambda -> p pi-, K_S0 -> pi+ pi-)
# ---------------------------------------------------------------------
alg_ch1 = Algorithm.new("JpsiLambdabarP_KsPi")
alg_ch1.set_header(["JpsiLambdabarP_KsPiAlg/JpsiLambdabarP_KsPi.h"])
       .set_constant({"ECMS" => [:double, 3.097]})
       .set_alias({"std::vector<double>" => "Vdouble"})
alg_ch1.note(:recoil_mass_window,
             "Lambda recoil mass |M_recoil(Lambda) - m_Lambdabar| < 0.020 GeV/c^2 "
             "(m_Lambdabar = 1.1157 GeV/c^2) required to tag Lambdabar-p annihilation events")
       .note(:annihilation_vertex,
             "Lambda-bar p annihilation vertex required in 3.0 <= R_xy <= 3.5 cm "
             "(annihilation inside the beam-pipe cooling-oil layer)")
       .note(:oil_momentum,
             "|p_oil| < 0.04 GeV/c required to reject Lambdabar annihilations on bound Au/Be/C nucleons")

sel_ch1 = Selection.new
  .select_track {
     cos_theta 0.93     # |cos(theta)| < 0.93
     nChrp  ">=2"       # >= 2 positively charged tracks
     nChrn  ">=1"       # >= 1 negatively charged track
     nTot   ">=3"       # >= 3 charged tracks in total
  }
  .pid(method: :probability) {
     prob_cut 0.001
     identify :proton, against: [:kaon, :pion]  # p / p-bar for Lambda -> p pi-
     identify :pion,   against: [:kaon, :proton]  # pions (Lambda, K_S0 and direct pi+)
     nprp ">=1"
  }
  # Lambda tag: p pi- pair with smallest |M - m_Lambda|; signal window |M - m_Lambda| < 0.003
  .invariant_mass_of(:prp, :pim).within(1.1127, 1.1187)
  .secondary_vertex_fit([:prp, :pim]) {
     build_virtual_particle(:Lambda).by_minimizing_mass_difference
     remove_used_particle_from_candidate_list
  }
  # K_S0 : pi+ pi- pair with smallest |M - m_K_S0|; window 0.015 GeV/c^2
  .invariant_mass_of(:pip, :pim).within(0.4826, 0.5126)
  .secondary_vertex_fit([:pip, :pim]) {
     build_virtual_particle(:K_S0).by_minimizing_mass_difference
     remove_used_particle_from_candidate_list
  }
  # full final state: 4-momentum conservation, Lambda and K_S0 mass constraints
  .kinematic_fit([:Lambda, :K_S0, :pip]) {
     nominal
     constrain_four_momentum
     chi2_cut 200
  }

alg_ch1.with_decay_card(decay_card_ch1).apply(sel_ch1)
alg_ch1.execute_on([jpsi_data, jpsi_incMC, exMC_ch1])

# ---------------------------------------------------------------------
# Channel 2 : Lambdabar p -> K+ pi+ pi-     (Lambda -> p pi-)
# ---------------------------------------------------------------------
alg_ch2 = Algorithm.new("JpsiLambdabarP_KPiPi")
alg_ch2.set_header(["JpsiLambdabarP_KPiPiAlg/JpsiLambdabarP_KPiPi.h"])
       .set_constant({"ECMS" => [:double, 3.097]})
       .set_alias({"std::vector<double>" => "Vdouble"})
alg_ch2.note(:recoil_mass_window,
             "Lambda recoil mass |M_recoil(Lambda) - m_Lambdabar| < 0.020 GeV/c^2 "
             "(m_Lambdabar = 1.1157 GeV/c^2) required to tag Lambdabar-p annihilation events")
       .note(:annihilation_vertex,
             "Lambda-bar p annihilation vertex required in 3.0 <= R_xy <= 3.5 cm "
             "(annihilation inside the beam-pipe cooling-oil layer)")
       .note(:oil_momentum,
             "|p_oil| < 0.04 GeV/c required to reject Lambdabar annihilations on bound Au/Be/C nucleons")

sel_ch2 = Selection.new
  .select_track {
     cos_theta 0.93
     nChrp  ">=2"
     nChrn  ">=2"
     nTot   ">=4"
  }
  .pid(method: :probability) {
     prob_cut 0.001
     identify :proton, against: [:kaon, :pion]
     identify :kaon,   against: [:pion, :proton]   # K+ from the annihilation
     identify :pion,   against: [:kaon, :proton]
     nprp ">=1"
     nkp  ">=1"
  }
  .invariant_mass_of(:prp, :pim).within(1.1127, 1.1187)
  .secondary_vertex_fit([:prp, :pim]) {
     build_virtual_particle(:Lambda).by_minimizing_mass_difference
     remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :kp, :pip, :pim]) {
     nominal
     constrain_four_momentum
     chi2_cut 200
  }

alg_ch2.with_decay_card(decay_card_ch2).apply(sel_ch2)
alg_ch2.execute_on([jpsi_data, jpsi_incMC, exMC_ch2])

# ---------------------------------------------------------------------
# Channel 3 : Lambdabar p -> K_S0 2pi+ pi-  (Lambda -> p pi-, K_S0 -> pi+ pi-)
# ---------------------------------------------------------------------
alg_ch3 = Algorithm.new("JpsiLambdabarP_Ks2PiPi")
alg_ch3.set_header(["JpsiLambdabarP_Ks2PiPiAlg/JpsiLambdabarP_Ks2PiPi.h"])
       .set_constant({"ECMS" => [:double, 3.097]})
       .set_alias({"std::vector<double>" => "Vdouble"})
alg_ch3.note(:recoil_mass_window,
             "Lambda recoil mass |M_recoil(Lambda) - m_Lambdabar| < 0.020 GeV/c^2 "
             "(m_Lambdabar = 1.1157 GeV/c^2) required to tag Lambdabar-p annihilation events")
       .note(:annihilation_vertex,
             "Lambda-bar p annihilation vertex required in 3.0 <= R_xy <= 3.5 cm "
             "(annihilation inside the beam-pipe cooling-oil layer)")
       .note(:oil_momentum,
             "|p_oil| < 0.04 GeV/c required to reject Lambdabar annihilations on bound Au/Be/C nucleons")

sel_ch3 = Selection.new
  .select_track {
     cos_theta 0.93
     nChrp  ">=3"
     nChrn  ">=2"
     nTot   ">=5"
  }
  .pid(method: :probability) {
     prob_cut 0.001
     identify :proton, against: [:kaon, :pion]
     identify :pion,   against: [:kaon, :proton]
     nprp ">=1"
  }
  .invariant_mass_of(:prp, :pim).within(1.1127, 1.1187)
  .secondary_vertex_fit([:prp, :pim]) {
     build_virtual_particle(:Lambda).by_minimizing_mass_difference
     remove_used_particle_from_candidate_list
  }
  .invariant_mass_of(:pip, :pim).within(0.4826, 0.5126)
  .secondary_vertex_fit([:pip, :pim]) {
     build_virtual_particle(:K_S0).by_minimizing_mass_difference
     remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :K_S0, :pip, :pip, :pim]) {
     nominal
     constrain_four_momentum
     chi2_cut 200
  }

alg_ch3.with_decay_card(decay_card_ch3).apply(sel_ch3)
alg_ch3.execute_on([jpsi_data, jpsi_incMC, exMC_ch3])

# ---------------------------------------------------------------------
# Channel 4 : Lambdabar p -> K+ 2pi+ 2pi-   (Lambda -> p pi-)
# ---------------------------------------------------------------------
alg_ch4 = Algorithm.new("JpsiLambdabarP_K2Pi2Pi")
alg_ch4.set_header(["JpsiLambdabarP_K2Pi2PiAlg/JpsiLambdabarP_K2Pi2Pi.h"])
       .set_constant({"ECMS" => [:double, 3.097]})
       .set_alias({"std::vector<double>" => "Vdouble"})
alg_ch4.note(:recoil_mass_window,
             "Lambda recoil mass |M_recoil(Lambda) - m_Lambdabar| < 0.020 GeV/c^2 "
             "(m_Lambdabar = 1.1157 GeV/c^2) required to tag Lambdabar-p annihilation events")
       .note(:annihilation_vertex,
             "Lambda-bar p annihilation vertex required in 3.0 <= R_xy <= 3.5 cm "
             "(annihilation inside the beam-pipe cooling-oil layer)")
       .note(:oil_momentum,
             "|p_oil| < 0.04 GeV/c required to reject Lambdabar annihilations on bound Au/Be/C nucleons")

sel_ch4 = Selection.new
  .select_track {
     cos_theta 0.93
     nChrp  ">=3"
     nChrn  ">=3"
     nTot   ">=6"
  }
  .pid(method: :probability) {
     prob_cut 0.001
     identify :proton, against: [:kaon, :pion]
     identify :kaon,   against: [:pion, :proton]
     identify :pion,   against: [:kaon, :proton]
     nprp ">=1"
     nkp  ">=1"
  }
  .invariant_mass_of(:prp, :pim).within(1.1127, 1.1187)
  .secondary_vertex_fit([:prp, :pim]) {
     build_virtual_particle(:Lambda).by_minimizing_mass_difference
     remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :kp, :pip, :pip, :pim, :pim]) {
     nominal
     constrain_four_momentum
     chi2_cut 200
  }

alg_ch4.with_decay_card(decay_card_ch4).apply(sel_ch4)
alg_ch4.execute_on([jpsi_data, jpsi_incMC, exMC_ch4])

# ---------------------------------------------------------------------
# Channel 5 : Lambdabar p -> K+ 3pi+ 3pi-   (Lambda -> p pi-)
# ---------------------------------------------------------------------
alg_ch5 = Algorithm.new("JpsiLambdabarP_K3Pi3Pi")
alg_ch5.set_header(["JpsiLambdabarP_K3Pi3PiAlg/JpsiLambdabarP_K3Pi3Pi.h"])
       .set_constant({"ECMS" => [:double, 3.097]})
       .set_alias({"std::vector<double>" => "Vdouble"})
alg_ch5.note(:recoil_mass_window,
             "Lambda recoil mass |M_recoil(Lambda) - m_Lambdabar| < 0.020 GeV/c^2 "
             "(m_Lambdabar = 1.1157 GeV/c^2) required to tag Lambdabar-p annihilation events")
       .note(:annihilation_vertex,
             "Lambda-bar p annihilation vertex required in 3.0 <= R_xy <= 3.5 cm "
             "(annihilation inside the beam-pipe cooling-oil layer)")
       .note(:oil_momentum,
             "|p_oil| < 0.04 GeV/c required to reject Lambdabar annihilations on bound Au/Be/C nucleons")

sel_ch5 = Selection.new
  .select_track {
     cos_theta 0.93
     nChrp  ">=4"
     nChrn  ">=4"
     nTot   ">=8"
  }
  .pid(method: :probability) {
     prob_cut 0.001
     identify :proton, against: [:kaon, :pion]
     identify :kaon,   against: [:pion, :proton]
     identify :pion,   against: [:kaon, :proton]
     nprp ">=1"
     nkp  ">=1"
  }
  .invariant_mass_of(:prp, :pim).within(1.1127, 1.1187)
  .secondary_vertex_fit([:prp, :pim]) {
     build_virtual_particle(:Lambda).by_minimizing_mass_difference
     remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :kp, :pip, :pip, :pip, :pim, :pim, :pim]) {
     nominal
     constrain_four_momentum
     chi2_cut 200
  }

alg_ch5.with_decay_card(decay_card_ch5).apply(sel_ch5)
alg_ch5.execute_on([jpsi_data, jpsi_incMC, exMC_ch5])

# ---------------------------------------------------------------------
# Channel 6 : Lambdabar p -> 2K+ K-          (Lambda -> p pi-)
# ---------------------------------------------------------------------
alg_ch6 = Algorithm.new("JpsiLambdabarP_2KK")
alg_ch6.set_header(["JpsiLambdabarP_2KKAlg/JpsiLambdabarP_2KK.h"])
       .set_constant({"ECMS" => [:double, 3.097]})
       .set_alias({"std::vector<double>" => "Vdouble"})
alg_ch6.note(:recoil_mass_window,
             "Lambda recoil mass |M_recoil(Lambda) - m_Lambdabar| < 0.020 GeV/c^2 "
             "(m_Lambdabar = 1.1157 GeV/c^2) required to tag Lambdabar-p annihilation events")
       .note(:annihilation_vertex,
             "Lambda-bar p annihilation vertex required in 3.0 <= R_xy <= 3.5 cm "
             "(annihilation inside the beam-pipe cooling-oil layer)")
       .note(:oil_momentum,
             "|p_oil| < 0.04 GeV/c required to reject Lambdabar annihilations on bound Au/Be/C nucleons")

sel_ch6 = Selection.new
  .select_track {
     cos_theta 0.93
     nChrp  ">=3"
     nChrn  ">=2"
     nTot   ">=5"
  }
  .pid(method: :probability) {
     prob_cut 0.001
     identify :proton, against: [:kaon, :pion]
     identify :kaon,   against: [:pion, :proton]   # two K+ and one K-
     identify :pion,   against: [:kaon, :proton]
     nprp ">=1"
     nkp  ">=2"
     nkm  ">=1"
  }
  .invariant_mass_of(:prp, :pim).within(1.1127, 1.1187)
  .secondary_vertex_fit([:prp, :pim]) {
     build_virtual_particle(:Lambda).by_minimizing_mass_difference
     remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :kp, :kp, :km]) {
     nominal
     constrain_four_momentum
     chi2_cut 200
  }

alg_ch6.with_decay_card(decay_card_ch6).apply(sel_ch6)
alg_ch6.execute_on([jpsi_data, jpsi_incMC, exMC_ch6])