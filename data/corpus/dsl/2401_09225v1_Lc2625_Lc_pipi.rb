# Paper: arXiv:2401.09225v1
# First measurements of absolute BF of Lambda_c(2625)+ -> Lambda_c+ pi+ pi- and upper limit on Lambda_c(2595)+ -> Lambda_c+ pi+ pi-
# Data at sqrt(s) = 4.918 and 4.950 GeV
# Ordinary analysis with partial reconstruction (missing anti-Lambda_c-)

### Dataset preparation ###
data_4918 = DatasetManager.real_data.find("707_4914")  # 4.918 GeV
data_4950 = DatasetManager.real_data.find("707_4946")  # 4.950 GeV
incMC_4918 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4950 = DatasetManager.inclusive_mc.find("707_4946")

# ============================================================
# Mode 1: Lambda_c+ -> p+ K- pi+  (direct charged tracks)
# ============================================================

decay_card_mode1 = <<~DECAYCARD
    Decay psi(4260)
    1.000 anti-Lambda_c- Lambda_c_star+ PHSP;
    Enddecay
    Decay Lambda_c_star+
    1.000 Lambda_c+ pi+ pi- PHSP;
    Enddecay
    Decay Lambda_c+
    1.000 p+ K- pi+ PHSP;
    Enddecay
    Decay anti-Lambda_c-
    1.000 anti-p- K+ pi- PHSP;
    Enddecay
    End
DECAYCARD

exMC_mode1_4918 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_Lc2625_Lc_pKpi_pipi_4918"
  config.related_dataset = data_4918
  config.events = 100000
  config.decay_card = decay_card_mode1
  config.cross_section = :default
end

exMC_mode1_4950 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_Lc2625_Lc_pKpi_pipi_4950"
  config.related_dataset = data_4950
  config.events = 100000
  config.decay_card = decay_card_mode1
  config.cross_section = :default
end

alg_mode1 = Algorithm.new("Lc2625LcPKPi")
alg_mode1.set_header(["Lc2625LcPKPiAlg/Lc2625LcPKPi.h"])
         .set_constant({"ECMS" => [:double, 4.918]})

event_selection_mode1 = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=4"   # p+, K+, pi+(from Lambda_c+), pi+(from Lc*), plus anti-p-, K+, pi- on missing side
    nChrn ">=2"   # K-, pi-(from Lambda_c+), pi-(from Lc*)
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
  end
  .partial_miss([1]) do    # miss anti-Lambda_c- and its daughters
    best_combination_by_mass :Lambda_c_star_plus, 2.62811
    require_recoil_mass 2.000, 2.800
  end

alg_mode1.note(:two_energy_points, "Analysis at 4.918 and 4.950 GeV; ECMS constant is approximate placeholder (4.918 GeV shown); multi-energy run handles per-run beam energy via MeasuredEcmsSvc")
alg_mode1.note(:tag_approach, "Lambda_c+ tagged via pK-pi+ hadronic mode; two more tag modes (pK_S0, Lambda pi+) in separate algorithms")
alg_mode1.note(:partial_rec_method, "Partial reconstruction with missing anti-Lambda_c-; recoil mass window applied; N_tag determined from recoil mass fit (in ROOT)")
alg_mode1.note(:signal_yield, "Signal yield N_sig from full reconstruction of Lambda_c(2625)+ -> Lambda_c+ pi+ pi-; Delta_M = M(Lc+ pi+ pi-) - M(Lc+) used for selection (ROOT analysis)")
alg_mode1.note(:custom_generator, "Signal MC uses KKMC generator with cross-section line shapes from BESIII measurements of e+e- -> anti-Lambda_c- Lambda_c(2595/2625)+")

alg_mode1.with_decay_card(decay_card_mode1).apply(event_selection_mode1)
alg_mode1.execute_on([data_4918, data_4950, incMC_4918, incMC_4950, exMC_mode1_4918, exMC_mode1_4950])

# ============================================================
# Mode 2: Lambda_c+ -> p+ K_S0  (K_S0 -> pi+ pi- vertex fit)
# ============================================================

decay_card_mode2 = <<~DECAYCARD
    Decay psi(4260)
    1.000 anti-Lambda_c- Lambda_c_star+ PHSP;
    Enddecay
    Decay Lambda_c_star+
    1.000 Lambda_c+ pi+ pi- PHSP;
    Enddecay
    Decay Lambda_c+
    1.000 p+ K_S0 PHSP;
    Enddecay
    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay
    Decay anti-Lambda_c-
    1.000 anti-p- K+ pi- PHSP;
    Enddecay
    End
DECAYCARD

exMC_mode2_4918 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_Lc2625_Lc_pKs0_pipi_4918"
  config.related_dataset = data_4918
  config.events = 100000
  config.decay_card = decay_card_mode2
  config.cross_section = :default
end

exMC_mode2_4950 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_Lc2625_Lc_pKs0_pipi_4950"
  config.related_dataset = data_4950
  config.events = 100000
  config.decay_card = decay_card_mode2
  config.cross_section = :default
end

alg_mode2 = Algorithm.new("Lc2625LcPKs")
alg_mode2.set_header(["Lc2625LcPKsAlg/Lc2625LcPKs.h"])
         .set_constant({"ECMS" => [:double, 4.918]})

event_selection_mode2 = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=4"
    nChrn ">=4"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .partial_miss([1]) do    # miss anti-Lambda_c- and its daughters
    best_combination_by_mass :Lambda_c_star_plus, 2.62811
    require_recoil_mass 2.000, 2.800
  end

alg_mode2.note(:two_energy_points, "Analysis at 4.918 and 4.950 GeV; ECMS constant is approximate placeholder")
alg_mode2.note(:tag_approach, "Lambda_c+ tagged via pK_S0 hadronic mode; K_S0 reconstructed via secondary vertex fit")
alg_mode2.note(:partial_rec_method, "Partial reconstruction with missing anti-Lambda_c-; recoil mass window applied")
alg_mode2.note(:signal_yield, "Signal yield from full reconstruction of Lambda_c(2625)+ -> Lambda_c+ pi+ pi-; Delta_M used for selection (ROOT)")
alg_mode2.note(:custom_generator, "Signal MC uses KKMC generator with cross-section line shapes from BESIII measurements")

alg_mode2.with_decay_card(decay_card_mode2).apply(event_selection_mode2)
alg_mode2.execute_on([data_4918, data_4950, incMC_4918, incMC_4950, exMC_mode2_4918, exMC_mode2_4950])

# ============================================================
# Mode 3: Lambda_c+ -> Lambda pi+  (Lambda -> p+ pi- vertex fit)
# ============================================================

decay_card_mode3 = <<~DECAYCARD
    Decay psi(4260)
    1.000 anti-Lambda_c- Lambda_c_star+ PHSP;
    Enddecay
    Decay Lambda_c_star+
    1.000 Lambda_c+ pi+ pi- PHSP;
    Enddecay
    Decay Lambda_c+
    1.000 Lambda0 pi+ PHSP;
    Enddecay
    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay
    Decay anti-Lambda_c-
    1.000 anti-p- K+ pi- PHSP;
    Enddecay
    End
DECAYCARD

exMC_mode3_4918 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_Lc2625_Lc_Lpi_pipi_4918"
  config.related_dataset = data_4918
  config.events = 100000
  config.decay_card = decay_card_mode3
  config.cross_section = :default
end

exMC_mode3_4950 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_Lc2625_Lc_Lpi_pipi_4950"
  config.related_dataset = data_4950
  config.events = 100000
  config.decay_card = decay_card_mode3
  config.cross_section = :default
end

alg_mode3 = Algorithm.new("Lc2625LcLambdaPi")
alg_mode3.set_header(["Lc2625LcLambdaPiAlg/Lc2625LcLambdaPi.h"])
         .set_constant({"ECMS" => [:double, 4.918]})

event_selection_mode3 = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=4"
    nChrn ">=4"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
  end
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .partial_miss([1]) do    # miss anti-Lambda_c- and its daughters
    best_combination_by_mass :Lambda_c_star_plus, 2.62811
    require_recoil_mass 2.000, 2.800
  end

alg_mode3.note(:two_energy_points, "Analysis at 4.918 and 4.950 GeV; ECMS constant is approximate placeholder")
alg_mode3.note(:tag_approach, "Lambda_c+ tagged via Lambda pi+ hadronic mode; Lambda reconstructed via secondary vertex fit")
alg_mode3.note(:partial_rec_method, "Partial reconstruction with missing anti-Lambda_c-; recoil mass window applied")
alg_mode3.note(:signal_yield, "Signal yield from full reconstruction of Lambda_c(2625)+ -> Lambda_c+ pi+ pi-; Delta_M used for selection (ROOT)")
alg_mode3.note(:custom_generator, "Signal MC uses KKMC generator with cross-section line shapes from BESIII measurements")

alg_mode3.with_decay_card(decay_card_mode3).apply(event_selection_mode3)
alg_mode3.execute_on([data_4918, data_4950, incMC_4918, incMC_4950, exMC_mode3_4918, exMC_mode3_4950])