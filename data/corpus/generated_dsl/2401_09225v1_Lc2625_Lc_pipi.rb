### Dataset description ###
# Real data (+ matching inclusive MC) at sqrt(s) = 4.918 GeV (BOSS 707, sample 4914)
# and sqrt(s) = 4.950 GeV (BOSS 707, sample 4946)
data_4918  = DatasetManager.real_data.find("707_4914")        # sqrt(s) = 4.918 GeV
data_4950  = DatasetManager.real_data.find("707_4946")        # sqrt(s) = 4.950 GeV
incMC_4918 = DatasetManager.inclusive_mc.find("707_4914")     # matching inclusive MC
incMC_4950 = DatasetManager.inclusive_mc.find("707_4946")     # matching inclusive MC

# Decay card 1 : e+e- -> anti-Lambda_c- Lambda_c(2625)+ ; Lambda_c+ -> p+ K- pi+
decay_card_pKpi = <<~DECAYCARD
    Decay psi(4260)
    1.0 anti-Lambda_c- Lambda_c(2625)+ PHSP;
    Enddecay

    Decay Lambda_c(2625)+
    1.0 Lambda_c+ pi+ pi- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0 p+ K- pi+ PHSP;
    Enddecay

    End
DECAYCARD

# Decay card 2 : Lambda_c+ -> p+ K_S0 (K_S0 -> pi+ pi-)
decay_card_pKS0 = <<~DECAYCARD
    Decay psi(4260)
    1.0 anti-Lambda_c- Lambda_c(2625)+ PHSP;
    Enddecay

    Decay Lambda_c(2625)+
    1.0 Lambda_c+ pi+ pi- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0 p+ K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card 3 : Lambda_c+ -> Lambda0 pi+ (Lambda0 -> p+ pi-)
decay_card_Lpi = <<~DECAYCARD
    Decay psi(4260)
    1.0 anti-Lambda_c- Lambda_c(2625)+ PHSP;
    Enddecay

    Decay Lambda_c(2625)+
    1.0 Lambda_c+ pi+ pi- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0 Lambda0 pi+ PHSP;
    Enddecay

    Decay Lambda0
    1.0 p+ pi- HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC : 100k events at EACH energy point, for EACH Lambda_c+ tag mode
exMC_pKpi = DatasetManager.create_exclusive_mc_for([data_4918, data_4950]) do |config|
  config.sample_name   = "sig_Lc2625_pKpi"    # -> sig_Lc2625_pKpi_707_4914 / _707_4946
  config.events        = 100000
  config.decay_card    = decay_card_pKpi
  config.cross_section = :default
end

exMC_pKS0 = DatasetManager.create_exclusive_mc_for([data_4918, data_4950]) do |config|
  config.sample_name   = "sig_Lc2625_pKS0"
  config.events        = 100000
  config.decay_card    = decay_card_pKS0
  config.cross_section = :default
end

exMC_Lpi = DatasetManager.create_exclusive_mc_for([data_4918, data_4950]) do |config|
  config.sample_name   = "sig_Lc2625_Lpi"
  config.events        = 100000
  config.decay_card    = decay_card_Lpi
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ---------------------------------------------------------------------------
# Mode I : Lambda_c+ -> p+ K- pi+
# ---------------------------------------------------------------------------
alg_pKpi = Algorithm.new("Lc2625pKpi")
alg_pKpi.set_header(["Lc2625pKpiAlg/Lc2625pKpi.h"])
        .set_constant({"ECMS" => [:double, 4.918]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:ecms_multi_point, "real data combine two energy points (4.918 and 4.950 GeV); the single ECMS constant is set to 4.918 GeV and must be re-set to 4.950 GeV when the second point is processed")

sel_pKpi = Selection.new
  .select_track {                       # charged track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm (transverse plane)
    nNet      "==0"                     # net charge zero
    nChrp     ">=4"                     # >= 4 positive tracks
    nChrn     ">=2"                     # >= 2 negative tracks
  }
  .select_photon {                      # photon selection (no minimum multiplicity)
    tdc_emc_start     0                 # EMC TDC in [0, 14]
    tdc_emc_end       14
    angle_to_track    10.0              # > 10 deg from any charged track
    energyThreshold_b 0.025             # E > 25 MeV in the barrel
    energyThreshold_e 0.050             # E > 50 MeV in the endcap
  }
  .pid(method: :probability) {          # PID by probability method
    prob_cut 0.001                      # probability > 0.001
    identify :proton, against: [:kaon, :pion]     # explicit p / K / pi separation
    identify :kaon,   against: [:pion, :proton]
    identify :pion,   against: [:kaon, :proton]
  }
  # partial reconstruction: the anti-Lambda_c- (recID 1) and its daughters are left
  # unreconstructed; pick the Lambda_c(2625)+ combination closest to its nominal mass
  # and require the recoil mass against the missing anti-Lambda_c-
  .partial_miss([1]) {
    best_combination_by_mass :Lambda_c2625, 2.62811   # M(Lambda_c+ pi+ pi-) closest to 2.62811 GeV
    require_recoil_mass 2.000, 2.800                  # recoil mass in [2.000, 2.800] GeV
  }

alg_pKpi.with_decay_card(decay_card_pKpi).apply(sel_pKpi)

# ---------------------------------------------------------------------------
# Mode II : Lambda_c+ -> p+ K_S0 (K_S0 -> pi+ pi-)
# ---------------------------------------------------------------------------
alg_pKS0 = Algorithm.new("Lc2625pKS0")
alg_pKS0.set_header(["Lc2625pKS0Alg/Lc2625pKS0.h"])
        .set_constant({"ECMS" => [:double, 4.918]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:ecms_multi_point, "real data combine two energy points (4.918 and 4.950 GeV); the single ECMS constant is set to 4.918 GeV and must be re-set to 4.950 GeV when the second point is processed")

sel_pKS0 = Selection.new
  .select_track {                       # charged track selection
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nNet      "==0"
    nChrp     ">=4"                     # >= 4 positive tracks
    nChrn     ">=4"                     # >= 4 negative tracks (secondary-vertex mode)
  }
  .select_photon {                      # photon selection (no minimum multiplicity)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {          # PID by probability method
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]     # explicit p / K / pi separation
    identify :kaon,   against: [:pion, :proton]
    identify :pion,   against: [:kaon, :proton]
  }
  .secondary_vertex_fit([:pip, :pim]) { # build K_S0 from the pi+ pi- pair
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # partial reconstruction: anti-Lambda_c- (recID 1) left unreconstructed
  .partial_miss([1]) {
    best_combination_by_mass :Lambda_c2625, 2.62811
    require_recoil_mass 2.000, 2.800
  }

alg_pKS0.with_decay_card(decay_card_pKS0).apply(sel_pKS0)

# ---------------------------------------------------------------------------
# Mode III : Lambda_c+ -> Lambda0 pi+ (Lambda0 -> p+ pi-)
# ---------------------------------------------------------------------------
alg_Lpi = Algorithm.new("Lc2625Lpi")
alg_Lpi.set_header(["Lc2625LpiAlg/Lc2625Lpi.h"])
       .set_constant({"ECMS" => [:double, 4.918]})
       .set_alias({"std::vector<double>" => "Vdouble"})
       .note(:ecms_multi_point, "real data combine two energy points (4.918 and 4.950 GeV); the single ECMS constant is set to 4.918 GeV and must be re-set to 4.950 GeV when the second point is processed")

sel_Lpi = Selection.new
  .select_track {                       # charged track selection
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nNet      "==0"
    nChrp     ">=4"                     # >= 4 positive tracks
    nChrn     ">=4"                     # >= 4 negative tracks (secondary-vertex mode)
  }
  .select_photon {                      # photon selection (no minimum multiplicity)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {          # PID by probability method
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]     # explicit p / K / pi separation
    identify :kaon,   against: [:pion, :proton]
    identify :pion,   against: [:kaon, :proton]
  }
  .secondary_vertex_fit([:prp, :pim]) { # build Lambda0 from the p+ pi- pair
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # partial reconstruction: anti-Lambda_c- (recID 1) left unreconstructed
  .partial_miss([1]) {
    best_combination_by_mass :Lambda_c2625, 2.62811
    require_recoil_mass 2.000, 2.800
  }

alg_Lpi.with_decay_card(decay_card_Lpi).apply(sel_Lpi)

### Execute on real data, inclusive MC and signal MC ###
root_files_pKpi = alg_pKpi.execute_on([data_4918, data_4950, incMC_4918, incMC_4950, *exMC_pKpi])
root_files_pKS0 = alg_pKS0.execute_on([data_4918, data_4950, incMC_4918, incMC_4950, *exMC_pKS0])
root_files_Lpi  = alg_Lpi.execute_on([data_4918, data_4950, incMC_4918, incMC_4950, *exMC_Lpi])

# NOTE: signal extraction (Delta M = M(Lambda_c+ pi+ pi-) - M(Lambda_c+)) is performed
# downstream in the ROOT analysis, not in BOSS.