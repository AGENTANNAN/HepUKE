# =====================================================================
# BESIII (proceedings): Study of charmonium decays at BESIII
#   Search for baryonic decays of psi(3770) and psi(4040):
#     Lambda Lambdabar pi+pi-, Lambda Lambdabar pi0, Lambda Lambdabar eta,
#     Sigma+ Sigmabar-, Sigma0 Sigmabar0, Xi- Xibar+, Xi0 Xibar0
#   Data: 2.9 fb^-1 at sqrt(s) = 3.773 GeV, 482 pb^-1 at 4.009 GeV and
#         67 pb^-1 of continuum data at 3.542, 3.554, 3.561, 3.600, 3.650 GeV
#
# Each final state needs a different track/photon multiplicity and a
# different reconstruction chain, so one Algorithm object is created per
# final state (Rule T1); each algorithm is executed on the psi(3770), the
# psi(4040) and the continuum datasets. The resonance energy is injected
# per dataset at run time (multi-energy / cross-section analysis), so ECMS
# is not fixed in set_constant.
# =====================================================================

### Dataset description ###
data_3773     = DatasetManager.real_data.find("712_3773")   # psi(3770), 2.9 fb^-1
data_4009     = DatasetManager.real_data.find("703_4009")   # psi(4040), 482 pb^-1
data_cont     = DatasetManager.real_data.find("709_3650")   # continuum 3.542-3.650 GeV, 67 pb^-1

datasets_baryonic = [data_3773, data_4009, data_cont]

# ---------------------------------------------------------------------
# Decay cards (EvtGen), top mother psi(3770) (the psi(4040) samples use
# the identical cards with psi(4040) as the top mother)
# ---------------------------------------------------------------------
# Mode A: psi(3770) -> Lambda Lambdabar pi+ pi-
decay_card_A = <<~DECAYCARD
    Decay psi(3770)
    1.0000 Lambda0 anti-Lambda0 pi+ pi- PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Mode B: psi(3770) -> Lambda Lambdabar pi0
decay_card_B = <<~DECAYCARD
    Decay psi(3770)
    1.0000 Lambda0 anti-Lambda0 pi0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode C: psi(3770) -> Lambda Lambdabar eta
decay_card_C = <<~DECAYCARD
    Decay psi(3770)
    1.0000 Lambda0 anti-Lambda0 eta PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode D: psi(3770) -> Sigma+ Sigmabar-, Sigma+ -> p pi0, Sigmabar- -> pbar pi0
decay_card_D = <<~DECAYCARD
    Decay psi(3770)
    1.0000 Sigma+ anti-Sigma- PHSP;
    Enddecay

    Decay Sigma+
    1.0000 p+ pi0 PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000 anti-p- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode E: psi(3770) -> Sigma0 Sigmabar0, Sigma0 -> Lambda gamma
decay_card_E = <<~DECAYCARD
    Decay psi(3770)
    1.0000 Sigma0 anti-Sigma0 PHSP;
    Enddecay

    Decay Sigma0
    1.0000 Lambda0 gamma PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 anti-Lambda0 gamma PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Mode F: psi(3770) -> Xi- Xibar+, Xi- -> Lambda pi-, Xibar+ -> Lambdabar pi+
decay_card_F = <<~DECAYCARD
    Decay psi(3770)
    1.0000 anti-Xi+ Xi- PHSP;
    Enddecay

    Decay Xi-
    1.0000 Lambda0 pi- PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 anti-Lambda0 pi+ PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Mode G: psi(3770) -> Xi0 Xibar0, Xi0 -> Lambda pi0, Xibar0 -> Lambdabar pi0
decay_card_G = <<~DECAYCARD
    Decay psi(3770)
    1.0000 anti-Xi0 Xi0 PHSP;
    Enddecay

    Decay Xi0
    1.0000 Lambda0 pi0 PHSP;
    Enddecay

    Decay anti-Xi0
    1.0000 anti-Lambda0 pi0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: one signal sample per mode, generated at the psi(3770) energy
exMC_A = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_LambdaLambdabar_pipi"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_A
  config.cross_section   = :default
end

exMC_B = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_LambdaLambdabar_pi0"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_B
  config.cross_section   = :default
end

exMC_C = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_LambdaLambdabar_eta"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_C
  config.cross_section   = :default
end

exMC_D = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_SigmaPlusSigmabarMinus"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_D
  config.cross_section   = :default
end

exMC_E = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_Sigma0Sigmabar0"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_E
  config.cross_section   = :default
end

exMC_F = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_XiMinusXibarPlus"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_F
  config.cross_section   = :default
end

exMC_G = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_Xi0Xibar0"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_G
  config.cross_section   = :default
end

# Exclusive MC produced at the psi(4040) energy (same cards, psi(4040) as
# the top mother); created with the batch helper over the 4.009 GeV point
exMCs_4040_A = DatasetManager.create_exclusive_mc_for([data_4009]) do |config|
  config.sample_name   = "exmc_4040_LambdaLambdabar_pipi"
  config.events        = 100_000
  config.decay_card    = decay_card_A.sub("Decay psi(3770)", "Decay psi(4040)")
  config.cross_section = :default
end

exMCs_4040_D = DatasetManager.create_exclusive_mc_for([data_4009]) do |config|
  config.sample_name   = "exmc_4040_SigmaPlusSigmabarMinus"
  config.events        = 100_000
  config.decay_card    = decay_card_D.sub("Decay psi(3770)", "Decay psi(4040)")
  config.cross_section = :default
end

exMCs_4040_F = DatasetManager.create_exclusive_mc_for([data_4009]) do |config|
  config.sample_name   = "exmc_4040_XiMinusXibarPlus"
  config.events        = 100_000
  config.decay_card    = decay_card_F.sub("Decay psi(3770)", "Decay psi(4040)")
  config.cross_section = :default
end

# ---------------------------------------------------------------------
### Event selection (BOSS) ###
# ---------------------------------------------------------------------

# =====================================================================
# Mode A -- psi(3770)/psi(4040) -> Lambda Lambdabar pi+ pi-
# =====================================================================
algA_name = "LambdaLambdabarPiPi"
algA = Algorithm.new(algA_name)
algA.set_header(["#{algA_name}Alg/#{algA_name}.h"])
    .set_alias({"std::vector<double>" => "Vdouble"})

selA = Selection.new
selA.select_track do
       cos_theta 0.93   # |cos(theta)| < 0.93
       Vz        10.0   # |Vz| < 10 cm
       Vr        1.0    # Vr < 1 cm
       nChrp     "==2"  # p and pi+
       nChrn     "==2"  # pi- and p-bar
       nNet      "==0"
     end
     .pid(method: :probability) do
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       nprp     ">=1"
       nprm     ">=1"
     end
     .secondary_vertex_fit([:prp, :pim]) do
       build_virtual_particle(:Lambda).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     end
     .secondary_vertex_fit([:prm, :pip]) do
       build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     end
     .assign({chrgp: :pip, chrgn: :pim})
     .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim]) do
       nominal
       constrain_four_momentum
       chi2_cut 200   # loose cut in BOSS; optimal cut applied in ROOT
     end

algA.with_decay_card(decay_card_A).apply(selA)

# =====================================================================
# Mode B -- psi(3770)/psi(4040) -> Lambda Lambdabar pi0, pi0 -> gamma gamma
# =====================================================================
algB_name = "LambdaLambdabarPi0"
algB = Algorithm.new(algB_name)
algB.set_header(["#{algB_name}Alg/#{algB_name}.h"])
    .set_alias({"std::vector<double>" => "Vdouble"})

selB = Selection.new
selB.select_track do
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==2"
       nChrn     "==2"
       nNet      "==0"
     end
     .select_photon do
       tdc_emc_start     0
       tdc_emc_end       14
       angle_to_track    10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam              ">=2"   # two photons from pi0 -> gamma gamma
     end
     .pid(method: :probability) do
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       nprp     ">=1"
       nprm     ">=1"
     end
     .secondary_vertex_fit([:prp, :pim]) do
       build_virtual_particle(:Lambda).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     end
     .secondary_vertex_fit([:prm, :pip]) do
       build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     end
     .assign({chrgp: :pip, chrgn: :pim})
     # pi0 reconstructed from the photon pair via a 1C Kalman mass constraint
     .kalman_kinematic_fit([:gamma, :gamma]) do
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
       chi2_cut 25
       npi0 ">=1"
     end
     .kinematic_fit([:Lambda, :Lambda_bar, :pi0]) do
       nominal
       constrain_four_momentum
       chi2_cut 200
     end

algB.with_decay_card(decay_card_B).apply(selB)

# =====================================================================
# Mode C -- psi(3770)/psi(4040) -> Lambda Lambdabar eta, eta -> gamma gamma
# =====================================================================
algC_name = "LambdaLambdabarEta"
algC = Algorithm.new(algC_name)
algC.set_header(["#{algC_name}Alg/#{algC_name}.h"])
    .set_alias({"std::vector<double>" => "Vdouble"})

selC = Selection.new
selC.select_track do
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==2"
       nChrn     "==2"
       nNet      "==0"
     end
     .select_photon do
       tdc_emc_start     0
       tdc_emc_end       14
       angle_to_track    10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam              ">=2"   # two photons from eta -> gamma gamma
     end
     .pid(method: :probability) do
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       nprp     ">=1"
       nprm     ">=1"
     end
     .secondary_vertex_fit([:prp, :pim]) do
       build_virtual_particle(:Lambda).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     end
     .secondary_vertex_fit([:prm, :pip]) do
       build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     end
     .assign({chrgp: :pip, chrgn: :pim})
     # eta reconstructed from the photon pair via a 1C Kalman mass constraint
     .kalman_kinematic_fit([:gamma, :gamma]) do
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
       chi2_cut 25
       neta ">=1"
     end
     .kinematic_fit([:Lambda, :Lambda_bar, :eta]) do
       nominal
       constrain_four_momentum
       chi2_cut 200
     end

algC.with_decay_card(decay_card_C).apply(selC)

# =====================================================================
# Mode D -- psi(3770)/psi(4040) -> Sigma+ Sigmabar-,
#           Sigma+ -> p pi0, Sigmabar- -> pbar pi0, pi0 -> gamma gamma
# =====================================================================
algD_name = "SigmaPlusSigmabarMinus"
algD = Algorithm.new(algD_name)
algD.set_header(["#{algD_name}Alg/#{algD_name}.h"])
    .set_alias({"std::vector<double>" => "Vdouble"})

selD = Selection.new
selD.select_track do
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==1"  # p from Sigma+
       nChrn     "==1"  # p-bar from Sigmabar-
       nNet      "==0"
     end
     .select_photon do
       tdc_emc_start     0
       tdc_emc_end       14
       angle_to_track    10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam              ">=4"   # two pi0 -> four photons
     end
     .pid(method: :probability) do
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       nprp     "==1"
       nprm     "==1"
     end
     # Two pi0 candidates from the four photons (1C Kalman mass constraint)
     .kalman_kinematic_fit([:gamma, :gamma]) do
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
       chi2_cut 25
       npi0 ">=2"
     end
     .kinematic_fit([:prp, :prm, :pi0, :pi0]) do
       nominal
       constrain_four_momentum
       chi2_cut 200
     end

algD.with_decay_card(decay_card_D).apply(selD)

# =====================================================================
# Mode E -- psi(3770)/psi(4040) -> Sigma0 Sigmabar0,
#           Sigma0 -> Lambda gamma, Lambda -> p pi-
# =====================================================================
algE_name = "Sigma0Sigmabar0"
algE = Algorithm.new(algE_name)
algE.set_header(["#{algE_name}Alg/#{algE_name}.h"])
    .set_alias({"std::vector<double>" => "Vdouble"})

selE = Selection.new
selE.select_track do
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==2"  # p and pi+
       nChrn     "==2"  # pi- and p-bar
       nNet      "==0"
     end
     .select_photon do
       tdc_emc_start     0
       tdc_emc_end       14
       angle_to_track    10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam              ">=2"   # one photon from each Sigma0 -> Lambda gamma
     end
     .pid(method: :probability) do
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       nprp     ">=1"
       nprm     ">=1"
     end
     # Lambda / Lambdabar secondary vertices
     .secondary_vertex_fit([:prp, :pim]) do
       build_virtual_particle(:Lambda).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     end
     .secondary_vertex_fit([:prm, :pip]) do
       build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     end
     .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) do
       nominal
       constrain_four_momentum
       invariant_mass_of(:Lambda).within(1.111, 1.121)   # Lambda mass window
       invariant_mass_of(:Lambda_bar).within(1.111, 1.121)
       chi2_cut 200
     end

algE.with_decay_card(decay_card_E).apply(selE)

# =====================================================================
# Mode F -- psi(3770)/psi(4040) -> Xi- Xibar+,
#           Xi- -> Lambda pi-, Xibar+ -> Lambdabar pi+
# =====================================================================
algF_name = "XiMinusXibarPlus"
algF = Algorithm.new(algF_name)
algF.set_header(["#{algF_name}Alg/#{algF_name}.h"])
    .set_alias({"std::vector<double>" => "Vdouble"})

selF = Selection.new
selF.select_track do
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==3"  # p, pi+ (from Lambdabar), pi+ (from Xibar+)
       nChrn     "==3"  # pi-, pi- (from Xi-), p-bar
       nNet      "==0"
     end
     .pid(method: :probability) do
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       nprp     ">=1"
       nprm     ">=1"
     end
     # Lambda -> p pi- vertex (the two pions from the Xi- remain unassigned)
     .secondary_vertex_fit([:prp, :pim]) do
       build_virtual_particle(:Lambda).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     end
     .secondary_vertex_fit([:prm, :pip]) do
       build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     end
     .assign({chrgp: :pip, chrgn: :pim})
     # The two extra pions and the Lambdas are combined in the 4C fit; the
     # best Xi- / Xibar+ assignment is the one with the smallest chi2
     .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim, :pim, :pip]) do
       nominal
       constrain_four_momentum
       invariant_mass_of(:Lambda).within(1.111, 1.121)
       invariant_mass_of(:Lambda_bar).within(1.111, 1.121)
       chi2_cut 200
     end

algF.with_decay_card(decay_card_F).apply(selF)

# =====================================================================
# Mode G -- psi(3770)/psi(4040) -> Xi0 Xibar0,
#           Xi0 -> Lambda pi0, Xibar0 -> Lambdabar pi0, pi0 -> gamma gamma
# =====================================================================
algG_name = "Xi0Xibar0"
algG = Algorithm.new(algG_name)
algG.set_header(["#{algG_name}Alg/#{algG_name}.h"])
    .set_alias({"std::vector<double>" => "Vdouble"})

selG = Selection.new
selG.select_track do
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==2"  # p and pi+
       nChrn     "==2"  # pi- and p-bar
       nNet      "==0"
     end
     .select_photon do
       tdc_emc_start     0
       tdc_emc_end       14
       angle_to_track    10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam              ">=4"   # two pi0 -> four photons
     end
     .pid(method: :probability) do
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       nprp     ">=1"
       nprm     ">=1"
     end
     .secondary_vertex_fit([:prp, :pim]) do
       build_virtual_particle(:Lambda).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     end
     .secondary_vertex_fit([:prm, :pip]) do
       build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     end
     # Two pi0 candidates from the four photons
     .kalman_kinematic_fit([:gamma, :gamma]) do
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
       chi2_cut 25
       npi0 ">=2"
     end
     .kinematic_fit([:Lambda, :Lambda_bar, :pi0, :pi0]) do
       nominal
       constrain_four_momentum
       invariant_mass_of(:Lambda).within(1.111, 1.121)
       invariant_mass_of(:Lambda_bar).within(1.111, 1.121)
       chi2_cut 200
     end

algG.with_decay_card(decay_card_G).apply(selG)

# ---------------------------------------------------------------------
# Execution: every algorithm runs on the psi(3770), psi(4040) and
# continuum datasets; the continuum datasets are used for the continuum
# subtraction (the dominant source of uncertainty).
# ---------------------------------------------------------------------
exMCs = [exMC_A, exMC_B, exMC_C, exMC_D, exMC_E, exMC_F, exMC_G]

algA.execute_on(datasets_baryonic + [exMC_A] + exMCs_4040_A)
algB.execute_on(datasets_baryonic + [exMC_B])
algC.execute_on(datasets_baryonic + [exMC_C])
algD.execute_on(datasets_baryonic + [exMC_D] + exMCs_4040_D)
algE.execute_on(datasets_baryonic + [exMC_E])
algF.execute_on(datasets_baryonic + [exMC_F] + exMCs_4040_F)
algG.execute_on(datasets_baryonic + [exMC_G])
