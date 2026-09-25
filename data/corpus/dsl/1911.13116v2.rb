# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Analysis uses 2.93 fb-1 e+e- collision data at sqrt(s) = 3.773 GeV (psi(3770))
psipp_data = DatasetManager.real_data.find("712_3773")
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")

#============================================================================#
# Mode 1: D+ -> Lambda e+   (Delta(B-L) = 0)
#============================================================================#

decay_card_lambda_eplus = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- VSS;
    Enddecay

    Decay D+
    1.0000 Lambda e+ PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay Lambda
    1.0000 p+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_lambda_eplus = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "BNV_Lambda_eplus_exMC"
  config.related_dataset = psipp_data
  config.events = 100000
  config.decay_card = decay_card_lambda_eplus
  config.cross_section = :default
end

alg_lambda_eplus = Algorithm.new("BNV_Lambda_eplus")
alg_lambda_eplus.set_header(["BNV_Lambda_eplusAlg/BNV_Lambda_eplus.h"])
  .set_constant({"ECMS" => [:double, 3.773]})

selection_lambda_eplus = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 20.0          # Lambda daughter tracks require |Vz| < 20 cm
    Vr 1.0           # positron requires |Vr| < 1 cm
    nChrp ">=2"      # at least 2 positive tracks: proton (from Lambda) + e+
    nChrn ">=1"      # at least 1 negative track: pi- (from Lambda)
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]  # CL_p > 0.001, CL_p > CL_K, CL_p > CL_pi
    nprp "==1"       # exactly one proton; also enforces no extra proton/anti-proton in event
  end
  .remove([:prp <= :chrgp])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :electron, against: [:pion, :kaon, :proton]   # CL_e > 0.001
    nep ">=1"
  end
  .remove([:ep <= :chrgp])
  .assign({chrgp: :pip, chrgn: :pim})   # remaining positive tracks = pi+, negative = pi-
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:Lambda, :ep]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

alg_lambda_eplus.note(:positron_pid_combined_prob,
  "Paper requires CL_e/(CL_e+CL_pi+CL_K+CL_p) > 0.8 and 0.8 < E/p < 1.2 for positron. " \
  "These cuts are applied in ROOT analysis; only CL_e > 0.001 is expressed in the DSL PID block.")
alg_lambda_eplus.note(:lambda_mass_window,
  "Paper applies Lambda mass window (1.110, 1.121) GeV/c^2 after vertex fit. Applied in ROOT analysis.")
alg_lambda_eplus.note(:photon_conversion_veto,
  "Paper vetoes events where the positron forms a photon-conversion pair (Delta_xy, cos(theta_eg), R_xy cuts). Applied in ROOT analysis.")
alg_lambda_eplus.note(:deltaE_mbc,
  "Paper applies DeltaE (-0.023, 0.022) GeV and MBC selection with ARGUS background fit. Applied in ROOT analysis.")

alg_lambda_eplus.with_decay_card(decay_card_lambda_eplus).apply(selection_lambda_eplus)
root_files_lambda_eplus = alg_lambda_eplus.execute_on([psipp_data, psipp_incMC, exMC_lambda_eplus])

#============================================================================#
# Mode 2: D+ -> Sigma0 e+   (Delta(B-L) = 0, Sigma0 -> gamma Lambda)
#============================================================================#

decay_card_sigma0_eplus = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- VSS;
    Enddecay

    Decay D+
    1.0000 Sigma0 e+ PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda PHSP;
    Enddecay

    Decay Lambda
    1.0000 p+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_sigma0_eplus = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "BNV_Sigma0_eplus_exMC"
  config.related_dataset = psipp_data
  config.events = 100000
  config.decay_card = decay_card_sigma0_eplus
  config.cross_section = :default
end

alg_sigma0_eplus = Algorithm.new("BNV_Sigma0_eplus")
alg_sigma0_eplus.set_header(["BNV_Sigma0_eplusAlg/BNV_Sigma0_eplus.h"])
  .set_constant({"ECMS" => [:double, 3.773]})

selection_sigma0_eplus = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 20.0
    Vr 1.0
    nChrp ">=2"      # proton (from Lambda) + e+
    nChrn ">=1"      # pi- (from Lambda)
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=1"       # photon from Sigma0 -> gamma Lambda
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"       # exactly one proton; no extra proton in event
  end
  .remove([:prp <= :chrgp])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :electron, against: [:pion, :kaon, :proton]
    nep ">=1"
  end
  .remove([:ep <= :chrgp])
  .assign({chrgp: :pip, chrgn: :pim})
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:Lambda, :gamma, :ep]) do
    invariant_mass_of(:Lambda, :gamma).constrain_to_nominal_mass_of(:Sigma0)  # 1C Sigma0 mass constraint
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

alg_sigma0_eplus.note(:positron_pid_combined_prob,
  "Paper requires CL_e/(CL_e+CL_pi+CL_K+CL_p) > 0.8 and 0.8 < E/p < 1.2 for positron. " \
  "These cuts are applied in ROOT analysis; only CL_e > 0.001 is expressed in the DSL PID block.")
alg_sigma0_eplus.note(:lambda_sigma0_mass_window,
  "Paper applies Lambda mass window (1.110, 1.121) GeV/c^2 and Sigma0 mass window (1.173, 1.200) GeV/c^2. Applied in ROOT analysis.")
alg_sigma0_eplus.note(:photon_conversion_veto,
  "Paper vetoes events where the positron forms a photon-conversion pair (Delta_xy, cos(theta_eg), R_xy cuts). Applied in ROOT analysis.")
alg_sigma0_eplus.note(:deltaE_mbc,
  "Paper applies DeltaE (-0.028, 0.024) GeV and MBC selection with ARGUS background fit. Applied in ROOT analysis.")

alg_sigma0_eplus.with_decay_card(decay_card_sigma0_eplus).apply(selection_sigma0_eplus)
root_files_sigma0_eplus = alg_sigma0_eplus.execute_on([psipp_data, psipp_incMC, exMC_sigma0_eplus])

#============================================================================#
# Mode 3: D+ -> Lambdabar e+   (Delta(B-L) = 2)
#============================================================================#

decay_card_lambdabar_eplus = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- VSS;
    Enddecay

    Decay D+
    1.0000 anti-Lambda0 e+ PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ PHSP;
    Enddecay

    End
DECAYCARD

exMC_lambdabar_eplus = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "BNV_Lambdabar_eplus_exMC"
  config.related_dataset = psipp_data
  config.events = 100000
  config.decay_card = decay_card_lambdabar_eplus
  config.cross_section = :default
end

alg_lambdabar_eplus = Algorithm.new("BNV_Lambdabar_eplus")
alg_lambdabar_eplus.set_header(["BNV_Lambdabar_eplusAlg/BNV_Lambdabar_eplus.h"])
  .set_constant({"ECMS" => [:double, 3.773]})

selection_lambdabar_eplus = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 20.0
    Vr 1.0
    nChrp ">=2"      # pi+ (from Lambdabar) + e+
    nChrn ">=1"      # anti-p (from Lambdabar)
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]  # identifies both prp and prm
    nprm "==1"       # exactly one anti-proton; no extra anti-proton in event
  end
  .remove([:prm <= :chrgn])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :electron, against: [:pion, :kaon, :proton]
    nep ">=1"
  end
  .remove([:ep <= :chrgp])
  .assign({chrgp: :pip, chrgn: :pim})
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:Lambda_bar, :ep]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

alg_lambdabar_eplus.note(:positron_pid_combined_prob,
  "Paper requires CL_e/(CL_e+CL_pi+CL_K+CL_p) > 0.8 and 0.8 < E/p < 1.2 for positron. " \
  "These cuts are applied in ROOT analysis; only CL_e > 0.001 is expressed in the DSL PID block.")
alg_lambdabar_eplus.note(:lambda_mass_window,
  "Paper applies Lambdabar mass window (1.110, 1.121) GeV/c^2 after vertex fit. Applied in ROOT analysis.")
alg_lambdabar_eplus.note(:photon_conversion_veto,
  "Paper vetoes events where the positron forms a photon-conversion pair (Delta_xy, cos(theta_eg), R_xy cuts). Applied in ROOT analysis.")
alg_lambdabar_eplus.note(:deltaE_mbc,
  "Paper applies DeltaE (-0.023, 0.022) GeV and MBC selection with ARGUS background fit. Applied in ROOT analysis.")

alg_lambdabar_eplus.with_decay_card(decay_card_lambdabar_eplus).apply(selection_lambdabar_eplus)
root_files_lambdabar_eplus = alg_lambdabar_eplus.execute_on([psipp_data, psipp_incMC, exMC_lambdabar_eplus])

#============================================================================#
# Mode 4: D+ -> Sigmabar0 e+   (Delta(B-L) = 2, Sigmabar0 -> gamma Lambdabar)
#============================================================================#

decay_card_sigmabar0_eplus = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- VSS;
    Enddecay

    Decay D+
    1.0000 anti-Sigma0 e+ PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 gamma anti-Lambda0 PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ PHSP;
    Enddecay

    End
DECAYCARD

exMC_sigmabar0_eplus = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "BNV_Sigmabar0_eplus_exMC"
  config.related_dataset = psipp_data
  config.events = 100000
  config.decay_card = decay_card_sigmabar0_eplus
  config.cross_section = :default
end

alg_sigmabar0_eplus = Algorithm.new("BNV_Sigmabar0_eplus")
alg_sigmabar0_eplus.set_header(["BNV_Sigmabar0_eplusAlg/BNV_Sigmabar0_eplus.h"])
  .set_constant({"ECMS" => [:double, 3.773]})

selection_sigmabar0_eplus = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 20.0
    Vr 1.0
    nChrp ">=2"      # pi+ (from Lambdabar) + e+
    nChrn ">=1"      # anti-p (from Lambdabar)
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=1"       # photon from Sigmabar0 -> gamma Lambdabar
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprm "==1"       # exactly one anti-proton; no extra anti-proton in event
  end
  .remove([:prm <= :chrgn])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :electron, against: [:pion, :kaon, :proton]
    nep ">=1"
  end
  .remove([:ep <= :chrgp])
  .assign({chrgp: :pip, chrgn: :pim})
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:Lambda_bar, :gamma, :ep]) do
    invariant_mass_of(:Lambda_bar, :gamma).constrain_to_nominal_mass_of(:Sigma_bar0)  # 1C Sigmabar0 mass constraint
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

alg_sigmabar0_eplus.note(:positron_pid_combined_prob,
  "Paper requires CL_e/(CL_e+CL_pi+CL_K+CL_p) > 0.8 and 0.8 < E/p < 1.2 for positron. " \
  "These cuts are applied in ROOT analysis; only CL_e > 0.001 is expressed in the DSL PID block.")
alg_sigmabar0_eplus.note(:lambda_sigma0_mass_window,
  "Paper applies Lambdabar mass window (1.110, 1.121) GeV/c^2 and Sigmabar0 mass window (1.173, 1.200) GeV/c^2. Applied in ROOT analysis.")
alg_sigmabar0_eplus.note(:photon_conversion_veto,
  "Paper vetoes events where the positron forms a photon-conversion pair (Delta_xy, cos(theta_eg), R_xy cuts). Applied in ROOT analysis.")
alg_sigmabar0_eplus.note(:deltaE_mbc,
  "Paper applies DeltaE (-0.028, 0.024) GeV and MBC selection with ARGUS background fit. Applied in ROOT analysis.")

alg_sigmabar0_eplus.with_decay_card(decay_card_sigmabar0_eplus).apply(selection_sigmabar0_eplus)
root_files_sigmabar0_eplus = alg_sigmabar0_eplus.execute_on([psipp_data, psipp_incMC, exMC_sigmabar0_eplus])