# BESIII Analysis: psi(3686) → Lambda Lambdabar eta'
# Paper: 2309.13790v1
# Data: psi(2S) at 3.686 GeV (BOSS 709) + psi(3770) continuum (BOSS 712)

# ============================================================
# Datasets
# ============================================================
psip_data   = DatasetManager.real_data.find("709_3686")
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")
psipp_data  = DatasetManager.real_data.find("712_3773")
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================
# Decay Card (shared by both modes up to eta' decay)
# ============================================================

# ============================================================
# Mode I: psi(3686) → Lambda Lambdabar eta', eta' → gamma pi+ pi-
# Final state: p pbar pi+ pi- pi+ pi- gamma
# ============================================================
decay_card_LambdaLambdabar_etap_modeI = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Lambda anti-Lambda- eta'  PHSP;
  Enddecay
  Decay eta'
  1.0000 gamma pi+ pi-  PHSP;
  Enddecay
  Decay Lambda
  1.0000 p+ pi-  PHSP;
  Enddecay
  Decay anti-Lambda-
  1.0000 anti-p- pi+  PHSP;
  Enddecay
  End
DECAYCARD

signal_mc_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_LambdaLambdabar_etap_gammapipi"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_LambdaLambdabar_etap_modeI
  config.cross_section   = :default
end

alg_I = Algorithm.new("LambdaLambdabarEtap_ModeI")
alg_I
  .set_header(["LambdaLambdabarEtap_ModeI/LambdaLambdabarEtap_ModeI.h"])

sel_I = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=3"
    nChrn ">=3"
    nTot ">5"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
    nprp ">=1"
    nprm ">=1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=1"
  end
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :pip, :pim, :pip, :pim]) do
    constrain_four_momentum
    chi2_cut 30
    nominal
  end

alg_I
  .note(:best_lambda_pair, "Best Lambda-Lambdabar pair selected by minimizing |M(p pi-) - m_Lambda| + |M(pbar pi+) - m_Lambda|")
  .note(:primary_vertex, "Primary vertex fit applied to remaining pi+ pi- tracks together with Lambda and Lambdabar")
  .note(:sigma0_veto, "Sigma0 veto: |M(gamma Lambda/Lambdabar) - m_Sigma0| > 10 MeV")
  .note(:jpsi_veto, "J/psi veto: |M(pi+ pi- gamma) - m_Jpsi| > 20 MeV for events where pi+pi- are from eta' decay")
  .note(:chi_c_veto, "chi_c vetoes applied")
  .note(:xi_veto, "Xi vetoes applied")
  .note(:sigma1385_veto, "Sigma(1385) vetoes applied")
  .note(:continuum_suppression, "psi(3770) data used for QED continuum background estimation")
  .with_decay_card(decay_card_LambdaLambdabar_etap_modeI)
  .apply(sel_I)

alg_I.execute_on([psip_data, psip_incMC, signal_mc_modeI])

# ============================================================
# Mode II: psi(3686) → Lambda Lambdabar eta', eta' → eta pi+ pi-, eta → gamma gamma
# Final state: p pbar pi+ pi- pi+ pi- gamma gamma
# ============================================================
decay_card_LambdaLambdabar_etap_modeII = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Lambda anti-Lambda- eta'  PHSP;
  Enddecay
  Decay eta'
  1.0000 eta pi+ pi-  PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma  PHSP;
  Enddecay
  Decay Lambda
  1.0000 p+ pi-  PHSP;
  Enddecay
  Decay anti-Lambda-
  1.0000 anti-p- pi+  PHSP;
  Enddecay
  End
DECAYCARD

signal_mc_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_LambdaLambdabar_etap_etapipi"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_LambdaLambdabar_etap_modeII
  config.cross_section   = :default
end

alg_II = Algorithm.new("LambdaLambdabarEtap_ModeII")
alg_II
  .set_header(["LambdaLambdabarEtap_ModeII/LambdaLambdabarEtap_ModeII.h"])

sel_II = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=3"
    nChrn ">=3"
    nTot ">5"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
    nprp ">=1"
    nprm ">=1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=2"
  end
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 40
    neta ">=1"
  end
  .kinematic_fit([:Lambda, :Lambda_bar, :eta, :pip, :pim, :pip, :pim]) do
    invariant_mass_of(:pip, :pim, :eta).constrain_to_nominal_mass_of("eta'")
    constrain_four_momentum
    chi2_cut 40
    nominal
  end

alg_II
  .note(:best_lambda_pair, "Best Lambda-Lambdabar pair selected by minimizing mass difference sum")
  .note(:primary_vertex, "Primary vertex fit applied to remaining pi+ pi- tracks together with Lambda and Lambdabar")
  .note(:sigma0_veto, "Sigma0 veto applied")
  .note(:jpsi_veto, "J/psi veto applied")
  .note(:chi_c_veto, "chi_c vetoes applied")
  .note(:xi_veto, "Xi vetoes applied")
  .note(:sigma1385_veto, "Sigma(1385) vetoes applied")
  .note(:eta_mass_window, "eta mass window applied in ROOT after kalman fit")
  .note(:continuum_suppression, "psi(3770) data used for QED continuum background estimation")
  .with_decay_card(decay_card_LambdaLambdabar_etap_modeII)
  .apply(sel_II)

alg_II.execute_on([psip_data, psip_incMC, signal_mc_modeII])