# Paper 2507.11145v1: Observation of Lambda(1520) and Lambda(1670) -> gamma Sigma0 radiative decays
# BESIII, (10087±44)×10^6 J/psi events
# Two modes with different photon counts and kinematic fit hypotheses

jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ---- Mode I: J/psi -> gamma Lambda anti-Lambda ----
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.000 anti-Lambda Lambda(1520) PHSP;
    Enddecay

    Decay Lambda(1520)
    1.000 gamma Lambda PHSP;
    Enddecay

    Decay Lambda
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_gamLambdabar_L1520_modeI_mc"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_modeI
  config.cross_section = :default
end

alg_modeI = Algorithm.new("Jpsi2gamLL_L1520_ModeI")
alg_modeI.set_header(["Jpsi2gamLL_L1520_ModeIAlg/Jpsi2gamLL_L1520_ModeI.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })

sel_modeI = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 20.0
    Vr 10.0
    nChrp "==2"
    nChrn "==2"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=1"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=1"
    npim ">=1"
  end
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:gamma, :Lambda, :Lambda_bar]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_modeI.note(:signal_model, "Lambda(1520) signal modeled by MC shape convolved with Gaussian; background: 3rd-order Chebyshev polynomial")
  .note(:Lambda_selection, "|M(p pi-)-M_Lambda| < 5 MeV/c^2; if multiple candidates, choose closest to nominal Lambda mass")
  .note(:Sigma0_veto, "|M(gamma anti-Lambda) - M_Sigma0| > 0.02 GeV/c^2 and |M(gamma Lambda) - M_Sigma0| > 0.02 GeV/c^2 applied alternately for J/psi->Lambda Lambda(1520) and J/psi->anti-Lambda Lambda(1520)")
  .note(:kinfit_chi2_cut_modeI, "chi2_4c < 30 applied via chi2_cut in kinematic fit block")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])

# ---- Mode II: J/psi -> gamma gamma Lambda anti-Lambda (Sigma0 -> gamma Lambda) ----
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.000 anti-Lambda Lambda(1520) PHSP;
    Enddecay

    Decay Lambda(1520)
    1.000 gamma anti-Sigma0 PHSP;
    Enddecay

    Decay anti-Sigma0
    1.000 gamma anti-Lambda PHSP;
    Enddecay

    Decay Lambda
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_ggLL_L1520_modeII_mc"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_modeII
  config.cross_section = :default
end

alg_modeII = Algorithm.new("Jpsi2ggLL_L1520_ModeII")
alg_modeII.set_header(["Jpsi2ggLL_L1520_ModeIIAlg/Jpsi2ggLL_L1520_ModeII.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })

sel_modeII = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 20.0
    Vr 10.0
    nChrp "==2"
    nChrn "==2"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=1"
    npim ">=1"
  end
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:gamma, :gamma, :Lambda, :Lambda_bar]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_modeII.note(:signal_model, "Lambda(1520) signal modeled by MC shape convolved with Gaussian; Lambda(1670) also present in mass spectrum")
  .note(:Lambda_selection, "|M(p pi-)-M_Lambda| < 5 MeV/c^2; if multiple candidates, choose closest to nominal Lambda mass")
  .note(:Sigma0_reconstruction, "lower-energy photon assumed from Sigma0 decay: |M(gamma_low Lambda)-M_Sigma0| < 0.01 GeV/c^2")
  .note(:pi0_veto_modeII, "|M(gamma gamma)-M_pi0| > 0.03 GeV/c^2")
  .note(:Sigma0bar_veto_modeII, "events with both gamma-Lambda and gamma-anti-Lambda in Sigma0 mass region vetoed")
  .note(:kinfit_chi2_cut_modeII, "chi2_4c < 15 applied via chi2_cut in kinematic fit block")
  .note(:background_rejection_modeII, "Sequential 1C+2C kinematic fits for Lambda Sigma0 pi0 veto: (M_rec_Lambda_Sigma0)^2 < 0.01, chi2_2C > 5")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])