# BESIII DSL: Amplitude analysis of Ds+ -> pi+ pi- pi+
# arXiv: 2108.10050v2
# e+e- -> Ds*+- Ds-+ at sqrt(s) = 4.178 GeV, Ds+ -> pi+ pi- pi+

### Dataset preparation ###
ds_data = DatasetManager.real_data.find("703_4180")
ds_incMC = DatasetManager.inclusive_mc.find("703_4180")

decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s-           PHSP;
    Enddecay

    Decay D_s*+
    1.000 D_s+ gamma           VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.000 pi+ pi- pi+          PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi-            PHSP;
    Enddecay

    End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "ds_pipipi_amp_exclusive_mc"
  config.related_dataset = ds_data
  config.events = 500_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg = Algorithm.new("DsPiPiPiAmp")
alg.set_header(["DsPiPiPiAmpAlg/DsPiPiPiAmp.h"])
alg.set_constant({"ECMS" => [:double, 4.178]})

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=1"
    nChrn ">=1"
    nTot "==3"
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
    identify :pion, against: [:kaon, :proton]
    npip ">=1"
    npim ">=1"
  end
  # Ds+ -> pi+ pi- pi+ kinematic fit with Ds+ mass constraint
  .kinematic_fit([:pip, :pip, :pim]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg.note(:ks_veto, "K_S0 veto: any pion from a reconstructed K_S0->pi+pi- is rejected; K_S0 mass window +/-12 MeV/c^2, decay length > 2x vertex resolution")
  .note(:pi0_veto, "Photon from pi0->gamma gamma rejected for Ds* tagging; pi0 mass window 0.125-0.145 GeV/c^2")
  .note(:ds_star_tag, "D_s*+ -> D_s+ gamma: direct candidates via |M_rec - m(D_s*+)| < 0.02 GeV/c^2; indirect candidates via DeltaM = M(Ds+gamma)-M(Ds+) in [0.135, 0.15] GeV/c^2")
  .note(:nn_classifier, "TMVA multilayer perceptron NN for background suppression; separate classifiers for direct and indirect Ds+ categories with different input variables")
  .note(:ds_mass_window, "Ds+ mass window [1.900, 2.0535] GeV/c^2; signal region |m(pi+pi-pi+) - m(Ds+)| < 12 MeV/c^2; signal purity ~80.6%")
  .note(:mass_constrained_fit, "Kinematic fit with Ds+ mass constraint applied; corrected four-momenta used for Dalitz-plot amplitude analysis")
  .note(:amplitude_analysis, "QMIPWA approach: 29 S-wave control points with cubic spline interpolation + Breit-Wigner P/D waves (rho(770), rho(1450), f2(1270)); GooFit GPU framework; 62 free parameters")
  .with_decay_card(decay_card)
  .apply(event_selection)

root_files = alg.execute_on([ds_data, ds_incMC, exMC])