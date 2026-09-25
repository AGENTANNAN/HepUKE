DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

data_4600 = DatasetManager.real_data.find("703_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")
all_data = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
all_incMC = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0 D_s*+ D_sJ- PHSP;
  Enddecay
  Decay D_s*+
  1.0 gamma D_s+ VSP_PWAVE 1.0 2.0 0.0;
  Enddecay
  Decay D_s+
  1.0 K+ K- pi+ PHSP;
  Enddecay
  Decay D_sJ-
  1.0 pi+ pi- pi0 PHSP;
  Enddecay
  End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "sig_DsStarDsJ"
  config.events        = 500_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

algorithm = Algorithm.new("DsStarDsJ_XS")
algorithm
  .set_header(["DsStarDsJ_XS/DsStarDsJ_XS.h"])
  .note(:partial_reconstruction,
    "DsJ- (Ds0*(2317)-, Ds1(2460)-, Ds1(2536)-) is reconstructed via the recoil mass of Ds*+. " \
    "The improved recoil mass is defined as RM(Ds*+) = M_recoil(gamma K+ K- pi+) + M(gamma K+ K- pi+) - m_Ds*+. " \
    "DsJ- signal yields are extracted from unbinned maximum-likelihood fits to the recoil mass distributions in ROOT, " \
    "with signal shape from MC simulation and ARGUS function for background. " \
    "Born cross sections are calculated as sigma_B = N_fit / (L_int * (1+delta) * (1+delta_vp) * epsilon_Ds*+).")
  .note(:kinematic_fit_2C,
    "A 2C kinematic fit is applied to gamma K+ K- pi+ constraining the masses to the nominal Ds+ (1.96834 GeV/c^2) " \
    "and Ds*+ (2.1122 GeV/c^2) masses. The chi2 of this fit is required to be less than 10 to suppress backgrounds. " \
    "This fit is performed before the partial reconstruction in the BOSS implementation.")
  .note(:submode_selection,
    "Ds+ -> phi pi+ submode: |M(K+K-) - m_phi| < 9 MeV/c^2. " \
    "Ds+ -> anti-K*0 K+ submode: |M(K-pi+) - m_anti-K*0| < 84 MeV/c^2. " \
    "Events must satisfy at least one submode. Both submodes contribute to the Ds*+ reconstruction.")
  .note(:background,
    "Generic MC samples of open-charm processes used to estimate background contributions. " \
    "No peaking backgrounds found in the signal region from generic MC studies using TopoAna. " \
    "Ds*+ mass sidebands used to verify absence of peaking structures from processes such as " \
    "e+e- -> Ds+ gamma DsJ- and e+e- -> D0 (K-pi+) K+ gamma Ds0*(2317)-.")
  .note(:cross_section_measurement,
    "Born cross sections measured for three DsJ- states at 7 energy points (4.600-4.700 GeV). " \
    "KKMC generator used for signal MC with ISR and beam-energy spread effects. " \
    "ISR radiative correction factor (1+delta) obtained from QED calculation with 1% accuracy. " \
    "Vacuum polarization factor delta_vp ~0.055 for all energy points. " \
    "Cross section line shape changed to first-order polynomial multiplied by sqrt(E_m - E_0) " \
    "for systematic uncertainty estimation. Bayesian upper limits at 90% C.L. set for energy points " \
    "with signal significance < 3sigma.")
  .with_decay_card(decay_card)

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==1"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    20.0
    nGam              ">=1"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion]
    identify :pip, against: [:kaon]
    nkp  "==1"
    nkm  "==1"
    npip "==1"
  end
  .partial_miss([2]) do
    best_combination_by_mass :D_s_star_p, 2.1122
    require_recoil_mass 2.0, 2.8
  end

algorithm.apply(event_selection)
algorithm.execute_on(all_data + all_incMC + sig_mc)