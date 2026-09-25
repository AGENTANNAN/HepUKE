# Paper: 2310.00720v1 - Lambda + 9Be -> Sigma+ + X
# J/psi -> Lambda anti-Lambda, anti-Lambda -> anti-p pi+, Sigma+ -> p pi0, pi0 -> gamma gamma
# Single-tag (anti-Lambda) + Double-tag (Sigma+) method for hyperon-nucleon scattering

DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card_jpsi_ll = <<~DECAYCARD
  Decay J/psi
  1.0000 anti-Lambda-0 Lambda0 PHSP;
  Enddecay
  Decay anti-Lambda-0
  1.0000 anti-p- pi+ PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

signal_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_lambda_sigma"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_jpsi_ll
  config.cross_section   = :default
end

algorithm = Algorithm.new("LambdaSigmaX", "00-00-01")

# Single-tag: anti-Lambda -> anti-p pi+ (vertex fit); Pi0 -> gamma gamma (kalman fit)
# Double-tag: Sigma+ -> p pi0 in the recoil side
event_selection = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nChrp       ">=1"   # need p (from Sigma+) and pi+ (from anti-Lambda)
    nChrn       ">=1"   # need anti-p (from anti-Lambda)
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .kinematic_fit([:Lambda_bar, :pi0, :prp]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

algorithm
  .set_header(["LambdaSigmaXAlg/LambdaSigmaX.h"])
  .set_constant(ECMS: 3.097)
  .note(:single_tag_Lambda, "Single-tag: Lambda_bar reconstructed via secondary vertex fit (anti-p pi+). Vertex fit chi2 < 200, decay length > 0. Mass window on anti-p pi+ invariant mass [1.111, 1.120] GeV/c2. RM(anti-p pi+) window [1.071, 1.153] GeV/c2. Best candidate selected by minimum vertex-fit chi2. N_ST determined by fitting RM distribution with double Gaussian + Chebyshev polynomial.")
  .note(:double_tag_Sigma, "Double-tag: Sigma+ -> p pi0. pi0 reconstructed via 1C kinematic fit on gamma gamma combinations. Sigma+ selected from p pi0 with M(p pi0) in [1.12, 1.25] GeV/c2, using candidate with maximum proton PID likelihood. RM(anti-p pi+ p) required to be negative (suppresses non-scattering background). Lambda veto in double-tag side: if M(p pi-) in [1.108, 1.124] GeV/c2, the event is rejected. N_DT extracted by fitting M(p pi0) with double Gaussian + 3rd-order Chebyshev polynomial.")
  .note(:effective_luminosity, "Effective luminosity L_Lambda calculated using target geometry (beam pipe + MDC inner wall, 7 layers), average path lengths per layer from MC, and cross-section ratio R_sigma for different nuclei normalized to Be. L_Lambda = (17.00 +/- 0.01) x 10^28 cm^-2.")
  .note(:cross_section_calc, "Cross section sigma(Be) = N_DT / (epsilon_sig * L_Lambda * B(Sigma+ -> p pi0)). Sigma+ -> p pi0 BF = (51.57 +/- 0.30)%. Single-tag efficiency 52.16%, double-tag efficiency 24.32%.")
  .with_decay_card(decay_card_jpsi_ll)
  .apply(event_selection)

algorithm.execute_on([jpsi_data, jpsi_incMC, signal_mc])