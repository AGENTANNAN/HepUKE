# BOSS DSL for 2401.00918v2: J/psi → gamma gamma phi PWA
# Data: (10087±44)×10^6 J/psi events at 3.097 GeV
# Final state: gamma gamma K+ K- (phi → K+K-)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card: J/psi → gamma gamma phi, phi → K+ K-
# PHSP generator for signal MC (uniform phase space for PWA)
decay_card = <<~DECAYCARD
  Decay jpsi
  1.0000 gamma gamma phi PHSP;
  Enddecay
  Decay phi
  1.0000 K+ K- VSS;
  Enddecay
  End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_gamma_phi"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# Event selection
event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
  end
  .pid(method: :probability) do
    identify :kaon, against: [:pion, :proton]
    nkp "==1"; nkm "==1"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  # Nominal 4C kinematic fit: gamma gamma K+ K-
  .kinematic_fit([:gamma, :gamma, :kp, :km]) do
    constrain_four_momentum
    chi2_cut 40
    nominal
  end
  # Competing hypothesis: 1-gamma 4C fit
  .kinematic_fit([:gamma, :kp, :km]) do
    use_track_index_from_nominal_kmfit
    constrain_four_momentum
  end
  # Competing hypothesis: 3-gamma 4C fit
  .kinematic_fit([:gamma, :gamma, :gamma, :kp, :km]) do
    use_track_index_from_nominal_kmfit
    constrain_four_momentum
  end
  # Competing hypothesis: 4-gamma 4C fit
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :kp, :km]) do
    use_track_index_from_nominal_kmfit
    constrain_four_momentum
  end

algorithm = Algorithm.new("JpsiGammaGammaPhi")
algorithm
  .set_header(["JpsiGammaGammaPhi/JpsiGammaGammaPhi.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })
  .note(:competing_hypothesis_veto,
    "Event retained only if chi2_4C(gamma gamma K+K-) < chi2_4C of any additional hypothesis " \
    "(gamma K+K-, 3gamma K+K-, 4gamma K+K-). Implemented as multiple kinematic_fit blocks " \
    "with stored chi2 values compared in ROOT analysis.")
  .note(:phi_mass_window,
    "Post-fit: phi candidates identified by |M(K+K-) - M_phi| < 0.005 GeV/c^2.")
  .note(:pi0_eta_etap_veto,
    "Post-fit: events with M(gamma gamma) in (0.11, 0.16), (0.46, 0.59), or (0.92, 0.99) GeV/c^2 " \
    "excluded to suppress pi0/eta/etap backgrounds.")
  .note(:dalitz_cut,
    "Post-fit: Dalitz cut M^2(gamma_high phi) < 9 GeV^2/c^4 to suppress residual phi-pi0/eta/etap backgrounds.")
  .note(:q_factor_method,
    "Q-factor (quality factor) method for non-phi background subtraction: " \
    "7 kinematic coordinates, n_c=200 nearest neighbors. Performed in ROOT analysis stage.")
  .note(:multidimensional_reweighting,
    "Multi-dimensional reweighting of phi pi0 pi0 PHSP MC for phi-related background subtraction. " \
    "Performed in ROOT analysis stage.")
  .note(:pwa,
    "Partial-wave analysis using GPUPWA framework. PWA is performed in ROOT analysis stage " \
    "with covariant tensor amplitudes; resonances with JPC=0++, 0-+, 1++, 1-+, 2++, 2-+ considered. " \
    "Not expressible at BOSS DSL level.")
  .with_decay_card(decay_card)
  .apply(event_selection)
  .execute_on([jpsi_data, jpsi_incMC, exMC])