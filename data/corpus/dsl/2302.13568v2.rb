# 2302.13568v2: Sigma+ -> p gamma at J/psi
# Conventional tag-and-probe (NOT TagAnalysis): ST Sigma_bar- -> anti-p pi0, DT Sigma+ -> p gamma
# 5C kinematic fit (4-momentum + ST pi0 mass constraint)
# Two charge-conjugate modes analyzed jointly

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 Sigma+ anti-Sigma-  PHSP;
  Enddecay
  Decay Sigma+
  1.0000 p+ gamma  PHSP;
  Enddecay
  Decay anti-Sigma-
  1.0000 anti-p- pi0  PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

decay_card_cp = <<~DECAYCARD
  Decay J/psi
  1.0000 Sigma+ anti-Sigma-  PHSP;
  Enddecay
  Decay Sigma+
  1.0000 p+ pi0  PHSP;
  Enddecay
  Decay anti-Sigma-
  1.0000 anti-p- gamma  PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_SigmaPlus_pgamma"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

algorithm = Algorithm.new("SigmaPlusToPGammaJpsi")
algorithm
  .set_header(["SigmaPlusToPGammaJpsiAlg/SigmaPlusToPGammaJpsi.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })
  .with_decay_card(decay_card_signal)

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=1"
    nChrn     ">=1"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=3"
  end
  .kinematic_fit([:prp, :prm, :gamma, :gamma, :gamma]) do
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 30
    nominal
  end
  .kinematic_fit([:prp, :prm, :gamma, :gamma, :gamma, :gamma]) do
    use_track_index_from_nominal_kmfit
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  end

algorithm
  .note(:st_approach, "ST: reconstruct anti-Sigma- -> anti-p pi0 via recoil mass M_rec; ST yield from fit to M_rec distribution. DT: Sigma+ -> p gamma searched in remaining tracks/photons")
  .note(:sigma_mass_window, "ST anti-Sigma- candidates: |M(anti-p pi0) - M(Sigma-)| < 13.5 MeV/c^2; best candidate chosen by minimizing mass difference")
  .note(:best_photon, "If multiple photon candidates on DT side, keep the one with minimum chi2_5C")
  .note(:decay_length_veto, "Secondary vertex fit on p for Sigma+ decay length: L/sigma_L > 1.5 to veto Delta(1232)+ -> p pi0 background; 93% background reduction for 78% signal retention")
  .note(:competing_hypothesis, "Veto Sigma+ -> p pi0 background: eliminate candidate if chi2(J/psi -> p anti-p pi0 gamma gamma) < chi2_5C(signal hypothesis)")
  .note(:charge_conjugate, "Charge-conjugate mode (anti-Sigma- -> anti-p gamma) analyzed with symmetric selection and simultaneous likelihood fit for branching fraction")
  .note(:p_momentum_signal_region, "DT signal region: 0.215 < P_p(Sigma+ rest frame) < 0.235 GeV/c; signal peak at 0.223 GeV/c")
  .note(:alpha_gamma_fit, "Decay asymmetry alpha_gamma determined from unbinned MLL fit to angular distribution; additional chi2 veto for purity")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])