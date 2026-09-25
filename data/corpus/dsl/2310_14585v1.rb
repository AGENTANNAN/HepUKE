# Paper: 2310.14585v1 - psi(3686) -> Sigma+ anti-Sigma- omega(phi)
# omega-mode: full reconstruction, Sigma+->p pi0, anti-Sigma-->anti-p pi0, omega->pi+ pi- pi0, 7C kinematic fit
# phi-mode: partial reconstruction, Sigma+->p pi0, anti-Sigma-->anti-p pi0, phi->K+ K-, one pi0 missing, 2C fit

DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for omega-mode: psi(3686) -> Sigma+ anti-Sigma- omega
decay_card_omega = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Sigma+ anti-Sigma- omega PHSP;
  Enddecay
  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay
  Decay anti-Sigma-
  1.0000 anti-p- pi0 PHSP;
  Enddecay
  Decay omega
  1.0000 pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Decay card for phi-mode: psi(3686) -> Sigma+ anti-Sigma- phi
decay_card_phi = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Sigma+ anti-Sigma- phi PHSP;
  Enddecay
  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay
  Decay anti-Sigma-
  1.0000 anti-p- pi0 PHSP;
  Enddecay
  Decay phi
  1.0000 K+ K- PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

signal_mc_omega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_SigmaSigma_omega"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_omega
  config.cross_section   = :default
end

signal_mc_phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_SigmaSigma_phi"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_phi
  config.cross_section   = :default
end

# ============================================================================
# omega-mode: psi(3686) -> Sigma+ anti-Sigma- omega, 7C kinematic fit
# Sigma+ -> p pi0, anti-Sigma- -> anti-p pi0, omega -> pi+ pi- pi0
# ============================================================================
algorithm_omega = Algorithm.new("SigmaSigmaOmega", "00-00-01")

selection_omega = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          2.0        # looser transverse for tracks from Sigma decay
    nChrp       ">=2"      # p (from Sigma+) and pi+ (from omega)
    nChrn       ">=2"      # anti-p (from anti-Sigma-) and pi- (from omega)
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
    identify :pion, against: [:kaon, :proton]
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=6"    # three pi0 -> six photons
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=3"
  end
  .kinematic_fit([:prp, :prm, :pip, :pim, :pi0, :pi0, :pi0]) do
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 45
    nominal
  end

algorithm_omega
  .set_header(["SigmaSigmaVAlg/SigmaSigmaVAlg.h"])
  .set_constant(ECMS: 3.686)
  .note(:mode, "omega-mode: psi(3686) -> Sigma+ anti-Sigma- omega. Sigma+ -> p pi0, anti-Sigma- -> anti-p pi0, omega -> pi+ pi- pi0. Three pi0s reconstructed from gamma gamma with 1C kinematic fit.")
  .note(:kinematic_fit, "7C kinematic fit: 4C (four-momentum) + 3x1C (pi0 mass constraints). chi2_7C < 45. Best combination by minimum chi2_7C.")
  .note(:Sigma_selection, "Sigma+ and anti-Sigma- selected by minimizing delta = sqrt((M(p pi0)-m_Sigma+)^2 + (M(anti-p pi0)-m_anti-Sigma-)^2 + (M(pi+ pi- pi0)-m_omega)^2). Sigma mass window [1176, 1197] MeV/c2.")
  .note(:veto, "J/psi veto: |RM(pi+ pi-) - m_J/psi| > 10 MeV/c2. |RM(pi+ pi0) - m_J/psi| > 15 MeV/c2. eta veto: |M(pi+ pi- pi0) - m_eta| > 13 MeV/c2. Lambda veto: |M(p pi-) - m_Lambda| > 6 MeV/c2, |M(anti-p pi+) - m_anti-Lambda| > 6 MeV/c2.")
  .note(:photon_count, "N(gamma)=6 selected via chi2_4C(6g p anti-p pi+ pi-) < chi2_4C(5g p anti-p pi+ pi-) and < chi2_4C(7g p anti-p pi+ pi-) to suppress wrong photon multiplicities.")
  .note(:result, "Observed with 13.8 sigma significance. BF = (1.90 +/- 0.18 +/- 0.21) x 10^-5.")
  .with_decay_card(decay_card_omega)
  .apply(selection_omega)

algorithm_omega.execute_on([psip_data, psip_incMC, signal_mc_omega])

# ============================================================================
# phi-mode: psi(3686) -> Sigma+ anti-Sigma- phi, partial reco (one pi0 missing)
# Sigma+ -> p pi0, anti-Sigma- -> anti-p pi0, phi -> K+ K-, one pi0 missing
# ============================================================================
algorithm_phi = Algorithm.new("SigmaSigmaPhi", "00-00-01")

selection_phi = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          2.0        # looser transverse for tracks from Sigma decay
    nChrp       ">=2"      # p (from Sigma+) and K+ (from phi)
    nChrn       ">=2"      # anti-p (from anti-Sigma-) and K- (from phi)
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
    identify :kaon, against: [:pion, :proton]
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=1"
  end
  .partial_miss([:pi0]) do
    best_combination_by_mass
  end
  .kinematic_fit([:prp, :prm, :kp, :km, :pi0]) do
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    nominal
  end

algorithm_phi
  .set_header(["SigmaSigmaVAlg/SigmaSigmaVAlg.h"])
  .set_constant(ECMS: 3.686)
  .note(:mode, "phi-mode: psi(3686) -> Sigma+ anti-Sigma- phi. Partial reconstruction: only one pi0 reconstructed (Sigma+ -> p pi0), the other pi0 (from anti-Sigma- -> anti-p pi0) treated as missing.")
  .note(:kinematic_fit, "2C kinematic fit on p anti-p K+ K- pi0 with missing mass constrained to nominal pi0 mass. chi2_2C < 20. Best combination by minimum chi2_2C.")
  .note(:Sigma_selection, "Sigma+ and anti-Sigma- selected by minimizing delta = sqrt((M(p pi0)-m_Sigma+)^2 + (M(anti-p pi0)-m_anti-Sigma-)^2). Sigma mass window [1176, 1200] MeV/c2.")
  .note(:veto, "J/psi veto: |RM(pi+ pi0) - m_J/psi| > 9 MeV/c2 (applied in ROOT).")
  .note(:partial_reco, "One pi0 missing. 2C fit constrains missing mass to pi0 mass and gamma gamma invariant mass to pi0 mass.")
  .note(:result, "Observed with 7.6 sigma significance. BF = (2.96 +/- 0.54 +/- 0.41) x 10^-6.")
  .note(:sideband, "2D Sigma+ anti-Sigma- sideband method for background estimation: signal region with both M(p pi0) and M(anti-p pi0) in Sigma signal window; sideband I with one in signal and one in sideband; sideband II with both in sidebands.")
  .with_decay_card(decay_card_phi)
  .apply(selection_phi)

algorithm_phi.execute_on([psip_data, psip_incMC, signal_mc_phi])