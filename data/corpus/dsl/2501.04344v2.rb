# Paper 2501.04344v2 — Study of electromagnetic Dalitz decay J/ψ → e⁺ e⁻ π⁰
# with (10087±44)×10⁶ J/ψ events
# Final state: e⁺ e⁻ γ γ (π⁰ → γγ)

### Dataset preparation ###
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card = <<~DECAYCARD
  Decay J/psi
  1.0000 e+ e- pi0 PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Signal exclusive MC
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_jpsi_ee_pi0"
  config.related_dataset = jpsi_data
  config.events = 500_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
algorithm = Algorithm.new("JpsiEEPizero")

algorithm.set_header(["JpsiEEPizeroAlg/JpsiEEPizero.h"])
          .set_constant("ECMS" => [:double, 3.097])

event_selection = Selection.new
  # Charged track selection: 2 oppositely charged tracks
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==1"
    nChrn "==1"
    nNet "==0"
  end
  # Photon selection: at least 2 photons
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  end
  # Lepton identification (high-momentum e±)
  .pid(method: :probability) do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  end
  # Reconstruct π⁰ → γγ via Kalman kinematic fit (1-C mass constraint)
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  # 4C kinematic fit with vertex constraint on e⁺e⁻ pair
  .kinematic_fit([:lp, :lm, :pi0]) do
    nominal
    vertex_fit([0, 1])
    constrain_four_momentum
    chi2_cut 200
  end

algorithm
  .note(:ep_over_p_cut, "E/p_{e±} > 0.8c required for tracks with momentum > 0.25 GeV/c to suppress pion contamination; applied after PID likelihood selection L(e±)>L(π±) and L(e±)>L(K±)")
  .note(:m_gammagamma_window, "diphoton invariant mass m_γγ required to be in [0.09, 0.18] GeV/c² for π⁰ signal region")
  .note(:gamma_conversion_veto, "δ_xy < 2 cm requirement suppresses J/ψ→γπ⁰ where photon converts to e⁺e⁻ in detector material; δ_xy = sqrt(R_x²+R_y²), eliminates ~98% conversion background with ~20% signal loss")
  .note(:two_photon_veto, "cosθ(e⁺) < 0.8 and cosθ(e⁻) > -0.8 to suppress two-photon process e⁺e⁻→e⁺e⁻π⁰ background; events from two-photon process accumulate at cosθ(e⁺) > 0.8 or cosθ(e⁻) < -0.8")
  .note(:electron_momentum_cut, "p_{e±} < 1.45 GeV/c to suppress radiative Bhabha background e⁺e⁻→γe⁺e⁻, especially at high m_{e⁺e⁻}")
  .note(:low_photon_energy_cut, "E_{γ_low} > 0.14 GeV for the lower-energy photon used in π⁰ reconstruction, to suppress non-peaking background")
  .note(:chi2_tight_cut, "4C kinematic fit χ² < 100 applied to suppress peaking background from J/ψ→γπ⁰π⁰ where one π⁰ decays via Dalitz mode π⁰→γe⁺e⁻")
  .note(:m_ee_binning, "signal yield extracted in 30 m_{e⁺e⁻} bins via extended ML fit to m_γγ; differential branching fraction and TFF |F(q²)|² extracted per bin")
  .note(:signal_extraction, "signal yield from extended ML fit to m_γγ: signal PDF = MC shape ⊗ Gaussian (free parameters, data-MC resolution); non-peaking background = 1st-order Chebyshev; peaking backgrounds (J/ψ→π⁺π⁻π⁰, γπ⁰, γπ⁰π⁰, ωπ⁰, two-photon) subtracted")
  .with_decay_card(decay_card)
  .apply(event_selection)

root_files = algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])