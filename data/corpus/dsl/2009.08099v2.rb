# Paper: 2009.08099v2 -- Resonant structures in e+e- -> omega eta and e+e- -> omega pi0 at 2.00-3.08 GeV
# BESIII energy scan at 22 c.m. energies. NOT a tag-based analysis.

# Energy points: 2000, 2050, 2100, 2125, 2150, 2175, 2200, 2232, 2309,
#               2386, 2396, 2500, 2644, 2646, 2700, 2800, 2900, 2950,
#               2981, 3000, 3020, 3080
rscan_energies = %w[
  2000 2050 2100 2125 2150 2175 2200 2232 2309
  2386 2396 2500 2644 2646 2700 2800 2900 2950
  2981 3000 3020 3080
]

rscan_datasets = rscan_energies.map { |e| DatasetManager.real_data.find("713_Rscan_#{e}") }
rscan_incMC   = rscan_energies.map { |e| DatasetManager.inclusive_mc.find("713_Rscan_#{e}") }

# ==============================
# Algorithm 1: e+e- -> omega eta, omega -> pi+ pi- pi0, pi0 -> gamma gamma, eta -> gamma gamma
# ==============================
decay_card_omega_eta = <<~DECAYCARD
    Decay vpho
    1.0000 omega eta PHSP;
    Enddecay
    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay
    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
End
DECAYCARD

exMC_omega_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "omega_eta_2125_exclusive_mc"
  config.related_dataset = DatasetManager.real_data.find("713_Rscan_2125")
  config.events = 200000
  config.decay_card = decay_card_omega_eta
  config.cross_section = :default
end

alg_omega_eta = Algorithm.new("OmegaEta", version: '00-00-01')
alg_omega_eta.set_header(["OmegaEtaAlg/OmegaEta.h"])
             .set_constant({"ECMS" => [:double, 2.125]})
             .note(:analysis, "e+e- -> omega eta at 22 c.m. energies 2.000-3.080 GeV. omega -> pi+pi-pi0, pi0->gamma gamma, eta->gamma gamma. Final state: pi+ pi- 4gamma.")
             .note(:track_selection, "Charged tracks: |cos theta| < 0.93, Vr < 1 cm, Vz < 10 cm. Exactly two oppositely charged pions. PID: highest likelihood for pi hypothesis.")
             .note(:photon_selection, "E_gamma > 25 MeV (barrel, |cos theta|<0.80), > 50 MeV (endcap, 0.86<|cos theta|<0.92). EMC time [0,700] ns.")
             .note(:kinematic_fit, "4C kinematic fit: e+e- -> pi+ pi- 4gamma. chi2_4C < 70. If more than 4 photons, combination with smallest chi2 retained. Veto: chi2_4C(pi+pi-4gamma) < chi2_4C(pi+pi-5gamma) to suppress omega pi0 pi0.")
             .note(:pi0_eta_pairing, "Best pi0-eta pair selected by minimizing chi2_alpha_beta = (M(g1g2)-m_alpha)^2/sigma12^2 + (M(g3g4)-m_beta)^2/sigma34^2. Require chi2_pi0eta < chi2_pi0pi0 and chi2_pi0eta < chi2_etaeta.")
             .note(:mass_windows, "|M(gammagamma)-m_pi0| < 0.020 GeV/c2; |M(gammagamma)-m_eta| < 0.030 GeV/c2. |E_gamma3 - E_gamma4|/p_eta < 0.9 veto for omega gamma_ISR and omega pi0 pi0 suppression.")
             .note(:signal_extraction, "Simultaneous unbinned ML fit to M(pi+pi-pi0) in eta signal and sideband regions. Signal: MC shape conv. Gaussian. Background: 2nd-order Chebychev. Peaking background from omega pi0 pi0 subtracted via eta sideband scaling.")
             .note(:cross_section, "sigma = N_sig / (L * epsilon * (1+delta) * B), B = B(omega->pi+pi-pi0) * B(pi0->gamma gamma) * B(eta->gamma gamma) = 34.7%")
             .with_decay_card(decay_card_omega_eta)

event_selection_omega_eta = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
  end
  .pid(method: :probability) do
    identify :pion, against: [:kaon, :proton]
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=4"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  end
  .kinematic_fit([:pip, :pim, :pi0, :eta]) do
    constrain_four_momentum
    chi2_cut 70
    nominal
  end

alg_omega_eta.with_decay_card(decay_card_omega_eta).apply(event_selection_omega_eta)
alg_omega_eta.execute_on(rscan_datasets + rscan_incMC + [exMC_omega_eta])

# ==============================
# Algorithm 2: e+e- -> omega pi0, omega -> pi+ pi- pi0, pi0 -> gamma gamma (x2)
# ==============================
decay_card_omega_pi0 = <<~DECAYCARD
    Decay vpho
    1.0000 omega pi0 PHSP;
    Enddecay
    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
End
DECAYCARD

exMC_omega_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "omega_pi0_2125_exclusive_mc"
  config.related_dataset = DatasetManager.real_data.find("713_Rscan_2125")
  config.events = 200000
  config.decay_card = decay_card_omega_pi0
  config.cross_section = :default
end

alg_omega_pi0 = Algorithm.new("OmegaPi0", version: '00-00-01')
alg_omega_pi0.set_header(["OmegaPi0Alg/OmegaPi0.h"])
             .set_constant({"ECMS" => [:double, 2.125]})
             .note(:analysis, "e+e- -> omega pi0 at 22 c.m. energies 2.000-3.080 GeV. omega -> pi+pi-pi0, both pi0 -> gamma gamma. Final state: pi+ pi- 4gamma.")
             .note(:pi0_pairing, "Two pi0 candidates. pi0-pi0 pair selection via chi2_pi0pi0 minimization. The pi0 whose pi+pi-pi0 invariant mass is closest to m_omega is the omega daughter (pi_omega0); the other is the bachelor pi0 (pi_bach0).")
             .note(:pi0_mass_window, "|M(gammagamma) - m_pi0| < 0.020 GeV/c2 for both pi0 candidates.")
             .note(:background, "Dominant background: e+e- -> pi+pi-pi0pi0 (same final state). Peaking background from mis-combination negligible.")
             .note(:signal_extraction, "Unbinned ML fit to M(pi+pi-pi_omega0). Signal: MC shape conv. Gaussian. Background: 2nd-order Chebychev. One-dimensional fit (peaking backgrounds negligible).")
             .note(:cross_section, "sigma = N_sig / (L * epsilon * (1+delta) * B), B = B(omega->pi+pi-pi0) * B^2(pi0->gamma gamma) = 87.1%")
             .with_decay_card(decay_card_omega_pi0)

event_selection_omega_pi0 = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
  end
  .pid(method: :probability) do
    identify :pion, against: [:kaon, :proton]
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=4"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  end
  .kinematic_fit([:pip, :pim, :pi0]) do
    constrain_four_momentum
    chi2_cut 70
    nominal
  end

alg_omega_pi0.with_decay_card(decay_card_omega_pi0).apply(event_selection_omega_pi0)
alg_omega_pi0.execute_on(rscan_datasets + rscan_incMC + [exMC_omega_pi0])