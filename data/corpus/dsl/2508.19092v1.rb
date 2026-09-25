# 2508.19092v1: First observation of psi(3686) -> omega eta eta
# omega -> pi+ pi- pi0, eta -> gamma gamma (two eta)
# Final state: pi+ pi- pi0 2eta -> pi+ pi- 6gamma
# psi(2S) data (BOSS 709_3686), 2.712 x 10^9 events

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.000 omega eta eta PHSP;
    Enddecay

    Alias another_eta eta

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay another_eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psi2S_omega_eta_eta_exclusive_mc"
  config.related_dataset = psip_data
  config.events = 200000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection ###
alg = Algorithm.new("Psi2SToOmegaEtaEta")
alg.set_header(["Psi2SToOmegaEtaEtaAlg/Psi2SToOmegaEtaEta.h"])
   .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==1"
  nChrn "==1"
  nNet "==0"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=6"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :pion, against: [:kaon, :proton]
  npip "==1"
  npim "==1"
end
# Reconstruct pi0 from photon pairs (1C Kalman fit)
.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0 ">=1"
end
# Reconstruct two eta from photon pairs (1C Kalman fit)
.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 200
  neta ">=2"
end
# 4C kinematic fit: psi(3686) -> pi+ pi- pi0 eta eta
.kinematic_fit([:pip, :pim, :pi0, :eta, :eta]) do
  nominal
  constrain_four_momentum
  chi2_cut 40
end
.note(:pi0_eta_selection,
  "pi0 and eta combinations selected by minimizing chi2_eta_eta_pi0 combining mass deviations. " \
  "pi0 mass window: |M(gamma gamma) - m_pi0| < 20 MeV/c^2. " \
  "eta mass windows: |M(gamma gamma) - m_eta| < 25 MeV/c^2 for both eta candidates.")
.note(:competing_hypothesis,
  "Additional 4C fits under psi(3686) -> pi+ pi- 7gamma and pi+ pi- 8gamma hypotheses. " \
  "Events with chi2_4C smaller than signal hypothesis are discarded. " \
  "chi2 comparisons: chi2_eta_eta_pi0 < chi2_pi0_pi0_eta, " \
  "chi2_eta_eta_pi0 < chi2_pi0_pi0_pi0, chi2_eta_eta_pi0 < chi2_eta_eta_eta.")
.note(:background_veto,
  "Two-pi0 veto: |M(gamma gamma) - m_pi0| > 0.03 GeV/c^2 for non-pi0 photon pairs. " \
  "psi(3686) -> pi pi J/psi veto: |M(pi pi)_recoil - M_J/psi| > 0.02 GeV/c^2. " \
  "psi(3686) -> X J/psi, J/psi -> omega eta veto: M(omega eta) < 3.0 GeV/c^2.")
.note(:signal_extraction,
  "omega signal from M(pi+ pi- pi0) distribution fit: " \
  "MC template convolved with Gaussian for signal, 2nd-order Chebyshev for combinatorial background, " \
  "non-eta background from 2D eta_H-eta_L sideband. " \
  "Continuum background estimated from psi(3770) data scaled to psi(2S) energy.")
.note(:interference,
  "Interference with continuum production evaluated from psi(3770) data (709_3773), " \
  "1/s cross-section scaling applied. Phase angle phi = +/- 90 deg for conservative estimate.")

alg.with_decay_card(decay_card).apply(event_selection)
alg.execute_on([psip_data, psip_incMC, exMC_signal])