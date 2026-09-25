# 2508.18601v1: Search for chi_c1 -> pi+ pi- eta_c via psi(3686) -> gamma chi_c1
# eta_c reconstructed in 16 exclusive decay modes (see Table 1)
# psi(2S) data (BOSS 709_3686), 2712.4 x 10^6 events

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 HELAMP 1.0 0.0 -0.333 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
    Enddecay

    Decay chi_c1
    1.000 pi+ pi- eta_c PHSP;
    Enddecay

    Decay eta_c
    1.000 p+ anti-p- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "chi_c1_pipi_etac_exclusive_mc"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection: representative charged-mode channel (ppbar as example) ###
alg = Algorithm.new("Chic1ToPiPiEtaC_Charged")
alg.set_header(["Chic1ToPiPiEtaC_ChargedAlg/Chic1ToPiPiEtaC_Charged.h"])
   .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=1"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon, :proton]
  nprp ">=1"
  nprm ">=1"
end
.remove([:prp <= :chrgp, :prm <= :chrgn])
.assign({chrgp: :pip, chrgn: :pim})
# Soft pion momentum: p_pi < 0.4 GeV/c
# (This is implemented as a post-selection filter)
.for_each(:pip) do
  where { p < 0.4 }
  remove
end
.for_each(:pim) do
  where { p < 0.4 }
  remove
end
# 4C kinematic fit: psi(3686) -> gamma pi+ pi- p anti-p
# RM(gamma pi+ pi-) window [2.80, 3.20] via invariant_mass_of in fit
.kinematic_fit([:gamma, :pip, :pim, :prp, :prm]) do
  nominal
  constrain_four_momentum
  invariant_mass_of(:gamma, :pip, :pim).between(2.80, 3.20)
  chi2_cut 42
end
.note(:eta_c_modes,
  "eta_c reconstructed in 16 exclusive decay modes (Table 1): " \
  "ppbar, K+K-pi+pi-, ppbar pi+pi-, 4pi, 4K, 6pi, K+K-4pi, " \
  "K_S K pi, K_S K 3pi, K+K-pi0, ppbar pi0, K+K-eta, pi+pi-eta, " \
  "2pi 2pi0, 4pi eta, 6pi pi0. Modes grouped into Charged (C), K_S (K), Neutral (N) " \
  "categories with different chi2_4C cuts: <42 (C), <36 (K), <23 (N). " \
  "K_S reconstructed via secondary vertex fit pi+pi- with mass window |M-m_K_S|<12 MeV/c^2. " \
  "pi0/eta: 1C kinematic fit on gamma gamma, chi2_1C < 200. " \
  "Multiple candidates: smallest chi2 = chi2_4C + chi2_1C + chi2_PID + chi2_vertex retained.")
.note(:background_veto,
  "psi(3686) -> pi+ pi- J/psi veto: |RM(pi+pi-) - m_J/psi| > 3 sigma_J/psi. " \
  "psi(3686) -> eta J/psi veto: |M(gamma pi+pi-) - m_eta| > 3 sigma_eta.")
.note(:signal_extraction,
  "chi_c1 selected in mass window [3.48, 3.54] GeV/c^2 (RM(gamma)). " \
  "eta_c yield from combined maximum-likelihood fit to RM(gamma pi+pi-) over all 16 modes. " \
  "Background: Johnson function for RM(pi+pi-) dip + 5th-order polynomial. " \
  "90% C.L. Bayesian upper limit with systematic uncertainties incorporated.")
.note(:competing_hypothesis,
  "Competing 4pi background hypothesis fit (non-nominal, no chi2_cut) for ROOT-level veto: " \
  "kinematic_fit([pip, pip, pim, pim]) { constrain_four_momentum }.")

alg.with_decay_card(decay_card).apply(event_selection)
alg.execute_on([psip_data, psip_incMC, exMC_signal])