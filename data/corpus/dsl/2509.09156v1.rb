# BESIII DSL: psi(3686) -> gamma eta(1405) -> gamma f0(980) pi0 -> gamma pi+ pi- pi0
# Observation of eta(1405) -> f0(980) pi0 in psi(3686) radiative decays
# arXiv:2509.09156v1

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta(1405) PHSP;
    Enddecay

    Decay eta(1405)
    1.000 f_0 pi0 PHSP;
    Enddecay

    Decay f_0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_gamma_eta1405_f0_pi0"
  config.related_dataset = psip_data
  config.events = 5_000_000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

alg = Algorithm.new("PsipToGammaEta1405ToGammaF0Pi0")
alg.set_header(["PsipToGammaEta1405ToGammaF0Pi0Alg/PsipToGammaEta1405ToGammaF0Pi0.h"])
   .set_constant({ "ECMS" => [:double, 3.686] })

event_selection = Selection.new
event_selection.select_track {
                 cos_theta 0.93
                 Vz 10.0
                 Vxy 1.0
                 nChrg ">=2"
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :pion, against: [:electron, :kaon, :proton]
                 npip ">=1"
                 npim ">=1"
               }
               .select_photon {
                 tdc_emc_start 0
                 tdc_emc_end 14
                 energyThreshold_b 0.025
                 energyThreshold_e 0.050
                 nGam ">=3"
               }
               .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 20
               }

alg.note(:radiative_photon, "Radiative photon is the one not paired into pi0; energy window applied in ROOT analysis")
   .note(:pi0_mass_window, "|M(gamma gamma) - m_pi0| < 0.025 GeV/c^2 for pi0 candidate selection")
   .note(:eta_veto, "M(gamma gamma) not in [0.52, 0.57] GeV/c^2 to veto eta background")
   .note(:jpsi_veto, "|M(pi+ pi- gamma) - m_J/psi| > 0.05 GeV/c^2 to suppress J/psi background")
   .note(:omega_veto, "|M(pi0 gamma) - m_omega| > 0.04 GeV/c^2 to suppress omega background")
   .note(:photon_angle_to_track, "Minimum angle of 10 degrees between photon candidates and charged tracks")
   .note(:partial_wave_analysis, "PWA performed on eta(1405) -> f0(980)pi0 -> pi+pi-pi0 using Breit-Wigner parametrization")
   .note(:systematic_uncertainties, "Track/photon efficiency, PID efficiency, kinematic fit, mass window, PWA model, and branching fraction uncertainties evaluated")
   .with_decay_card(decay_card_signal)
   .apply(event_selection)

alg.execute_on([psip_data, psip_incMC, exMC_signal])