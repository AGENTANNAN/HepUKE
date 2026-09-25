# eta' -> e+ e- omega  (via J/psi -> gamma eta')
# omega -> pi+ pi- pi0 ; pi0 -> gamma gamma
# Final state: 4 charged tracks (e+ e- pi+ pi-) + 3 photons (radiative + 2 from pi0)
# Data: ~10.087e9 J/psi events

### Datasets ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma  eta'                           HELAMP 1.0 0.0 1.0 0.0;
    Enddecay

    Decay eta'
    1.0000 e+  e-  omega                         PHSP;
    Enddecay

    Decay omega
    1.0000 pi+  pi-  pi0                         OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma  gamma                          PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "etap_to_eeomega_signal_mc"
  config.related_dataset = jpsi_data
  config.events          = 5000000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'exMC_etap_eeomega_config')

alg = Algorithm.new("EtapEEomega")
alg.set_header(["EtapEEomegaAlg/EtapEEomega.h"])
   .set_constant({"ECMS" => [:double, 3.097]})

sel = Selection.new
sel.select_track {
       cos_theta 0.93          # |cos(theta)| <= 0.93
       Vr        2.0           # |Vxy| < 2 cm
       Vz        20.0          # |Vz| < 20 cm
       nChrp     "==2"         # e+, pi+
       nChrn     "==2"         # e-, pi-
       nNet      "==0"
       nTot      "==4"
     }
    .select_photon {
       tdc_emc_start     0
       tdc_emc_end       14    # 0 - 700 ns
       angle_to_track    10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam              ">=3" # radiative gamma + 2 photons from pi0
     }
    .pid(method: :probability) {
       prob_cut 0.0
       identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.6
       identify :pion, against: [:kaon, :electron]
       npip ">=1"
       npim ">=1"
       nlp  ">=1"   # e+
       nlm  ">=1"   # e-
     }
    .kalman_kinematic_fit([:gamma, :gamma]) {
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
       chi2_cut 200
       npi0 ">=1"
     }
    # Main 4C kinematic fit: J/psi -> gamma e+ e- pi+ pi- pi0
    .kinematic_fit([:gamma, :lp, :lm, :pip, :pim, :pi0]) {
       nominal
       constrain_four_momentum
       chi2_cut 70    # loose combined chi^2_4C + PID; ROOT-side tighter cut
     }
    # Background hypothesis test: J/psi -> pi+ pi- pi+ pi- gamma gamma gamma
    # (e/pi mis-ID). Store chi2 for later comparison; competing-hypothesis fit.
    .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim, :pip, :pim]) {
       constrain_four_momentum
     }

alg
  .note(:radiative_gamma_energy_cut,
        "The most energetic photon in the event is required to have E > 1 GeV (J/psi -> gamma eta' radiative photon at 1.4 GeV).")
  .note(:combined_chisq_selection,
        "Best combination selected by minimum chi^2_(4C+PID) = chi^2_4C + sum_i chi^2_PID(i); require chi^2 < 70.")
  .note(:e_pi_mis_ID_veto,
        "Reject events whose 4C chi^2 under J/psi -> pi+pi-pi+pi- gamma gamma gamma is smaller than under the signal hypothesis (to suppress e/pi mis-ID).")
  .note(:pi0_mass_window,
        "|M(gamma gamma) - M(pi0)_PDG| < 0.015 GeV/c^2 for the two lowest-energy photons (pi0 candidate).")
  .note(:photon_conversion_veto,
        "Photon-conversion veto using conversion vertex Rxy: retain events with Rxy < 2 cm; for events with Rxy >= 2 cm, reject those satisfying both cos(theta_eg) > 0 AND -1.0 < Delta_xy < 0.8 cm.")
  .note(:pi0_from_gamma_ee_veto,
        "To suppress pi0 -> gamma e+ e-, require |M(gamma1 gamma2) - M(pi0)| < |M(gamma1 e+ e-) - M(pi0)| AND |M(gamma1 gamma2) - M(pi0)| < |M(gamma2 e+ e-) - M(pi0)|.")
  .note(:eta_gamma_ee_veto,
        "Veto 0.52 < M(gamma3 e+ e-) < 0.58 GeV/c^2 to suppress eta -> gamma e+ e- background; also require M(gamma1 e+ e-) < 0.5 GeV/c^2.")
  .note(:omega_mass_window,
        "M(pi+ pi- pi0 e+ e-) [i.e. eta' candidate] required to be in (0.90, 1.02) GeV/c^2 to suppress eta' -> eta pi+pi-, eta -> pi+pi-pi0.")
  .note(:etap_signal_region_for_TFF,
        "For TFF fit only: further refine to M(pi+pi-pi0 e+e-) in [0.94, 0.98] GeV/c^2 and M(pi+pi-pi0 e+e-)-M(pi+pi-pi0) in [0.14, 0.20] GeV/c^2.")
  .with_decay_card(decay_card_signal)
  .apply(sel)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
