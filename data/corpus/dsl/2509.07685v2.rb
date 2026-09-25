# BESIII DSL: e+ e- -> e+ e- pi0 at psi(3770)
# pi0 transition form factor measurement via two-photon fusion
# Single-tag method: detect one scattered electron
# arXiv:2509.07685v2

psip3770_data = DatasetManager.real_data.find("712_3773")
psip3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 e+ e- pi0 CONEXC;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip3770_ee_pi0_TFF"
  config.related_dataset = psip3770_data
  config.events = 5_000_000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

alg = Algorithm.new("TwoPhotonFusionPi0TFF")
alg.set_header(["TwoPhotonFusionPi0TFFAlg/TwoPhotonFusionPi0TFF.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })

event_selection = Selection.new
event_selection.select_track {
                 cos_theta 0.93
                 Vz 10.0
                 Vxy 1.0
                 nChrg ">=1"
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :electron, against: [:pion, :kaon]
                 nep ">=1"
                 nem ">=1"
               }
               .select_photon {
                 tdc_emc_start 0
                 tdc_emc_end 14
                 energyThreshold_b 0.025
                 energyThreshold_e 0.050
                 nGam ">=2"
               }
               .kinematic_fit([:ep, :em, :gamma, :gamma]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 200
               }

alg.note(:two_photon_fusion, "e+e- -> e+e- pi0 via two-photon fusion process with single-tag method")
   .note(:electron_tag, "Scattered electron detected; tag energy/momentum used for Q^2 binning")
   .note(:energy_deposit_ratio, "E/p > 0.8 required for electron identification")
   .note(:photon_angle_to_track, "Minimum angle of 10 degrees between photon candidates and charged tracks")
   .note(:undetected_cos_theta, "cos_theta_undetected > 0.99 to suppress backgrounds")
   .note(:hadronic_cos_theta, "|cos_theta_H| < 0.8 in phi*f0 center-of-mass frame")
   .note(:R_gamma, "R_gamma <= 0.15 to suppress non-signal photons")
   .note(:extra_energy, "Sigma_E_extra <= 0.17 GeV to ensure event cleanliness")
   .note(:best_combination, "Best pi0 candidate chosen by minimum p_t* in gamma-gamma rest frame")
   .note(:tff_form_factor, "pi0 transition form factor F(Q^2) extracted from Q^2-differential cross section")
   .note(:q2_binning, "Q^2 bin boundaries: 0.0, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5 GeV^2/c^2")
   .with_decay_card(decay_card_signal)
   .apply(event_selection)

alg.execute_on([psip3770_data, psip3770_incMC, exMC_signal])