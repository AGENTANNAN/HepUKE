# DSL for BESIII paper 2408.03531v2: Measurement of psi(2S)->gamma pi0 and e+e- -> gamma pi0 form factor
# Final state: gamma gamma gamma (pi0 -> gamma gamma + radiative photon)
# Data: 2.7B psi(2S) events, 7.9 fb-1 psi(3773), 0.8 fb-1 off-resonance (3.650, 3.682 GeV)

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")
psi3773_data = DatasetManager.real_data.find("712_3773")
psi3773_incMC = DatasetManager.inclusive_mc.find("712_3773")
off_res_3650 = DatasetManager.real_data.find("708_3650")
off_res_3650_incMC = DatasetManager.inclusive_mc.find("708_3650")
off_res_3682 = DatasetManager.real_data.find("709_3682")
off_res_3682_incMC = DatasetManager.inclusive_mc.find("709_3682")

decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma pi0                    HELAMP 0 0 1 0 1 0 0 0 1 0 1 0 0 0;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                  PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "eeto_gamma_pi0"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection ###
# No charged tracks expected in signal; veto charged tracks
alg = Algorithm.new("Psi2S_GammaPi0")
alg.set_header(["Psi2S_GammaPi0Alg/Psi2S_GammaPi0.h"])
  .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz 10.0
                  Vr 1.0
                  nTot "==0"
                }
               .select_photon {
                  # Only barrel EMC photons
                  energyThreshold_b 0.025
                  # nGam "==3" applied via exact photon count
                }
               .for_each(:gamma) {
                  # Require exactly 3 photons; reject events with any endcap photons
                  define(:is_barrel) { cos_theta_gamma.abs < 0.8 }
                  where { is_barrel }
                }
               # 4C kinematic fit on 3 photons, no mass constraint on pi0
               .kinematic_fit([:gamma, :gamma, :gamma]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 40
                }

alg.note(:exactly_three_photons, "require exactly 3 photons; all in barrel EMC (|cos(theta)| < 0.8)")
   .note(:radiative_photon, "highest energy photon is the radiative photon; other two are pi0 daughters")
   .note(:helicity_angle_cut, "cos(theta_hel) < 0.7 to suppress e+e- -> gamma gamma gamma QED background")
   .note(:mdc_hits_suppression, "MDC hits between radial lines connecting IP and EMC showers < 8 (suppress gamma-conversion)")
   .note(:no_photon_isolation, "photon selection uses only barrel EMC; standard timing and shower-quality cuts applied")
   .note(:alternative_nominal_pi0, "M(gamma gamma) not constrained in BOSS; pi0 mass window and signal extraction done in ROOT")
   .with_decay_card(decay_card)
   .apply(event_selection)

root_files = alg.execute_on([psip_data, psip_incMC, exMC_signal, psi3773_data, psi3773_incMC,
                              off_res_3650, off_res_3650_incMC,
                              off_res_3682, off_res_3682_incMC])