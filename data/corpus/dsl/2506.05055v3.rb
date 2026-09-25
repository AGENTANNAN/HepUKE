### Dataset preparation ###
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma pi0 pi0 pi0 PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_gamma_3pi0_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection ###
alg = Algorithm.new("JpsiGamma3Pi0")
alg.set_header(["JpsiGamma3Pi0Alg/JpsiGamma3Pi0.h"])
  .set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
event_selection.select_track { nTot "==0" }
  .select_photon {
    energyThreshold_b 0.025; energyThreshold_e 0.050
    tdc_emc_start 0; tdc_emc_end 14
    nGam ">=7"
  }
  # 1C kinematic fit to reconstruct pi0 from gamma pairs (chi2_1c < 10)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 10
    npi0 ">=3"
  }
  # 7C kinematic fit: J/psi -> gamma pi0 pi0 pi0 (4C E-p + 3 pi0 mass constraints)
  .kinematic_fit([:gamma, :pi0, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 40
  }
  # Competing hypothesis: J/psi -> gamma eta pi0 pi0 (for eta background veto in ROOT) — Rule T2
  .kinematic_fit([:gamma, :pi0, :pi0, :gamma, :gamma]) {
    constrain_four_momentum
  }

alg.with_decay_card(decay_card).apply(event_selection)
  .note(:background_veto, "events with |M(gamma pi0) - M(omega)| < 0.06 GeV/c^2 vetoed to suppress omega background; also veto |M(gamma_r gamma) - M(pi0)| < 0.02 GeV/c^2 for radiative photon mis-combination")
  .note(:pwa_method, "partial-wave analysis performed on M(pi0 pi0 pi0) < 1.6 GeV/c^2 using GPUPWA framework; intermediate resonances: eta(1405), f1(1285), f1(1420), f1(1510), plus non-resonant 0-+ and 1++ PHSP")
  .note(:f0_980_window, "events required to have at least one pi0 pi0 pair in f0(980) mass window [0.89, 1.09] GeV/c^2; narrow f0(980) width observed with detector resolution convolution (Gaussian sigma = 9.6 MeV)")
  .execute_on([jpsi_data, jpsi_incMC, exMC])