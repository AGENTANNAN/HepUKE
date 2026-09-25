# 1902.04862v1: J/psi -> omega eta' pi+ pi- and search for X(1835) -> eta' pi+ pi-
# Data: J/psi, 1.31 x 10^9 events
# omega -> pi+ pi- pi0 (pi0->gamma gamma), eta' -> eta pi+ pi- (eta->gamma gamma)
# Final state: 6 charged pions + 4 photons, 4C kinematic fit

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_jpsi = DatasetManager.real_data.find("708_3097")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")

# Signal: J/psi -> omega eta' pi+ pi-
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.000 omega eta' pi+ pi- PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  Decay eta'
  1.000 eta pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_OmegaEtapPiPi"
  config.related_dataset = data_jpsi
  config.events = 100000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

# Signal MC for X(1835): J/psi -> omega X(1835), X(1835) -> eta' pi+ pi-
decay_card_X1835 = <<~DECAYCARD
  Decay J/psi
  1.000 omega X_1835 PHSP;
  Enddecay

  Decay X_1835
  1.000 eta' pi+ pi- PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  Decay eta'
  1.000 eta pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_X1835 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_OmegaX1835_EtapPiPi"
  config.related_dataset = data_jpsi
  config.events = 100000
  config.decay_card = decay_card_X1835
  config.cross_section = :default
end

alg = Algorithm.new("OmegaEtapPiPi")
alg.set_header(["OmegaEtapPiPiAlg/OmegaEtapPiPi.h"])
    .set_constant({"ECMS" => [:double, 3.097]})

sel = Selection.new
sel.select_track {
      cos_theta 0.93
      Vz 100.0
      Vr 10.0
      nChrp "==3"
      nChrn "==3"
      nNet "==0"
    }
    .select_photon {
      tdc_emc_start 0
      tdc_emc_end 14
      angle_to_track 5.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam ">=4"
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon]
      npip "==3"
      npim "==3"
    }
    .remove([:pip <= :chrgp, :pim <= :chrgn])
    .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 25
      neta ">=1"
    }
    .kinematic_fit([:pip, :pim, :pip, :pim, :pip, :pim, :pi0, :eta]) {
      nominal
      constrain_four_momentum
      chi2_cut 60
    }
alg.with_decay_card(decay_card_signal).apply(sel)

alg.note(:omega_etap_pairing,
  "Best omega and eta' candidate chosen by minimizing " \
  "sqrt((M(pi+pi-pi0)-m_omega)^2 + (M(eta pi+pi-)-m_etap )^2); " \
  "omega mass window |M(pi+pi-pi0)-m_omega| < 22 MeV/c^2; " \
  "eta' mass window |M(eta pi+pi-)-m_etap | < 12 MeV/c^2; " \
  "all applied in ROOT analysis")
alg.note(:X1835_upper_limit,
  "X(1835) search: background-subtracted eta' pi+ pi- mass spectrum fitted with " \
  "BW (M=1836.5, Gamma=190 MeV) or Flatte function; " \
  "upper limit B(J/psi->omega X(1835), X(1835)->eta' pi+ pi-) < 6.2 x 10^-5 at 90% CL")
alg.note(:two_dimensional_fit,
  "2D fit to pi+pi-pi0 vs eta pi+pi- performed to extract signal yield; " \
  "omega parameterized by BW convoluted with double Gaussian, " \
  "eta' by double Gaussian; backgrounds by third-order polynomials")

alg.execute_on([data_jpsi, incMC_jpsi, exMC_signal, exMC_X1835])