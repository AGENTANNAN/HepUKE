# First study of J/psi -> gamma eta pi0 (eta -> gamma gamma, pi0 -> gamma gamma)
# with (223.7 +/- 1.4) x 10^6 J/psi events at BESIII.
# Final state topology: five photons, no charged tracks.

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi dataset at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Matching inclusive MC sample

# Signal: J/psi -> gamma eta pi0 -> 5 gamma (phase space)
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta pi0      PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma        PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma        PHSP;
  Enddecay

  End
DECAYCARD

# Intermediate-resonance signal samples: J/psi -> gamma a0(980), a0(980) -> eta pi0
decay_card_a0 = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma a_00         PHSP;
  Enddecay

  Decay a_00
  1.0000 eta pi0            PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma        PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma        PHSP;
  Enddecay

  End
DECAYCARD

# Intermediate-resonance signal samples: J/psi -> gamma a2(1320), a2(1320) -> eta pi0
decay_card_a2 = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma a_20         PHSP;
  Enddecay

  Decay a_20
  1.0000 eta pi0            PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma        PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma        PHSP;
  Enddecay

  End
DECAYCARD

# Background: J/psi -> eta omega (eta -> gamma gamma, omega -> gamma pi0)
decay_card_etaomega = <<~DECAYCARD
  Decay J/psi
  1.0000 eta omega          PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma        PHSP;
  Enddecay

  Decay omega
  1.0000 gamma pi0          PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma        PHSP;
  Enddecay

  End
DECAYCARD

# Background: J/psi -> eta phi (eta -> gamma gamma, phi -> gamma pi0)
decay_card_etaphi = <<~DECAYCARD
  Decay J/psi
  1.0000 eta phi            PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma        PHSP;
  Enddecay

  Decay phi
  1.0000 gamma pi0          PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma        PHSP;
  Enddecay

  End
DECAYCARD

# Background: J/psi -> gamma eta', eta' -> 2 pi0 eta  (via eta' -> pi0 pi0 eta)
decay_card_etap_2pi0eta = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'         PHSP;
  Enddecay

  Decay eta'
  1.0000 pi0 pi0 eta        PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma        PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma        PHSP;
  Enddecay

  End
DECAYCARD

# Background: J/psi -> gamma eta', eta' -> gamma omega  (with omega -> gamma pi0)
decay_card_etap_gammaomega = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'         PHSP;
  Enddecay

  Decay eta'
  1.0000 gamma omega        PHSP;
  Enddecay

  Decay omega
  1.0000 gamma pi0          PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma        PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_eta_pi0_5gamma"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

exMC_a0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_a0_980_eta_pi0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_a0
  config.cross_section   = :default
end

exMC_a2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_a2_1320_eta_pi0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_a2
  config.cross_section   = :default
end

exMC_etaomega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_eta_omega"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_etaomega
  config.cross_section   = :default
end

exMC_etaphi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_eta_phi"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_etaphi
  config.cross_section   = :default
end

exMC_etap_2pi0eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etap_2pi0eta"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_etap_2pi0eta
  config.cross_section   = :default
end

exMC_etap_gammaomega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etap_gammaomega"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_etap_gammaomega
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Single signal channel: five photons, no charged tracks, 4C kinematic fit under
# the e+e- -> 5 gamma hypothesis.  Only one Algorithm is needed.
alg_name = "JpsiGammaEtaPi0"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
event_selection
  .select_track {
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
     nChrp     "==0"
     nChrn     "==0"
  }
  .select_photon {
     tdc_emc_start     0
     tdc_emc_end       14
     angle_to_track    10.0
     energyThreshold_b 0.025
     energyThreshold_e 0.050
     nGam              ">=5"
  }
  # 4C kinematic fit: energy-momentum conservation under e+e- -> 5 gamma.
  # The fit iterates over all photon combinations and keeps the one with the
  # smallest chi2_4C.
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma]) {
     nominal
     constrain_four_momentum
     chi2_cut 200
  }

algorithm
  .note(:track_definition,
        "signal candidates require no reconstructed charged particles in the event.")
  .note(:emc_timing,
        "EMC timing information used to suppress electronic noise and energy " \
        "depositions unrelated to the event; photon candidates required to be " \
        "within 50 ns relative to the most energetic shower. " \
        "(select_photon tdc_emc_start/tdc_emc_end cover the EMC timing window.)")
  .note(:photon_energy_thresholds,
        "25 MeV minimum deposited energy in the EMC barrel (|cos theta| < 0.8) and " \
        "50 MeV in the end caps (0.86 < |cos theta| < 0.92), where theta is the " \
        "polar angle of the electromagnetic shower.")
  .note(:chi2_cut_published,
        "Published selection requires chi2_4C < 30; kept loose (200) in BOSS and " \
        "applied at the optimal value in ROOT.")
  .note(:eta_pi0_assignment,
        "among the five photons, the eta and pi0 photon pairs are chosen by " \
        "minimising Delta = sqrt((M(g1 g2) - m_eta)^2 + (M(g3 g4) - m_pi0)^2) " \
        "over all photon combinations; performed in ROOT on the four-momenta " \
        "updated by the 4C fit.")
  .note(:eta_pi0_mass_regions,
        "eta signal region |M(gamma gamma) - m_eta| < 0.024 GeV/c^2; pi0 signal " \
        "region |M(gamma gamma) - m_pi0| < 0.015 GeV/c^2; pi0 sidebands " \
        "0.030 < |M(gamma gamma) - m_pi0| < 0.045 GeV/c^2.")
  .note(:background_veto,
        "background with two pi0 in the final state (e.g. J/psi -> gamma pi0 pi0) " \
        "suppressed by rejecting events for which any photon combination satisfies " \
        "Delta_pi0 = sqrt((M(g1 g2) - m_pi0)^2 + (M(g3 g4) - m_pi0)^2) < 0.05 GeV/c^2.")
  .note(:background_veto_omega,
        "omega background from J/psi -> omega eta (omega -> gamma pi0) and " \
        "J/psi -> phi eta (phi -> gamma pi0) rejected by requiring " \
        "|M(gamma pi0) - m_omega| > 0.07 GeV/c^2.")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

### Execute the algorithm on data, inclusive MC and all exclusive MC samples ###
root_files = algorithm.execute_on([
  jpsi_data,
  jpsi_incMC,
  exMC_signal,
  exMC_a0,
  exMC_a2,
  exMC_etaomega,
  exMC_etaphi,
  exMC_etap_2pi0eta,
  exMC_etap_gammaomega,
])
