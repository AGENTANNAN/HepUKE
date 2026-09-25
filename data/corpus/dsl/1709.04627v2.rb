# -*- coding: utf-8 -*-
# BESIII BOSS DSL — arXiv:1709.04627v2
# Study of the matrix element for the decays eta' -> eta pi+ pi- and eta' -> eta pi0 pi0
# Data: 1.31e9 J/psi events (BOSS sample 708_3097)

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi data sample
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # J/psi inclusive MC (background study)

# Decay card — charged mode: J/psi -> gamma eta', eta' -> eta pi+ pi-, eta -> gamma gamma
decay_card_charged = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'  PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi-  PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card — neutral mode: J/psi -> gamma eta', eta' -> eta pi0 pi0,
# eta -> gamma gamma and pi0 -> gamma gamma (two identical pi0 daughters share one Decay block)
decay_card_neutral = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'  PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi0 pi0  PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples for both signal modes
exMC_charged = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "etap_eta_pip_pim_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_charged
  config.cross_section   = :default
end

exMC_neutral = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "etap_eta_pi0_pi0_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_neutral
  config.cross_section   = :default
end

### Event selection — Mode I: eta' -> eta pi+ pi- (charged decay mode) ###
alg_name_charged = "EtapToEtaPiPiCharged"
alg_charged = Algorithm.new(alg_name_charged)
alg_charged.set_header(["#{alg_name_charged}Alg/#{alg_name_charged}.h"])
           .set_constant({"ECMS" => [:double, 3.097]})

charged_selection = Selection.new
  .select_track {                    # two oppositely charged tracks
      cos_theta 0.93                 # |cos(theta)| < 0.93
      Vz        10.0                 # |Vz| < 10 cm along the beam direction
      Vr        1.0                  # Vr < 1 cm in the plane perpendicular to the beam
      nChrp     "==1"                # one positive track
      nChrn     "==1"                # one negative track
      nNet      "==0"                # opposite charge
  }
  .select_photon {                   # at least three photons (1 radiative + 2 from eta)
      tdc_emc_start     0            # EMC cluster timing window: 0 <= T <= 700 ns
      tdc_emc_end       14
      energyThreshold_b 0.025        # E(barrel) > 25 MeV
      energyThreshold_e 0.050        # E(endcap) > 50 MeV
      angle_to_track    10.0         # photon-track opening angle > 10 degrees
      nGam              ">=3"        # at least three photons
  }
  .pid(method: :probability) {       # pion identification for the two charged tracks
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]
      npip "==1"
      npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit: eta -> gamma gamma
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 25
      neta ">=1"
  }
  .kinematic_fit([:gamma, :eta, :pip, :pim]) {  # 6C: 4C + m(eta) + m(eta')
      nominal
      constrain_four_momentum
      invariant_mass_of(:eta, :pip, :pim).constrain_to_nominal_mass_of(:etap)
      chi2_cut 200                    # loose BOSS-level cut; published chi2_6C < 100 applied in ROOT
  }

alg_charged
  .note(:radiation_photon_assignment, "the photon with the maximum energy in the event is
    assumed to be the radiative photon from J/psi -> gamma eta'; the remaining photons are
    used for eta -> gamma gamma. Implemented via the combination ordering of the kinematic fit.")
  .note(:efficiency_curve, "momentum-dependent tracking efficiency correction for charged
    pions (from J/psi -> p pbar pi+ pi- control sample) and momentum-dependent eta
    reconstruction efficiency correction (from J/psi -> gamma eta pi+ pi- control sample)
    are applied; the associated uncertainties are estimated by alternative fits.")
  .note(:kinematic_fit_variation, "systematic check performed with a 4C kinematic fit
    (energy-momentum conservation only) instead of the nominal 6C fit.")
  .with_decay_card(decay_card_charged)
  .apply(charged_selection)

### Event selection — Mode II: eta' -> eta pi0 pi0 (neutral decay mode) ###
alg_name_neutral = "EtapToEtaPi0Pi0Neutral"
alg_neutral = Algorithm.new(alg_name_neutral)
alg_neutral.set_header(["#{alg_name_neutral}Alg/#{alg_name_neutral}.h"])
           .set_constant({"ECMS" => [:double, 3.097]})

neutral_selection = Selection.new
  .select_track {                    # no charged track
      nTot "==0"
  }
  .select_photon {                   # at least seven photons (1 radiative + 2 (eta) + 4 (2 pi0))
      tdc_emc_start     0            # EMC cluster timing (window referenced to the most
      tdc_emc_end       14           # energetic photon, -500 <= T <= 500 ns)
      energyThreshold_b 0.025        # E(barrel) > 25 MeV
      energyThreshold_e 0.050        # E(endcap) > 50 MeV
      nGam              ">=7"        # at least seven photons
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {    # 1C fit: pi0 -> gamma gamma (two pi0 candidates)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 "==2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {    # 1C fit: eta -> gamma gamma
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 25
      neta ">=1"
  }
  .kinematic_fit([:gamma, :eta, :pi0, :pi0]) { # 8C: 4C + 3 photon-pair masses + m(eta')
      nominal
      constrain_four_momentum
      invariant_mass_of(:eta, :pi0, :pi0).constrain_to_nominal_mass_of(:etap)
      chi2_cut 200                   # loose BOSS-level cut; published chi2_8C < 100 applied in ROOT
  }

alg_neutral
  .note(:radiation_photon_assignment, "the photon with the largest energy in the event is
    assumed to be the radiative photon originating from J/psi -> gamma eta'; the remaining
    clusters are paired into pi0/eta -> gamma gamma candidates.")
  .note(:pi0_decay_angle, "pi0 miscombination suppressed by requiring |cos(theta_decay)| < 0.95,
    where theta_decay is the polar angle of one decay photon in the gamma gamma rest frame
    with respect to the pi0 flight direction.")
  .note(:photon_miscombination, "fraction of events with wrongly matched photon pairs is 2.7%;
    systematic uncertainty from miscombination evaluated with truth-tagged MC events.")
  .note(:efficiency_curve, "momentum-dependent pi0 reconstruction efficiency correction
    (derived from the J/psi -> pi+ pi- pi0 control sample) is applied, in addition to the
    eta and photon reconstruction efficiency corrections.")
  .note(:kinematic_fit_variation, "systematic check performed with a 6C kinematic fit
    (energy-momentum conservation plus the eta' mass constraint) instead of the nominal 8C fit.")
  .note(:background_veto, "peaking background eta' -> pi0 pi0 pi0 and the flat contribution from
    J/psi -> omega eta (omega -> gamma pi0, eta -> pi0 pi0 pi0) contribute about 0.9% and are
    neglected in the nominal Dalitz-plot-parameter determination.")
  .with_decay_card(decay_card_neutral)
  .apply(neutral_selection)

### Execution ###
root_files_charged = alg_charged.execute_on([jpsi_data, jpsi_incMC, exMC_charged])
root_files_neutral = alg_neutral.execute_on([jpsi_data, jpsi_incMC, exMC_neutral])
