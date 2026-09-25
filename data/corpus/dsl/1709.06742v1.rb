# -*- coding: utf-8 -*-
# BESIII BOSS DSL — arXiv:1709.06742v1
# Improved measurements of two-photon widths of the chi_cJ states and helicity analysis
# for chi_c2 -> gamma gamma
# Data: 448.1e6 psi(3686) events (BOSS sample 709_3686)

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) data sample
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # psi(3686) inclusive MC

# Decay cards for the three signal hypotheses psi(3686) -> gamma chi_cJ, chi_cJ -> gamma gamma.
# chi_c1 -> gamma gamma is forbidden by the Landau-Yang theorem; it is generated
# hypothetically to extract the upper limit.
decay_card_chi_c0 = <<~DECAYCARD
    Decay psi(3686)
    1.0000 gamma chi_c0  PHSP;
    Enddecay

    Decay chi_c0
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chi_c2 = <<~DECAYCARD
    Decay psi(3686)
    1.0000 gamma chi_c2  PHSP;
    Enddecay

    Decay chi_c2
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chi_c1 = <<~DECAYCARD
    Decay psi(3686)
    1.0000 gamma chi_c1  PHSP;
    Enddecay

    Decay chi_c1
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples (1.2 million events per signal channel as in the paper)
exMC_chi_c0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3686_to_gamma_chi_c0_gg_exclusive_mc"
  config.related_dataset = psip_data
  config.events          = 1_200_000
  config.decay_card      = decay_card_chi_c0
  config.cross_section   = :default
end

exMC_chi_c2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3686_to_gamma_chi_c2_gg_exclusive_mc"
  config.related_dataset = psip_data
  config.events          = 1_200_000
  config.decay_card      = decay_card_chi_c2
  config.cross_section   = :default
end

exMC_chi_c1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3686_to_gamma_chi_c1_gg_exclusive_mc"
  config.related_dataset = psip_data
  config.events          = 1_200_000
  config.decay_card      = decay_card_chi_c1
  config.cross_section   = :default
end

### Event selection — common to chi_c0, chi_c1 and chi_c2 (identical final state) ###
# Selection chain: no charged track, exactly three photons, 4C kinematic fit
event_selection_common = Selection.new
  .select_track {                  # no charged track is required
      nTot "==0"
  }
  .select_photon {                 # three photon candidates
      energyThreshold_b 0.070      # E(gamma) > 70 MeV (barrel)
      energyThreshold_e 0.070      # E(gamma) > 70 MeV (endcap)
      nGam              "==3"      # exactly three photon candidates
  }
  .kinematic_fit([:gamma, :gamma, :gamma]) {   # 4C fit to the initial e+e- four-momentum
      nominal
      constrain_four_momentum
      chi2_cut 200                 # loose BOSS-level cut; published chi2_4C <= 80 applied in ROOT
  }

### Algorithm per signal channel ###
alg_name_chi_c0 = "Psi3686ToGammaChiC0ToGG"
alg_chi_c0 = Algorithm.new(alg_name_chi_c0)
alg_chi_c0.set_header(["#{alg_name_chi_c0}Alg/#{alg_name_chi_c0}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

alg_name_chi_c2 = "Psi3686ToGammaChiC2ToGG"
alg_chi_c2 = Algorithm.new(alg_name_chi_c2)
alg_chi_c2.set_header(["#{alg_name_chi_c2}Alg/#{alg_name_chi_c2}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

alg_name_chi_c1 = "Psi3686ToGammaChiC1ToGG"
alg_chi_c1 = Algorithm.new(alg_name_chi_c1)
alg_chi_c1.set_header(["#{alg_name_chi_c1}Alg/#{alg_name_chi_c1}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

# Common BOSS-side notes for all three channels
[alg_chi_c0, alg_chi_c2, alg_chi_c1].each do |alg|
  alg.note(:photon_angle_cut, "each photon candidate is required to have |cos(theta)| < 0.75,
    where theta is the angle of the photon with respect to the positron beam direction; this
    suppresses the continuum background e+e- -> gamma gamma(gamma) whose two energetic photons
    populate the forward and backward regions.")
  alg.note(:radiation_photon_assignment, "the photon with the smallest energy among the three
    candidates, E(gamma_1), is taken as the radiative photon; its spectrum is used to extract
    the chi_c0/chi_c2/chi_c1 signal yields.")
  alg.note(:background_shape, "the non-peaking continuum background shape
    f_bg = p0 + p1 E + p2 E^2 + p3 E^a is validated with the psi(3770) data sample
    (BOSS 712_3773, sqrt(s) = 3.773 GeV) and the off-resonance data sample taken at
    sqrt(s) = 3.65 GeV (BOSS 709_3650, 48 pb^-1).")
end

# chi_c2 helicity amplitude analysis notes
alg_chi_c2
  .note(:helicity_generation, "the signal MC for psi(3686) -> gamma chi_c2, chi_c2 -> gamma gamma
    is generated in a pure helicity-two process, because the helicity-zero component is
    negligible (f_0/2 = (0.0 +/- 0.6 +/- 1.2) x 10^-2); a 2% helicity-zero fraction is included
    to assign the 0.2% systematic uncertainty.")
  .note(:efficiency_curve, "the unbinned maximum-likelihood fit of Eq. (4) to the angular
    distribution uses the acceptance-corrected averages a_n_bar of the 12 angular factors
    computed from a phase-space MC sample after all selection criteria.")
  .note(:background_veto, "backgrounds are subtracted using sideband regions
    0.07 < E(gamma_1) < 0.09 GeV (lower) and 0.16 < E(gamma_1) < 0.19 GeV (upper); the
    chi_c0 -> gamma gamma contamination in the chi_c2 signal region (0.11 < E(gamma_1) < 0.14 GeV)
    is 0.044% and is neglected.")

alg_chi_c0.with_decay_card(decay_card_chi_c0).apply(event_selection_common.dup)
alg_chi_c2.with_decay_card(decay_card_chi_c2).apply(event_selection_common.dup)
alg_chi_c1.with_decay_card(decay_card_chi_c1).apply(event_selection_common.dup)

### Execution ###
root_files_chi_c0 = alg_chi_c0.execute_on([psip_data, psip_incMC, exMC_chi_c0])
root_files_chi_c2 = alg_chi_c2.execute_on([psip_data, psip_incMC, exMC_chi_c2])
root_files_chi_c1 = alg_chi_c1.execute_on([psip_data, psip_incMC, exMC_chi_c1])
