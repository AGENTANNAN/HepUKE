# -*- coding: utf-8 -*-
# BESIII BOSS DSL — arXiv:1709.10236v4
# Measurement of e+e- -> Lambda Lambda at sqrt(s) = 2.2324 GeV (threshold) and at
# sqrt(s) = 2.400, 2.800, 3.080 GeV + determination of the Lambda effective form factor.
#
# At threshold the nucleons/antinucleons are too soft to be detected, so the final state
# cannot be fully reconstructed: mode I reconstructs the two low-momentum pions and mode II
# reconstructs anti-n + pi0 from anti-Lambda -> anti-n pi0, in both cases the (anti-)Lambda
# recoil is left untagged (partial reconstruction, no kinematic fit).
# At 2.400 / 2.800 / 3.080 GeV the full final state p pbar pi+ pi- is reconstructed from
# Lambda -> p pi- and anti-Lambda -> anti-p pi+ with secondary vertex fits.

### Dataset description ###
# Threshold point sqrt(s) = 2.2324 GeV (1.0 MeV above the Lambda Lambda mass threshold)
data_2232  = DatasetManager.real_data.find("713_Rscan_2232")
incMC_2232 = DatasetManager.inclusive_mc.find("713_Rscan_2232")

# Higher c.m. energy points
data_2400  = DatasetManager.real_data.find("713_Rscan_2396")
incMC_2400 = DatasetManager.inclusive_mc.find("713_Rscan_2396")
data_2800  = DatasetManager.real_data.find("713_Rscan_2800")
incMC_2800 = DatasetManager.inclusive_mc.find("713_Rscan_2800")
data_3080  = DatasetManager.real_data.find("713_Rscan_3080")
incMC_3080 = DatasetManager.inclusive_mc.find("713_Rscan_3080")

# ---------------------------------------------------------------------------
# Decay cards
# ---------------------------------------------------------------------------
# Threshold, mode I: e+e- -> Lambda anti-Lambda, Lambda -> p+ pi-, anti-Lambda -> anti-p- pi+.
# Only the two pions are reconstructed; the (anti-)proton and both (anti-)Lambdas are not
# detected. At threshold the process is generated flat in phase space (G_E = G_M).
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.0 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0 p+ pi- PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0 anti-p- pi+ PHSP;
    Enddecay

    End
DECAYCARD

# Threshold, mode II: anti-Lambda -> anti-n0 pi0 is reconstructed (n = anti-neutron,
# pi0 -> gamma gamma); the other Lambda decays inclusively and is not reconstructed.
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.0 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0 p+ pi- PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0 anti-n0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Higher c.m. energies (2.400 / 2.800 / 3.080 GeV): the signal is generated with ConExc,
# which includes the higher-order processes with one radiative photon and the vacuum
# polarization; ConExc mode 2 = Lambda anti-Lambda (2.23 - 5.00 GeV). 'Particle vpho' is
# omitted so that each scan point gets its own sqrt(s) injected automatically.
decay_card_conexc = <<~DECAYCARD
    Decay vpho
    1 ConExc 2;
    Enddecay

    Decay vhdr
    1 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0 p+ pi- PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0 anti-p- pi+ PHSP;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive MC samples
# ---------------------------------------------------------------------------
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "ee_to_LambdaLambdabar_threshold_modeI_exclusive_mc"
  config.related_dataset = data_2232
  config.events          = 200_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "ee_to_LambdaLambdabar_threshold_modeII_exclusive_mc"
  config.related_dataset = data_2232
  config.events          = 200_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# One ConExc signal sample per higher-energy scan point (same card at every point)
exMCs_highE = DatasetManager.create_exclusive_mc_for([data_2400, data_2800, data_3080]) do |config|
  config.sample_name   = "ee_to_LambdaLambdabar_conexc_exclusive_mc"
  config.events        = 200_000
  config.decay_card    = decay_card_conexc
  config.cross_section = :default
end

# ---------------------------------------------------------------------------
# Event selection — mode I at sqrt(s) = 2.2324 GeV
# ---------------------------------------------------------------------------
alg_name_modeI = "LambdaLambdaThresholdModeI"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 2.2324]})

selection_modeI = Selection.new
  .select_track {                    # exactly two good charged tracks of opposite charge
      cos_theta 0.93                 # |cos(theta)| < 0.93
      Vz        10.0                 # |Vz| < 10 cm
      Vr        1.0                  # Vxy < 1 cm
      nChrp     "==1"
      nChrn     "==1"
      nNet      "==0"
  }
  .pid(method: :probability) {       # both tracks must be pions (dE/dx + TOF)
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]
      npip "==1"
      npim "==1"
  }
  .partial_rec([3, 6]) do            # recID 3 = pi- (from Lambda), recID 6 = pi+ (from anti-Lambda)
    # The proton, the antiproton and both Lambda's are left untagged; the recoil of the
    # two-pion system is the untagged (anti-)Lambda pair. The antiproton is not detected
    # directly but identified through its annihilation with the detector material.
  end

alg_modeI
  .note(:background_veto, "the antiproton annihilates with the detector material (mostly the beam
    pipe) and is identified via the distribution of Vr, the largest Vxy among all tracks in the
    event apart from the two good pions; the signal peak at Vr ~ 3 cm (the IP-to-beam-pipe
    distance) is fitted with the MC signal shape and a background shape taken from the pion
    momentum sideband [0.15, 0.18] GeV/c.")
  .note(:efficiency_curve, "the charged pion momentum is required to be within [0.08, 0.11] GeV/c,
    as expected from Lambda (anti-Lambda) decay; the not-expressible momentum-window cut is
    applied at the ROOT level.")
  .note(:tracking_efficiency, "low-momentum tracking efficiency (12.3%) and PID efficiency (1.0%)
    differences between data and MC are evaluated with the control sample
    J/psi -> p pbar pi+ pi-.")
  .note(:partial_reconstruction, "mode I does not fully reconstruct the final state; the number of
    signal events is extracted from a fit to the Vr distribution, not from a kinematic fit.")
  .note(:radiative_correction, "the radiative correction factor 1 + delta accounts for the energy
    spread of the beams (0.48 MeV, scaled from the J/psi peak) and for ISR photon emission, and
    for vacuum polarization.")
  .with_decay_card(decay_card_modeI)
  .apply(selection_modeI)

# ---------------------------------------------------------------------------
# Event selection — mode II at sqrt(s) = 2.2324 GeV
# ---------------------------------------------------------------------------
alg_name_modeII = "LambdaLambdaThresholdModeII"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 2.2324]})

selection_modeII = Selection.new
  .select_track {                    # at most one good charged track
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nTot      "<=1"
  }
  .select_photon {                   # at least three neutral candidates
      energyThreshold_b 0.025        # E > 25 MeV in the barrel
      energyThreshold_e 0.050        # E > 50 MeV in the endcap
      angle_to_track    10.0         # no charged track within 10 degrees
      nGam              ">=3"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit: pi0 -> gamma gamma
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 20                    # only events with chi2_1C < 20 are accepted
      npi0 ">=1"
  }
  .partial_rec([5, 6]) do            # recID 5 = anti-n0 (most energetic shower), recID 6 = pi0
    # The Lambda recoiling against the anti-n pi0 system is left untagged.
  end

alg_modeII
  .note(:background_veto, "the most energetic shower is assumed to be the anti-neutron (it
    annihilates in the EMC producing secondary particles with a total energy up to 2 m_n), the
    remaining showers are combined into pi0 candidates; the angle between the pi0 and the
    anti-neutron momentum directions is required to be larger than 140 degrees and the energy
    asymmetry |E(gamma_1) - E(gamma_2)| / p(pi0) must be less than 0.95.")
  .note(:bdt_selection, "a boosted decision tree (BDT) built on eight EMC-related variables
    (energy deposition within a 40 degree cone, deposited energy, energy-seed deposit, number of
    hits within a 40 degree cone, number of hits, lateral moment, second moment, deposition
    shape) is used to separate the signal from the hadronic and beam-associated backgrounds; an
    optimal requirement on the BDT output is applied. The signal yield is extracted from an
    unbinned maximum-likelihood fit to p(pi0) (signal: MC shape convoluted with a Gaussian,
    background: linear).")
  .note(:background_veto, "the dominant backgrounds are inclusive hadronic final states with
    multiple pi0's and beam-associated events, studied with a dedicated non-colliding-beam data
    sample; the anti-neutron and pi0 selection efficiencies and the BDT requirement carry
    2.2%/2.3%/4.8% systematic uncertainties.")
  .note(:efficiency_curve, "the signal MC is generated with a PHSP generator; the efficiency is
    corrected for the ISR / energy-spread effects through the radiative correction factor.")
  .with_decay_card(decay_card_modeII)
  .apply(selection_modeII)

# ---------------------------------------------------------------------------
# Event selection — sqrt(s) = 2.400, 2.800 and 3.080 GeV (full reconstruction)
# ---------------------------------------------------------------------------
alg_name_highE = "LambdaLambdaFullReco"
alg_highE = Algorithm.new(alg_name_highE)
alg_highE.set_header(["#{alg_name_highE}Alg/#{alg_name_highE}.h"])
# ECMS is not set: this algorithm runs at three different c.m. energies (cross-section scan).

selection_highE = Selection.new
  .select_track {                    # four good charged tracks
      cos_theta 0.93
      Vz        30.0                 # |Vz| < 30 cm
      Vr        10.0                 # Vxy < 10 cm
      nChrp     "==2"
      nChrn     "==2"
      nNet      "==0"
  }
  .pid(method: :probability) {       # one proton-antiproton pair and one pion pair
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      identify :pion, against: [:kaon, :proton]
      nprp "==1"
      nprm "==1"
      npip "==1"
      npim "==1"
  }
  .secondary_vertex_fit([:prp, :pim]) {        # Lambda -> p pi-
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {        # anti-Lambda -> anti-p pi+
      build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Lambda, :Lambda_bar]) {     # 4C fit to the initial e+e- four-momentum
      nominal
      constrain_four_momentum
      chi2_cut 200
  }

alg_highE
  .note(:mass_window, "the Lambda (anti-Lambda) candidates are selected with the mass window
    |M(p pi) - M_Lambda| < 0.01 GeV/c^2, where M(p pi) is obtained from the track parameters
    after the secondary vertex fit; the window is applied before the 4C kinematic fit.")
  .note(:opening_angle_cut, "c.m. energy dependent requirements on the opening angle between
    Lambda and anti-Lambda in the c.m. system: theta(Lambda anti-Lambda) > 170, 176 and 178
    degrees at sqrt(s) = 2.400, 2.800 and 3.080 GeV, respectively.")
  .note(:efficiency_curve, "the detection efficiency depends on the unknown |G_E/G_M| ratio
    through the Lambda angular distribution; it is evaluated with MC samples weighted by
    (1 + cos^2 theta) and (1 - cos^2 theta) (the |G_E| = 0 and |G_M| = 0 extremes) and the
    nominal efficiency is their average. The associated systematic uncertainty is 10.8% - 12.7%.")
  .note(:background_veto, "non-Lambda background is studied from the two-dimensional sideband
    1.084 < M(p pi) < 1.104 GeV/c^2; the Lambda peaking background from
    e+e- -> Sigma0 anti-Sigma0, Lambda anti-Sigma0, Xi0 anti-Xi0 is found to be negligible.")
  .note(:radiative_correction, "the radiative correction factor 1 + delta (ISR photon emission,
    beam energy spread and vacuum polarization) is evaluated from the ConExc generator log.")
  .note(:event_counting, "the number of observed signal events is taken as the number of entries
    in the range 0.98 < M(Lambda anti-Lambda)/sqrt(s) < 1.02.")
  .with_decay_card(decay_card_conexc)
  .apply(selection_highE)

# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------
root_files_modeI  = alg_modeI.execute_on([data_2232, incMC_2232, exMC_modeI])
root_files_modeII = alg_modeII.execute_on([data_2232, incMC_2232, exMC_modeII])
root_files_highE  = alg_highE.execute_on([data_2400, incMC_2400, data_2800, incMC_2800,
                                          data_3080, incMC_3080] + exMCs_highE)
