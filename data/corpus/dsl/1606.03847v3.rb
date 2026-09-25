# Amplitude analysis of the isospin-violating decays
# eta' -> pi+ pi- pi0 and eta' -> pi0 pi0 pi0 produced in J/psi -> gamma eta'
# (BESIII, 1.31e9 J/psi events).  The two eta' decay modes have different
# final states and different selection chains, so (Rule T1) each mode gets its
# own Algorithm object and Selection chain.

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi dataset at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Matching inclusive MC sample

# Decay card for Mode I: J/psi -> gamma eta', eta' -> pi+ pi- pi0
decay_card_etap_pipimpi0 = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'       PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- pi0      PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma      PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for Mode II: J/psi -> gamma eta', eta' -> pi0 pi0 pi0
decay_card_etap_3pi0 = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'       PHSP;
  Enddecay

  Decay eta'
  1.0000 pi0 pi0 pi0      PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma      PHSP;
  Enddecay

  End
DECAYCARD

exMC_etap_pipimpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etap_pipimpi0"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_etap_pipimpi0
  config.cross_section   = :default
end

exMC_etap_3pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etap_3pi0"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_etap_3pi0
  config.cross_section   = :default
end

### Event selection (BOSS): Mode I  J/psi -> gamma eta', eta' -> pi+ pi- pi0 ###
alg_name_modeI = "JpsiGammaEtaPrimePiPiPi0"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

# Two charged tracks of opposite charge and at least three photon candidates.
# The photon with the largest energy in the event is the radiative photon from
# J/psi -> gamma eta'; the best photon assignment is the one minimising the
# chi2 of the kinematic fit below.
selection_modeI = Selection.new
selection_modeI
  .select_track {
     cos_theta 0.93
     Vz        10.0
     Vr        1.0
     nChrp     "==1"
     nChrn     "==1"
     nNet      "==0"
  }
  .select_photon {
     tdc_emc_start     0
     tdc_emc_end       14
     angle_to_track    10.0
     energyThreshold_b 0.025
     energyThreshold_e 0.050
     nGam              ">=3"
  }
  .pid(method: :probability) {
     prob_cut 0.001
     identify :pion, against: [:kaon, :proton]
     npip "==1"
     npim "==1"
  }
  # pi0 reconstruction from photon pairs, M(gamma gamma) constrained to m_pi0
  # (1C Kalman fit, chi2 < 25).  The best photon pairing is taken automatically.
  .kalman_kinematic_fit([:gamma, :gamma]) {
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
     chi2_cut 25
     npi0 ">=1"
  }
  # Nominal 6C kinematic fit: 4C energy-momentum conservation, the
  # M(gamma gamma) = m_pi0 constraint (already imposed by the Kalman fit above)
  # and M(pi+ pi- pi0) = m_eta'.  The combination with the smallest chi2_6C is
  # retained automatically by the fit.
  .kinematic_fit([:gamma, :pip, :pim, :pi0]) {
     nominal
     constrain_four_momentum
     invariant_mass_of(:pip, :pim, :pi0).constrain_to_nominal_mass_of(:etap)
     chi2_cut 200
  }
  # Competing hypothesis A (Rule T2): 4C fit under the J/psi -> pi+ pi- gamma gamma
  # gamma signal hypothesis; the chi2 is stored for the ROOT-level probability
  # ordering veto against the background hypotheses.
  .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {
     constrain_four_momentum
  }
  # Competing hypothesis B (Rule T2): 4C fit under the J/psi -> pi+ pi- gamma gamma
  # background hypothesis (two-photon final state).
  .kinematic_fit([:gamma, :gamma, :pip, :pim]) {
     constrain_four_momentum
  }

alg_modeI
  .note(:radiative_photon_assignment,
        "the radiative photon from J/psi -> gamma eta' is the photon candidate with " \
        "the maximum energy in the event; the remaining two photons form the pi0. " \
        "When more than three photons are present, the combination with the smallest " \
        "chi2_6C is retained (handled by the chi2-minimising combination selection of " \
        "the kinematic fit).")
  .note(:background_veto,
        "omega background from J/psi -> omega pi+ pi- suppressed by rejecting " \
        "events with |M(gamma pi0) - m_omega| < 0.05 GeV/c^2, where gamma is the " \
        "radiative photon and pi0 the reconstructed pi+pi-pi0 daughter.")
  .note(:photon_multiplicity_veto,
        "Four-photon background hypothesis J/psi -> pi+ pi- gamma gamma gamma gamma " \
        "tested with a 4C fit; the event is discarded if the 4C fit probability of " \
        "the three-photon signal hypothesis is not larger than that of the two- and " \
        "four-photon hypotheses (probability ordering applied in the ROOT analysis).")
  .note(:chi2_cut_published,
        "Published selection requires chi2_6C < 25; kept loose (200) in BOSS and " \
        "applied at the optimal value in ROOT.")
  .with_decay_card(decay_card_etap_pipimpi0)
  .apply(selection_modeI)

### Event selection (BOSS): Mode II  J/psi -> gamma eta', eta' -> pi0 pi0 pi0 ###
alg_name_modeII = "JpsiGammaEtaPrime3Pi0"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})

# No charged tracks and at least seven photon candidates.  One photon is the
# radiative photon from J/psi -> gamma eta'; the remaining six photons form the
# three pi0 candidates.
selection_modeII = Selection.new
selection_modeII
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
     nGam              ">=7"
  }
  # pi0 reconstruction: each photon pair mass-constrained to m_pi0
  # (1C Kalman fit, chi2 < 25); at least three pi0 candidates required.
  .kalman_kinematic_fit([:gamma, :gamma]) {
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
     chi2_cut 25
     npi0 ">=3"
  }
  # eta reconstruction, needed for the J/psi -> gamma eta pi0 pi0 veto fit below.
  .kalman_kinematic_fit([:gamma, :gamma]) {
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
     chi2_cut 25
     neta ">=1"
  }
  # Nominal 8C kinematic fit: 4C energy-momentum conservation, three
  # M(gamma gamma) = m_pi0 constraints (from the Kalman fit above) and
  # M(pi0 pi0 pi0) = m_eta'.  If more than one combination survives, the one
  # with the smallest chi2_8C is retained by the fit itself.
  .kinematic_fit([:gamma, :pi0, :pi0, :pi0]) {
     nominal
     constrain_four_momentum
     invariant_mass_of(:pi0, :pi0, :pi0).constrain_to_nominal_mass_of(:etap)
     chi2_cut 200
  }
  # Competing hypothesis (Rule T2): 7C fit under the J/psi -> gamma eta pi0 pi0
  # hypothesis (4C + m_eta + 2 x m_pi0, the mass constraints coming from the
  # Kalman fits above).  The event is discarded later in ROOT if this hypothesis
  # has a larger probability than the signal hypothesis.
  .kinematic_fit([:gamma, :eta, :pi0, :pi0]) {
     constrain_four_momentum
  }

alg_modeII
  .note(:background_veto,
        "two vetoes applied: events containing at least one gamma gamma pair with " \
        "invariant mass in the eta signal region (0.52, 0.59) GeV/c^2 are rejected, " \
        "and events with |M(gamma pi0) - m_omega| < 0.05 GeV/c^2 are rejected to " \
        "suppress J/psi -> omega pi0 pi0.")
  .note(:background_veto_eta_pi0pi0,
        "J/psi -> gamma eta pi0 pi0 background suppressed by discarding events for " \
        "which the probability of the 7C fit under this hypothesis exceeds that of " \
        "the signal hypothesis (probability ordering applied in the ROOT analysis).")
  .note(:pi0_decay_angle,
        "pi0 candidates required to satisfy |cos theta_decay| < 0.95, where " \
        "theta_decay is the polar angle of a photon in the gamma gamma rest frame; " \
        "used to suppress pi0 mis-combinations.")
  .note(:chi2_cut_published,
        "Published selection requires chi2_8C < 70; kept loose (200) in BOSS and " \
        "applied at the optimal value in ROOT.")
  .with_decay_card(decay_card_etap_3pi0)
  .apply(selection_modeII)

### Execute the algorithms on data, inclusive MC and signal exclusive MC ###
root_files_modeI = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_etap_pipimpi0])
root_files_modeII = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_etap_3pi0])
