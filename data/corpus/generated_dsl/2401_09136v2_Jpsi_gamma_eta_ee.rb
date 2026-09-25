### Dataset description ###
# J/psi peak data at sqrt(s) = 3.097 GeV (708_3097)
jpsi_data  = DatasetManager.real_data.find("708_3097")      # real data at the J/psi peak
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

# Decay card for J/psi -> gamma eta, eta -> gamma e+ e-
# The eta -> gamma e+ e- decay is simulated by a custom generator based on theoretical
# amplitudes; PHSP is used here as the decay-card placeholder.
decay_card_eta = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for J/psi -> gamma eta', eta' -> gamma e+ e- (custom generator, PHSP placeholder)
decay_card_etap = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples, 100k events each
exMC_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_eta_gammaee"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_eta
  config.cross_section   = :default
end

exMC_etap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etap_gammaee"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_etap
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name_eta = "JpsiGammaEtaToGammaEE"
alg_eta = Algorithm.new(alg_name_eta)
alg_eta.set_header(["#{alg_name_eta}Alg/#{alg_name_eta}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})

alg_name_etap = "JpsiGammaEtaPrimeToGammaEE"
alg_etap = Algorithm.new(alg_name_etap)
alg_etap.set_header(["#{alg_name_etap}Alg/#{alg_name_etap}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})

# Common selection chain shared by the eta and eta' modes
event_selection = Selection.new
event_selection.select_track {        # charged track selection
                  cos_theta 0.93      # |cos(theta)| < 0.93
                  Vz        10.0      # |Vz| < 10 cm
                  Vr        1.0       # Vr < 1 cm
                  nChrp     "==1"     # exactly one positive track
                  nChrn     "==1"     # exactly one negative track
                  nNet      "==0"     # net charge zero
                }
               .select_photon {       # photon selection
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0     # min angle to nearest charged track (deg)
                  energyThreshold_b 0.025    # 25 MeV (barrel)
                  energyThreshold_e 0.050    # 50 MeV (endcap)
                  nGam              ">=2"    # at least two photons
                }
               .pid(method: :probability) {  # lepton identification
                  # tracks with p > 1.0 GeV treated as leptons; lepton with EMC energy > 0.6 GeV -> electron, else muon
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  nlp "==1"   # exactly one positive lepton
                  nlm "==1"   # exactly one negative lepton
                }
               .kinematic_fit([:gamma, :gamma, :lp, :lm]) {   # 4C fit to gamma gamma e+ e-
                  nominal                          # nominal fit
                  constrain_four_momentum          # 4-momentum constraint
                  # veto J/psi -> e+ e- eta/eta' with eta/eta' -> gamma gamma
                  invariant_mass_of(:gamma, :gamma).out_of(0.517, 0.577)  # +-0.03 GeV around m(eta)
                  invariant_mass_of(:gamma, :gamma).out_of(0.928, 0.988)  # +-0.03 GeV around m(eta')
                  chi2_cut 100                     # chi2 < 100
                }

# eta mode algorithm
alg_eta.note(:primary_vertex_fit,
             "events are required to pass a primary-vertex fit before charged-track selection; the fitted primary vertex is used for the track Vz/Vr cuts")
       .note(:photon_conversion_rejection,
             "photon conversions in the detector are rejected when R_xy > 2 cm and cos(theta_eg) > 0 and |Delta_xy| < 0.8 cm")
       .note(:background_veto,
             "QED background suppressed: the low-energy (Dalitz) photon is required to have E > 0.15 GeV and angle to the e+- > 10 deg, and the J/psi radiative photon is required to have angle to the e+- > 20 deg; the highest-energy photon is assigned as the J/psi radiative photon and the remaining photon as the eta Dalitz photon")
       .note(:electron_likelihood_ratio_cut,
             "electron likelihood-ratio cut L(e)/(L(e)+L(pi)) > 0.5 applied to the identified leptons for the eta mode")
       .note(:custom_generator,
             "eta -> gamma e+ e- simulated with a custom generator based on theoretical amplitudes; PHSP is used as the decay-card placeholder")
       .with_decay_card(decay_card_eta)
       .apply(event_selection.dup)

# eta' mode algorithm (same chain except the electron PID likelihood-ratio cut)
alg_etap.note(:primary_vertex_fit,
              "events are required to pass a primary-vertex fit before charged-track selection; the fitted primary vertex is used for the track Vz/Vr cuts")
        .note(:photon_conversion_rejection,
              "photon conversions in the detector are rejected when R_xy > 2 cm and cos(theta_eg) > 0 and |Delta_xy| < 0.8 cm")
        .note(:background_veto,
              "QED background suppressed: the low-energy (Dalitz) photon is required to have E > 0.15 GeV and angle to the e+- > 10 deg, and the J/psi radiative photon is required to have angle to the e+- > 20 deg; the highest-energy photon is assigned as the J/psi radiative photon and the remaining photon as the eta' Dalitz photon")
        .note(:electron_likelihood_ratio_cut,
              "electron likelihood-ratio cut L(e)/(L(e)+L(pi)) > 0.95 applied to the identified leptons for the eta' mode")
        .note(:custom_generator,
              "eta' -> gamma e+ e- simulated with a custom generator based on theoretical amplitudes; PHSP is used as the decay-card placeholder")
        .with_decay_card(decay_card_etap)
        .apply(event_selection.dup)

# Execute the algorithms on real data, inclusive MC and the corresponding exclusive MC
root_files_eta  = alg_eta.execute_on([jpsi_data, jpsi_incMC, exMC_eta])
root_files_etap = alg_etap.execute_on([jpsi_data, jpsi_incMC, exMC_etap])