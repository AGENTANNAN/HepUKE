# ============================================================================
# J/psi -> gamma eta' eta'  at 3.097 GeV  (two independent eta' decay modes)
# BOSS-side specification: datasets + event selection up to the final kmfit.
# ============================================================================

### Datasets ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # Real J/psi data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # Inclusive MC (10.09e9 J/psi events)

### Decay cards (EvtGen) ###
# Mode I: J/psi -> gamma eta' eta', BOTH eta' -> eta pi+ pi-, eta -> gamma gamma
#         final state 5 gamma 2 pi+ 2 pi-
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' eta' PHSP;
    Enddecay

    Decay eta'
    1.000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: J/psi -> gamma eta' eta', one eta' -> gamma pi+ pi-, the other eta' -> eta pi+ pi-
#          eta -> gamma gamma ; final state 4 gamma 2 pi+ 2 pi-
# (alias separates the two eta' with different decay chains in the same card)
decay_card_modeII = <<~DECAYCARD
    Alias another_etap eta'

    Decay J/psi
    1.0000 gamma eta' another_etap PHSP;
    Enddecay

    Decay eta'
    1.000 gamma pi+ pi- PHSP;
    Enddecay

    Decay another_etap
    1.000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC (500k events each) ###
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3097_g_etap_etap_modeI"
    config.related_dataset = jpsi_data
    config.events          = 500000
    config.decay_card      = decay_card_modeI
    config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3097_g_etap_etap_modeII"
    config.related_dataset = jpsi_data
    config.events          = 500000
    config.decay_card      = decay_card_modeII
    config.cross_section   = :default
end

# ============================================================================
# Mode I  (5 gamma 2 pi+ 2 pi-)
# ============================================================================
alg_name_modeI = "JpsiGammaEtapEtapModeI"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})

sel_modeI = Selection.new
    .select_track {                 # Charged-track quality cuts
        cos_theta 0.93              # |cos(theta)| < 0.93
        Vz        10.0              # |Vz| < 10 cm
        Vr        1.0               # Vr < 1 cm
        nChrp     "==2"             # exactly two positive tracks
        nChrn     "==2"             # exactly two negative tracks
        nNet      "==0"             # net charge zero
    }
    .select_photon {                # Photon selection
        tdc_emc_start     0         # EMC TDC 0-14
        tdc_emc_end       14
        angle_to_track    10.0      # angle to nearest charged track > 10 deg
        energyThreshold_b 0.025     # barrel threshold 25 MeV
        energyThreshold_e 0.050     # endcap threshold 50 MeV
        nGam              ">=5"     # at least five photons (Mode I final state)
    }
    .pid(method: :probability) {    # PID: probability method
        prob_cut 0.001              # prob > 0.001
        identify :pion, against: [:kaon, :proton]   # all charged tracks as pions
        npip "==2"                  # 2 pi+
        npim "==2"                  # 2 pi-
    }
    # 1C Kalman fit to reconstruct eta -> gamma gamma (at least two eta, chi2 < 25)
    .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=2"
    }
    # Nominal 6C fit : 4C four-momentum constraint + the two (Kalman-reconstructed,
    # mass-constrained) eta; loose chi2 cut (tight cut optimised in ROOT)
    .kinematic_fit([:gamma, :eta, :eta, :pip, :pim, :pip, :pim]) {
        nominal
        constrain_four_momentum
        # best eta' pair selected by minimising sum |M(eta pi+pi-) - m_eta'|^2 (ROOT side);
        # the corresponding mass window is |M(eta pi+pi-) - m_eta'| < 0.01 GeV
        invariant_mass_of(:eta, :pip, :pim).within(0.94778, 0.96778)
        chi2_cut 85                 # chi2 < 85
    }
    # Competing 4C fits for the photon-number veto (4 gamma / 5 gamma / 6 gamma hypotheses)
    .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :pip, :pim]) {
        constrain_four_momentum
    }
    .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma, :pip, :pim, :pip, :pim]) {
        constrain_four_momentum
    }
    .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma, :gamma, :pip, :pim, :pip, :pim]) {
        constrain_four_momentum
    }

alg_modeI
    .note(:efficiency_curve, "the two eta' are chosen as the eta pi+ pi- pair minimising the sum of |M(eta pi+pi-) - m_eta'|^2; the optimisation is performed at the ROOT level on the stored candidates")
    .note(:background_veto, "photon-number veto: accept the 5-gamma hypothesis only if chi2_5gamma < chi2_4gamma and chi2_5gamma < chi2_6gamma; the competing 4C fits store chi2_4gamma and chi2_6gamma, the comparison is applied in the ROOT analysis")
    .with_decay_card(decay_card_modeI)
    .apply(sel_modeI)

root_files_modeI = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])

# ============================================================================
# Mode II  (4 gamma 2 pi+ 2 pi-)
# ============================================================================
alg_name_modeII = "JpsiGammaEtapEtapModeII"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_modeII = Selection.new
    .select_track {                 # Charged-track quality cuts
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"
        nChrn     "==2"
        nNet      "==0"
    }
    .select_photon {                # Photon selection
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=4"     # at least four photons (Mode II final state)
    }
    .pid(method: :probability) {    # PID: probability method
        prob_cut 0.001
        identify :pion, against: [:kaon, :proton]
        npip "==2"
        npim "==2"
    }
    # 1C Kalman fit to reconstruct eta -> gamma gamma (at least one eta, chi2 < 25)
    .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=1"
    }
    # Nominal 5C fit : 4C four-momentum constraint + one eta mass constraint;
    # loose chi2 cut (tight cut optimised in ROOT)
    .kinematic_fit([:gamma, :gamma, :eta, :pip, :pim, :pip, :pim]) {
        nominal
        constrain_four_momentum
        # eta' from eta pi+ pi- : |M - m_eta'| < 0.01 GeV
        invariant_mass_of(:eta, :pip, :pim).within(0.94778, 0.96778)
        # eta' from gamma pi+ pi- : |M - m_eta'| < 0.02 GeV
        invariant_mass_of(:gamma, :pip, :pim).within(0.93778, 0.97778)
        # M(pi+ pi-) in [0.4, 0.85] GeV for the gamma pi+ pi- eta'
        invariant_mass_of(:pip, :pim).within(0.4, 0.85)
        # pi0 veto : |M(gamma gamma) - m_pi0| > 0.02 GeV
        invariant_mass_of(:gamma, :gamma).out_of(0.11498, 0.15498)
        chi2_cut 75                 # chi2 < 75
    }
    # Competing 4C fits for the photon-number veto (3 gamma / 4 gamma / 5 gamma hypotheses)
    .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim, :pip, :pim]) {
        constrain_four_momentum
    }
    .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :pip, :pim]) {
        constrain_four_momentum
    }
    .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma, :pip, :pim, :pip, :pim]) {
        constrain_four_momentum
    }

alg_modeII
    .note(:efficiency_curve, "the eta' pair is chosen by minimising the weighted sum of |M - m_eta'|^2 over the gamma pi+ pi- and eta pi+ pi- combinations; the minimisation is performed at the ROOT level")
    .note(:background_veto, "pi0 veto is required for all gamma gamma pairs EXCLUDING the pair forming the eta; the pair-exclusion cannot be expressed in the fit block, the gamma gamma mass and the eta mass constraint are stored for the ROOT-level veto")
    .note(:background_veto, "photon-number veto: accept the 4-gamma hypothesis by comparing chi2_4gamma against chi2_3gamma and chi2_5gamma from the competing 4C fits, applied in the ROOT analysis")
    .note(:efficiency_curve, "M(eta' eta') < 3.0 GeV is required; eta' are composite candidates, the cut is applied on the stored eta' four-momenta in ROOT")
    .with_decay_card(decay_card_modeII)
    .apply(sel_modeII)

root_files_modeII = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])