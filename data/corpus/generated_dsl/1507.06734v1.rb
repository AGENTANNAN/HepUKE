# Core DSL classes and dependencies are loaded automatically at execution time.

### Dataset description ###
# J/psi real data at sqrt(s) = 3.097 GeV (1.31e9 J/psi events) and its inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the omega-gamma mode:
#   J/psi -> gamma eta', eta' -> omega gamma, omega -> pi+ pi- pi0 (Dalitz), pi0 -> gamma gamma
decay_card_omega_gamma = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.000 omega gamma PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the omega-e+e- mode:
#   J/psi -> gamma eta', eta' -> omega e+ e- (phase space), omega -> pi+ pi- pi0 (Dalitz), pi0 -> gamma gamma
decay_card_omega_ee = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.000 omega e+ e- PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the omega-gamma mode (200k events)
exMC_omega_gamma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etaprime_omega_gamma"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_omega_gamma
  config.cross_section   = :default
end

# Exclusive MC for the omega-e+e- mode (600k events)
exMC_omega_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etaprime_omega_ee"
  config.related_dataset = jpsi_data
  config.events          = 600_000
  config.decay_card      = decay_card_omega_ee
  config.cross_section   = :default
end

### Event selection (BOSS) ###

## ---------------------------------------------------------------------------
## Mode I: J/psi -> gamma eta', eta' -> omega gamma,
##         omega -> pi+ pi- pi0, pi0 -> gamma gamma  (final state gamma gamma gamma pi+ pi-)
## ---------------------------------------------------------------------------
alg_name_omega_gamma = "JpsiGammaEtaPOmegaGamma"
alg_omega_gamma = Algorithm.new(alg_name_omega_gamma)
alg_omega_gamma.set_header(["#{alg_name_omega_gamma}Alg/#{alg_name_omega_gamma}.h"])
               .set_constant({"ECMS" => [:double, 3.097]})   # sqrt(s) = 3.097 GeV

sel_omega_gamma = Selection.new
    .select_track {                    # Charged track quality cuts
        cos_theta 0.93                 # |cos(theta)| < 0.93
        Vz        20.0                 # |Vz| < 20 cm
        Vr        2.0                  # Vr < 2 cm
        nChrp     "==1"                # Exactly one positive track
        nChrn     "==1"                # Exactly one negative track
        nNet      "==0"                # Net charge zero
    }
    .select_photon {                   # Photon selection
        tdc_emc_start     0            # EMC time window start
        tdc_emc_end       14           # EMC time window end (700 ns units)
        angle_to_track    10.0         # At least 10 deg from the nearest charged track
        energyThreshold_b 0.025        # > 25 MeV in the barrel
        energyThreshold_e 0.050        # > 50 MeV in the endcap
        nGam              ">=4"        # At least four photons
    }
    .assign({:chrgp => :pip, :chrgn => :pim})   # No PID: take the two tracks as pi+ pi-
    .kalman_kinematic_fit([:gamma, :gamma]) {   # Reconstruct pi0 from a photon pair (1-C fit)
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25                # chi2 < 25
        npi0     ">=1"             # At least one pi0 candidate
    }
    .kinematic_fit([:gamma, :pi0, :pip, :pim]) {   # 4C fit to gamma gamma gamma pi+ pi-
        nominal                    # Nominal fit; best combination chosen by smallest chi2
        constrain_four_momentum    # 4-momentum conservation against the CMS energy
        chi2_cut 200               # Loose cut; the tight chi2_4C < 80 is applied in ROOT
    }

alg_omega_gamma
    .note(:signal_selection, "Final selection applied in ROOT after the nominal 4C fit is not " \
                             "expressible in the selection chain: the radiative photon is " \
                             "required to be the highest-energy photon in the event with " \
                             "E > 1.0 GeV; the pi0 photon pair must satisfy " \
                             "|M(gamma gamma) - m_pi0| < 0.015 GeV; and the final chi2_4C < 80.")
    .with_decay_card(decay_card_omega_gamma).apply(sel_omega_gamma)

alg_omega_gamma.execute_on([jpsi_data, jpsi_incMC, exMC_omega_gamma])


## ---------------------------------------------------------------------------
## Mode II: J/psi -> gamma eta', eta' -> omega e+ e-,
##          omega -> pi+ pi- pi0, pi0 -> gamma gamma  (final state gamma gamma gamma pi+ pi- e+ e-)
## ---------------------------------------------------------------------------
alg_name_omega_ee = "JpsiGammaEtaPOmegaEE"
alg_omega_ee = Algorithm.new(alg_name_omega_ee)
alg_omega_ee.set_header(["#{alg_name_omega_ee}Alg/#{alg_name_omega_ee}.h"])
           .set_constant({"ECMS" => [:double, 3.097]})   # sqrt(s) = 3.097 GeV

sel_omega_ee = Selection.new
    .select_track {                    # Same track quality cuts, four charged tracks
        cos_theta 0.93                 # |cos(theta)| < 0.93
        Vz        20.0                 # |Vz| < 20 cm
        Vr        2.0                  # Vr < 2 cm
        nChrp     "==2"                # Exactly two positive tracks
        nChrn     "==2"                # Exactly two negative tracks
        nNet      "==0"                # Net charge zero
    }
    .select_photon {                   # Same photon quality cuts, at least three photons
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=3"        # At least three photons
    }
    .pid(method: :chi2_sum) {          # Combinatorial PID: minimise the combined TOF + dE/dx chi2
        chi_min_cut 4                  # chi2_PID < 4
        identify :pion, :electron      # positive tracks -> {pi+, e+}, negative -> {pi-, e-}; bijection per charge
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # Reconstruct pi0 from a photon pair (1-C fit)
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25                # chi2 < 25
        npi0     ">=1"             # At least one pi0 candidate
    }
    .kinematic_fit([:gamma, :pi0, :pip, :pim, :ep, :em]) {   # 4C fit to gamma gamma gamma pi+ pi- e+ e-
        nominal                    # Nominal fit; best combination chosen by smallest chi2
        constrain_four_momentum    # 4-momentum conservation against the CMS energy
        chi2_cut 200               # Loose cut; the tight chi2_4C < 80 is applied in ROOT
    }

alg_omega_ee
    .note(:signal_selection, "Final selection applied in ROOT after the nominal 4C fit: the " \
                             "radiative photon is required to be the highest-energy photon " \
                             "with E > 1.0 GeV; the pi0 photon pair must satisfy " \
                             "|M(gamma gamma) - m_pi0| < 0.015 GeV; the final chi2_4C < 80; and " \
                             "a gamma-conversion veto requires the e+ e- vertex R_xy < 2 cm.")
    .note(:pid_correction_method, "The event is ultimately selected by the smallest combined " \
                                  "chi2_4C+PID; the combinatorial PID chi2 is added to the 4C " \
                                  "fit chi2 and the final pi/e hypothesis choice is made in the " \
                                  "ROOT analysis.")
    .with_decay_card(decay_card_omega_ee).apply(sel_omega_ee)

alg_omega_ee.execute_on([jpsi_data, jpsi_incMC, exMC_omega_ee])