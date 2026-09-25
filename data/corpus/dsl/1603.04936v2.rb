# ============================================================================
# 1603.04936v2  Search for the radiative decays h_c -> gamma eta' and h_c -> gamma eta
#   psi' -> pi0 h_c, h_c -> gamma eta'(eta), with
#     eta' -> pi+ pi- eta (eta -> gamma gamma)
#     eta' -> gamma pi+ pi-
#     eta  -> gamma gamma
#     eta  -> pi+ pi- pi0 (pi0 -> gamma gamma)
# Baseline: 4.48e8 psi' events (2009 + 2012)
# ============================================================================

### Dataset description ###
psip_data   = DatasetManager.real_data.find("709_3686")     # psi(2S) data, 2009 + 2012
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")  # generic psi' inclusive MC
cont_data   = DatasetManager.real_data.find("709_3650")     # continuum data at sqrt(s) = 3.65 GeV (~44 pb^-1)

# ---------------------------------------------------------------------------
# Decay cards (one per independent reconstruction chain)
# ---------------------------------------------------------------------------
# Chain I : psi' -> pi0 h_c, h_c -> gamma eta', eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_chainI = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c HELAMP 1 0 1 0 1 0 1 0;
    Enddecay

    Decay h_c
    1.0000 gamma eta' HELAMP 1 0 1 0 1 0 1 0;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Chain II : psi' -> pi0 h_c, h_c -> gamma eta', eta' -> gamma pi+ pi-
decay_card_chainII = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c HELAMP 1 0 1 0 1 0 1 0;
    Enddecay

    Decay h_c
    1.0000 gamma eta' HELAMP 1 0 1 0 1 0 1 0;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Chain III : psi' -> pi0 h_c, h_c -> gamma eta, eta -> gamma gamma  (all-neutral final state)
decay_card_chainIII = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c HELAMP 1 0 1 0 1 0 1 0;
    Enddecay

    Decay h_c
    1.0000 gamma eta HELAMP 1 0 1 0 1 0 1 0;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Chain IV : psi' -> pi0 h_c, h_c -> gamma eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_chainIV = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c HELAMP 1 0 1 0 1 0 1 0;
    Enddecay

    Decay h_c
    1.0000 gamma eta HELAMP 1 0 1 0 1 0 1 0;
    Enddecay

    Decay eta
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive MC samples for the four signal chains
# ---------------------------------------------------------------------------
exMC_chainI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_pi0hc_hc_gamma_etap_2piEta_2gamma"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_chainI
  config.cross_section   = :default
end

exMC_chainII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_pi0hc_hc_gamma_etap_gamma2pi"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_chainII
  config.cross_section   = :default
end

exMC_chainIII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_pi0hc_hc_gamma_eta_2gamma"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_chainIII
  config.cross_section   = :default
end

exMC_chainIV = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_pi0hc_hc_gamma_eta_2piPi0"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_chainIV
  config.cross_section   = :default
end

# ===========================================================================
# Chain I : psi' -> pi0 h_c, h_c -> gamma eta', eta' -> pi+ pi- eta, eta -> gamma gamma
#   Final state: pi+ pi- (1+1) + 5 photons (pi0: 2, eta: 2, radiative: 1)
# ===========================================================================
alg_name_I = "HcToGammaEtaPrimeTo2PiEtaTo2Gamma"
alg_I = Algorithm.new(alg_name_I)
alg_I.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})

sel_I = Selection.new
sel_I.select_track {
        cos_theta 0.93        # |cos(theta)| < 0.93 within the MDC fiducial volume
        Vz        10.0        # |Vz| < 10 cm from the IP along the beam direction
        Vr        1.0         # point of closest approach < 1 cm in the radial direction
        nChrp     "==1"
        nChrn     "==1"
        nNet      "==0"
      }
     .select_photon {
        tdc_emc_start     0        # EMC time within 0 <= t <= 700 ns
        tdc_emc_end       14
        angle_to_track    10.0     # at least 10 degrees away from the nearest charged track
        energyThreshold_b 0.025    # > 25 MeV in the EMC barrel (|cos(theta)| < 0.8)
        energyThreshold_e 0.050    # > 50 MeV in the EMC endcaps (0.86 < |cos(theta)| < 0.92)
        nGam              ">=5"
      }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon, :proton]
        npip "==1"
        npim "==1"
      }
     .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      }
     .kalman_kinematic_fit([:gamma, :gamma]) {   # eta -> gamma gamma
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=1"
      }
     .kinematic_fit([:gamma, :pip, :pim, :pi0, :eta]) {
        nominal                     # 6C: 4-momentum + m(pi0) + m(eta) (the two mass
                                    # constraints are carried by the Kalman fits above)
        vertex_fit([1, 2])          # the two charged tracks must originate from the IP
        constrain_four_momentum
        chi2_cut 120                # chi2_6C < 120
      }

alg_I
  .note(:comb_least_chi2,
        "All possible photon combinations are looped over; the combination with the " \
        "least chi2_6C of the kinematic fit is selected.")
  .note(:fom_optimisation,
        "All selection criteria were optimised by maximising the figure of merit " \
        "S/sqrt(S+B), where S(B) is the number of signal (background) events in the " \
        "signal region.")
  .note(:signal_region,
        "eta' signal region: [M_eta' - 12, M_eta' + 12] MeV/c^2; sidebands " \
        "[M_eta' - 60, M_eta' - 36] and [M_eta' + 36, M_eta' + 60] MeV/c^2.")

alg_I.with_decay_card(decay_card_chainI).apply(sel_I)

# ===========================================================================
# Chain II : psi' -> pi0 h_c, h_c -> gamma eta', eta' -> gamma pi+ pi-
#   Final state: pi+ pi- (1+1) + 4 photons (pi0: 2, radiative: 1, eta': 1)
# ===========================================================================
alg_name_II = "HcToGammaEtaPrimeToGamma2Pi"
alg_II = Algorithm.new(alg_name_II)
alg_II.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})

sel_II = Selection.new
sel_II.select_track {
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
         nGam              ">=4"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :pion, against: [:kaon, :proton]
         npip "==1"
         npim "==1"
       }
       .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
         chi2_cut 25
         npi0 ">=1"
       }
       .kinematic_fit([:gamma, :gamma, :pip, :pim, :pi0]) {
         nominal                     # 5C: 4-momentum + m(pi0)
         vertex_fit([2, 3])          # the two charged tracks must originate from the IP
         constrain_four_momentum
         chi2_cut 50                 # chi2_5C < 50
       }

alg_II
  .note(:comb_least_chi2,
        "All possible photon combinations are looped over; the combination with the " \
        "least chi2_5C of the kinematic fit is selected.")
  .note(:radiative_photon_choice,
        "Of the two non-pi0 photons the one with the larger energy is taken as the " \
        "radiative photon from h_c; the other one belongs to eta' -> gamma pi+ pi-.")
  .note(:background_veto,
        "Continuum background studied on the 44 pb^-1 data set at sqrt(s) = 3.65 GeV; " \
        "no continuum event survives the selection.")
  .note(:signal_region,
        "eta' signal region: [M_eta' - 12, M_eta' + 12] MeV/c^2.")

alg_II.with_decay_card(decay_card_chainII).apply(sel_II)

# ===========================================================================
# Chain III : psi' -> pi0 h_c, h_c -> gamma eta, eta -> gamma gamma
#   Final state: 5 photons only (pi0: 2, eta: 2, radiative: 1), no charged track
# ===========================================================================
alg_name_III = "HcToGammaEtaTo2Gamma"
alg_III = Algorithm.new(alg_name_III)
alg_III.set_header(["#{alg_name_III}Alg/#{alg_name_III}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

sel_III = Selection.new
sel_III.select_track {
          cos_theta 0.93
          Vz        10.0
          Vr        1.0
          nChrp     "==0"
          nChrn     "==0"
        }
        .select_photon {
          tdc_emc_start     0
          tdc_emc_end       14        # applied only when charged particles are present
          angle_to_track    10.0
          energyThreshold_b 0.025
          energyThreshold_e 0.050
          nGam              ">=5"
        }
        .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          chi2_cut 25
          npi0 ">=1"
        }
        .kalman_kinematic_fit([:gamma, :gamma]) {   # eta -> gamma gamma
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
          chi2_cut 25
          neta ">=1"
        }
        .kinematic_fit([:gamma, :pi0, :eta]) {
          nominal                   # 6C: 4-momentum + m(pi0) + m(eta)
          constrain_four_momentum
          chi2_cut 200
        }

alg_III
  .note(:comb_least_chi2,
        "All possible photon combinations are looped over; the combination with the " \
        "least chi2_6C of the kinematic fit is selected.")
  .note(:eta_sideband_5c_fit,
        "For the selected five photons a second, 5C kinematic fit is performed " \
        "(final-state 4-momentum constrained to the total initial e+e- four-momentum " \
        "plus a pi0 mass constraint) so that the eta sideband can be used to verify " \
        "the signal; the paper requires chi2_5C < 35 for that cross-check sample. " \
        "The fit used for the nominal four-momenta is the 6C fit above.")
  .note(:signal_region,
        "eta (-> gamma gamma) signal region: [M_eta - 25, M_eta + 25] MeV/c^2, " \
        "resolution ~8 MeV/c^2; sidebands [M_eta - 100, M_eta - 50] and " \
        "[M_eta + 50, M_eta + 100] MeV/c^2.")

alg_III.with_decay_card(decay_card_chainIII).apply(sel_III)

# ===========================================================================
# Chain IV : psi' -> pi0 h_c, h_c -> gamma eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma
#   Final state: pi+ pi- (1+1) + 5 photons (two pi0's: 2+2, radiative: 1)
# ===========================================================================
alg_name_IV = "HcToGammaEtaTo2PiPi0"
alg_IV = Algorithm.new(alg_name_IV)
alg_IV.set_header(["#{alg_name_IV}Alg/#{alg_name_IV}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})

sel_IV = Selection.new
sel_IV.select_track {
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
         nGam              ">=5"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :pion, against: [:kaon, :proton]
         npip "==1"
         npim "==1"
       }
       .kalman_kinematic_fit([:gamma, :gamma]) {   # both pi0 -> gamma gamma candidates
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
         chi2_cut 25
         npi0 ">=2"
       }
       .kinematic_fit([:gamma, :pip, :pim, :pi0, :pi0]) {
         nominal                   # 6C: 4-momentum + the masses of the two pi0's
         vertex_fit([1, 2])        # the two charged tracks must originate from the IP
         constrain_four_momentum
         chi2_cut 120              # chi2_6C < 120
       }

alg_IV
  .note(:comb_least_chi2,
        "All possible photon combinations are looped over; the combination with the " \
        "least chi2_6C of the kinematic fit is selected.")
  .note(:signal_region,
        "eta (-> pi+ pi- pi0) signal region: [M_eta - 12, M_eta + 12] MeV/c^2, " \
        "resolution ~3 MeV/c^2; sidebands [M_eta - 48, M_eta - 24] and " \
        "[M_eta + 24, M_eta + 48] MeV/c^2.")

alg_IV.with_decay_card(decay_card_chainIV).apply(sel_IV)

# ===========================================================================
# Execute all four algorithms
# ===========================================================================
root_files_I   = alg_I.execute_on([psip_data, psip_incMC, cont_data, exMC_chainI])
root_files_II  = alg_II.execute_on([psip_data, psip_incMC, cont_data, exMC_chainII])
root_files_III = alg_III.execute_on([psip_data, psip_incMC, cont_data, exMC_chainIII])
root_files_IV  = alg_IV.execute_on([psip_data, psip_incMC, cont_data, exMC_chainIV])
