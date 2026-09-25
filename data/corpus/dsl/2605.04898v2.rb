### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi real data (10.087B events)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # inclusive MC for J/psi

# Decay card for Mode I: J/psi -> gamma eta, eta -> e+ e- e+ e-
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta                              PHOTOS HELAMP 1.0 0.0 -1.0 0.0;
    Enddecay
    Decay eta
    1.000 e+ e- e+ e-                            PHSP;
    Enddecay
    End
DECAYCARD

# Decay card for Mode II: J/psi -> gamma eta', eta' -> pi+ pi- eta, eta -> e+ e- e+ e-
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta'                             PHOTOS HELAMP 1.0 0.0 -1.0 0.0;
    Enddecay
    Decay eta'
    1.000 pi+ pi- eta                            PHSP;
    Enddecay
    Decay eta
    1.000 e+ e- e+ e-                            PHSP;
    Enddecay
    End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exclusive_jpsi_gamma_eta_eeee"
    config.related_dataset = jpsi_data
    config.events         = 2000000
    config.decay_card     = decay_card_modeI
    config.cross_section  = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exclusive_jpsi_gamma_etap_pipieta_eeee"
    config.related_dataset = jpsi_data
    config.events         = 2000000
    config.decay_card     = decay_card_modeII
    config.cross_section  = :default
end

### Event selection: Mode I - J/psi -> gamma eta, eta -> e+ e- e+ e- ###
alg_name_modeI = "JpsiGammaEtaEEEE"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_modeI = Selection.new
sel_modeI.select_track {
            cos_theta 0.93          # |cos(theta)| <= 0.93
            Vz        10.0
            Vr        10.0
            nChrp     "==2"         # 2 e+
            nChrn     "==2"         # 2 e-
            nNet      "==0"
          }
         .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    15.0  # angle to nearest charged track > 15 deg
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam              ">=1"
          }
         .pid(method: :probability) {
            prob_cut 0.001
            identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.0,
                                           treat_as_electron_if_energy_above: 0.0
            # electrons (L(e) > L(pi)) - identify_high_momentum_leptons handles e/mu
          }
         # 4C fit: J/psi -> gamma e+ e- e+ e- (nominal)
         .kinematic_fit([:gamma, :ep, :ep, :em, :em]) {
            nominal
            constrain_four_momentum
            chi2_cut 200            # tight chi2_4C+PID < 60 applied in ROOT
         }

alg_modeI.note(:pid_electron,
               "Electron PID: L(e) > L(pi). Combined chi2_4C+PID = chi2_4C + sum_i chi2_PID(i) " \
               "minimized to choose best combination per event; chi2_4C+PID(eeee) < 60 in ROOT.")
         .note(:gamma_conversion_veto,
               "Photon conversion veto using the photon-conversion-finder package: " \
               "veto Phi_ee < 70 deg and 2 cm < R_xy < 8 cm (mode I).")
         .note(:miscombination_veto,
               "Reject events with 10 deg < theta_ee^1 < 30 deg when 10 deg < theta_ee^2 < 60 deg " \
               "to suppress pi -> e misidentification background (mode I).")
         .note(:helix_correction,
               "Track helix parameter correction applied before kinematic fit; " \
               "half of efficiency difference w/ and w/o correction is the systematic.")

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
root_files_modeI = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])

### Event selection: Mode II - J/psi -> gamma eta', eta' -> pi+ pi- eta, eta -> e+ e- e+ e- ###
alg_name_modeII = "JpsiGammaEtapPiPiEtaEEEE"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})

sel_modeII = Selection.new
sel_modeII.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        10.0
             nChrp     "==3"         # 2 e+ + pi+
             nChrn     "==3"         # 2 e- + pi-
             nNet      "==0"
           }
          .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             angle_to_track    15.0
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             nGam              ">=1"
           }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :pion, against: [:kaon, :proton]  # pi+ / pi-
             identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.0,
                                            treat_as_electron_if_energy_above: 0.0
             npip ">=1"
             npim ">=1"
          }
          # 5C fit: J/psi -> gamma pi+ pi- e+ e- e+ e-, with pi+ pi- eta system constrained to eta' mass.
          # First step: reconstruct eta from e+ e- e+ e- (no mass constraint tokenised — the analysis
          # constrains pi+ pi- eta system to eta' mass in a 5C fit)
          .kinematic_fit([:gamma, :pip, :pim, :ep, :ep, :em, :em]) {
             nominal
             constrain_four_momentum
             invariant_mass_of(:pip, :pim, :ep, :ep, :em, :em).constrain_to_nominal_mass_of(:etap)
             chi2_cut 200            # tight chi2_5C < 60 applied in ROOT
          }

alg_modeII.note(:pid_electron,
                "Electron PID: L(e) > L(pi); pion PID: L(pi) > L(e). " \
                "Best combination per event chosen by min chi2_5C; chi2_5C < 60 in ROOT.")
          .note(:gamma_conversion_veto,
                "Photon conversion veto (mode II): Phi_ee < 40 deg and 2 cm < R_xy < 8 cm.")
          .note(:miscombination_veto,
                "Mode II: reject events with 20 deg < theta_ee^1 < 60 deg when " \
                "40 deg < theta_ee^2 < 80 deg to suppress pi -> e misidentification.")
          .note(:helix_correction,
                "Track helix parameter correction applied before kinematic fit.")

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
root_files_modeII = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])
