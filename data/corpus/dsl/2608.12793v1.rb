### Dataset description ###
psip3770_data  = DatasetManager.real_data.find("712_3773")     # 20.3 fb-1 at sqrt(s) = 3.773 GeV
psip3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # Inclusive MC at psi(3770)

# ------------------------------------------------------------------------------
# Signal decay cards: e+e- -> e+e- gamma gamma* -> e+e- eta'
# Mode I:  eta' -> pi+ pi- gamma
# Mode II: eta' -> pi+ pi- eta, eta -> gamma gamma
# ------------------------------------------------------------------------------
decay_card_signal_modeI = <<~DECAYCARD
    Decay psi(3770)
    1.0000  e+   e-   eta'                        PHSP;
    Enddecay

    Decay eta'
    1.0000  pi+  pi-  gamma                       ETA_DALITZ;
    Enddecay

    End
DECAYCARD

decay_card_signal_modeII = <<~DECAYCARD
    Decay psi(3770)
    1.0000  e+   e-   eta'                        PHSP;
    Enddecay

    Decay eta'
    1.0000  pi+  pi-  eta                         ETAPRIME_DALITZ;
    Enddecay

    Decay eta
    1.0000  gamma  gamma                          PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples
exMC_signal_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_ee_etap_pipigamma"
  config.related_dataset = psip3770_data
  config.events          = 500000
  config.decay_card      = decay_card_signal_modeI
  config.cross_section   = :default
end

exMC_signal_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_ee_etap_pipieta"
  config.related_dataset = psip3770_data
  config.events          = 500000
  config.decay_card      = decay_card_signal_modeII
  config.cross_section   = :default
end

# ==============================================================================
# Mode I:  eta' -> pi+ pi- gamma
# ==============================================================================
alg_modeI = Algorithm.new("EtapTFFModeI")
alg_modeI.set_header(["EtapTFFModeIAlg/EtapTFFModeI.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
sel_modeI.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nTot      "==3"          # exactly three charged tracks (net charge = +/-1)
           }
          .select_photon {
             energyThreshold_b 0.025  # >25 MeV in barrel  (|cos theta| < 0.80)
             energyThreshold_e 0.050  # >50 MeV in end cap (0.86 < |cos theta| < 0.92)
             tdc_emc_start     0
             tdc_emc_end       14     # shower time in [0, 700] ns
             angle_to_track    10.0   # >10 deg to closest charged-track extrapolation
             nGam              ">=1"  # >= 1 good photon for mode I
           }
          .pid(method: :probability) {
             identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.5,
                                            treat_as_electron_if_energy_above: 0.8
             prob_cut 0.001
             identify :pion, against: [:kaon, :proton]  # pion PID: L(pi) largest
             npip "==1"                # exactly one pi+
             npim "==1"                # exactly one pi-
             # remaining track = tagged e+ or e- (identified via E/p > 0.8 by lepton block)
           }
          # 1C global fit: gamma pi+ pi- e_tag e_miss  (missing lepton with e mass hypothesis)
          .kinematic_fit([:gamma, :pip, :pim, :ep, :em]) {
             nominal
             constrain_four_momentum
             miss_track_of(:em)        # untagged lepton missing (charge conjugate handled internally)
             chi2_cut 200              # loose cut in BOSS; tight chi2_1C < 70 applied in ROOT
           }

alg_modeI
  .note(:tag_configuration,
        "Single-tag two-photon fusion: exactly three charged tracks with net " \
        "charge +/-1 matching the tagged lepton charge. The untagged lepton " \
        "escapes acceptance and is treated as a missing e+/e- in the 1C fit.")
  .note(:untagged_lepton_angle,
        "cos(theta_miss) > 0.99 applied on the fitted missing-lepton direction " \
        "in ROOT to ensure Q^2_miss ~ 0 (quasi-real photon).")

alg_modeI.with_decay_card(decay_card_signal_modeI).apply(sel_modeI)
root_files_modeI = alg_modeI.execute_on([psip3770_data, psip3770_incMC, exMC_signal_modeI])

# ==============================================================================
# Mode II:  eta' -> pi+ pi- eta,  eta -> gamma gamma
# ==============================================================================
alg_modeII = Algorithm.new("EtapTFFModeII")
alg_modeII.set_header(["EtapTFFModeIIAlg/EtapTFFModeII.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
sel_modeII.select_track {
              cos_theta 0.93
              Vz        10.0
              Vr        1.0
              nTot      "==3"
            }
           .select_photon {
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              tdc_emc_start     0
              tdc_emc_end       14
              angle_to_track    10.0
              nGam              ">=2"  # >= 2 good photons for eta -> gamma gamma
            }
           .pid(method: :probability) {
              identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.5,
                                             treat_as_electron_if_energy_above: 0.8
              prob_cut 0.001
              identify :pion, against: [:kaon, :proton]
              npip "==1"
              npim "==1"
            }
           # 1C Kalman fit to reconstruct eta from photon pairs
           .kalman_kinematic_fit([:gamma, :gamma]) {
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
              chi2_cut 200
              neta ">=1"
            }
           # 2C global fit: pip pim eta e_tag e_miss with M(gg) = m_eta constraint
           .kinematic_fit([:gamma, :gamma, :pip, :pim, :ep, :em]) {
              nominal
              constrain_four_momentum
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
              miss_track_of(:em)
              chi2_cut 200              # loose cut in BOSS; tight chi2_2C < 200 (mode II) in ROOT
            }

alg_modeII
  .note(:tag_configuration,
        "Single-tag two-photon fusion; three charged tracks with net charge " \
        "+/-1 = tag charge; untagged lepton is missing (electron-mass hypothesis).")
  .note(:untagged_lepton_angle,
        "cos(theta_miss) > 0.99 applied in ROOT on fitted variables.")

alg_modeII.with_decay_card(decay_card_signal_modeII).apply(sel_modeII)
root_files_modeII = alg_modeII.execute_on([psip3770_data, psip3770_incMC, exMC_signal_modeII])
