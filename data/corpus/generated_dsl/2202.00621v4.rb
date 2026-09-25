### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi(3097) real data (10.09e9 J/psi)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC

# ---- Decay cards (EvtGen format) ----
# Mode I : J/psi -> gamma eta eta',  eta' -> eta pi+ pi-,  eta -> gamma gamma   (5 gamma pi+ pi-)
decay_card_modeI = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta eta' PHSP;
  Enddecay

  Decay eta'
  1.000 eta pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: J/psi -> gamma eta eta',  eta' -> gamma pi+ pi-, eta -> gamma gamma  (4 gamma pi+ pi-)
decay_card_modeII = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta eta' PHSP;
  Enddecay

  Decay eta'
  1.000 gamma pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- Exclusive MC: 500k events for each eta' mode ----
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_eta_etap_to_eta_pipi"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end
exMC_modeI.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_eta_etap_to_gamma_pipi"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end
exMC_modeII.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) — Mode I: eta' -> eta pi+ pi- ###
alg_name_I = "JpsiGammaEtaEtapEtaPiPi"
alg_modeI = Algorithm.new(alg_name_I)
alg_modeI.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})   # sqrt(s) = 3.097 GeV
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
sel_modeI.select_track {              # charged track selection
           cos_theta 0.93             # |cos(theta)| < 0.93
           Vz        10.0             # |Vz| < 10 cm
           Vr        1.0              # Vr < 1 cm
           nChrp     "==1"            # exactly one positive track
           nChrn     "==1"            # exactly one negative track
           nNet      "==0"            # net charge zero
         }
        .select_photon {              # photon selection
           tdc_emc_start     0        # EMC TDC window 0-14
           tdc_emc_end       14
           angle_to_track    10.0     # > 10 deg from nearest charged track
           energyThreshold_b 0.025    # barrel 25 MeV
           energyThreshold_e 0.050    # endcap 50 MeV
           nGam              ">=5"    # at least 5 photons (5 gamma pi+ pi-)
         }
        .pid(method: :probability) {  # pion PID (no lepton identification)
           prob_cut 0.001
           identify :pion, against: [:kaon, :proton]   # pi+ and pi- vs K and p
           npip "==1"
           npim "==1"
         }
        .kalman_kinematic_fit([:gamma, :gamma]) {      # reconstruct eta from gamma gamma
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
           chi2_cut 25
           neta ">=2"                                  # at least two eta in this mode
         }
        .kinematic_fit([:gamma, :eta, :eta, :pip, :pim]) {  # 5C fit to gamma eta eta pi+ pi-
           nominal
           constrain_four_momentum                     # four-momentum conservation
           invariant_mass_of(:eta, :pip, :pim).within(0.90, 1.01)  # eta' -> eta pi+ pi- window (best combination by chi2)
           chi2_cut 200                                # loose cut; tight cut applied in ROOT
         }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

### Event selection (BOSS) — Mode II: eta' -> gamma pi+ pi- ###
alg_name_II = "JpsiGammaEtaEtapGammaPiPi"
alg_modeII = Algorithm.new(alg_name_II)
alg_modeII.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
sel_modeII.select_track {
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
            nGam              ">=4"    # at least 4 photons (4 gamma pi+ pi-)
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :pion, against: [:kaon, :proton]
            npip "==1"
            npim "==1"
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {     # reconstruct eta from gamma gamma
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
            chi2_cut 25
            neta ">=1"                                  # at least one eta in this mode
          }
          .kinematic_fit([:gamma, :gamma, :eta, :pip, :pim]) {  # 4C fit to gamma gamma eta pi+ pi-
            nominal
            constrain_four_momentum                     # four-momentum conservation
            invariant_mass_of(:gamma, :pip, :pim).within(0.90, 1.01)  # eta' -> gamma pi+ pi- window (best combination by chi2)
            chi2_cut 200                                # loose cut; tight cut applied in ROOT
          }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

# ---- Execute on real data, inclusive MC and the two signal exclusive MC samples ----
root_files_modeI  = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])
root_files_modeII = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])