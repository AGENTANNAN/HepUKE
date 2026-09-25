# ============================================================================
#  BESIII  J/psi(3097) -> gamma eta
#  Absolute eta decay branching fractions from four exclusive eta decay modes
# ============================================================================

### Dataset description ###
# Real data: 1.0087 x 10^10 J/psi(3097) events
jpsi_data  = DatasetManager.real_data.find("708_3097")
# Corresponding inclusive MC at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ------------------------- Decay cards (EvtGen) -----------------------------
# J/psi -> gamma eta ; eta -> gamma gamma
decay_card_eta_gg = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# J/psi -> gamma eta ; eta -> pi0 pi0 pi0
decay_card_eta_3pi0 = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta PHSP;
  Enddecay

  Decay eta
  1.000 pi0 pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# J/psi -> gamma eta ; eta -> pi+ pi- pi0
decay_card_eta_pipimpi0 = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta PHSP;
  Enddecay

  Decay eta
  1.000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# J/psi -> gamma eta ; eta -> pi+ pi- gamma
decay_card_eta_pipimg = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta PHSP;
  Enddecay

  Decay eta
  1.000 pi+ pi- gamma PHSP;
  Enddecay

  End
DECAYCARD

# ------------------- Exclusive signal MC (1M events each) -------------------
exMC_eta_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_eta_gg"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_eta_gg
  config.cross_section   = :default
end

exMC_eta_3pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_eta_3pi0"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_eta_3pi0
  config.cross_section   = :default
end

exMC_eta_pipimpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_eta_pipimpi0"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_eta_pipimpi0
  config.cross_section   = :default
end

exMC_eta_pipimg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_eta_pipimg"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_eta_pipimg
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ======================= Mode I : eta -> gamma gamma ========================
alg_gg_name = "JpsiEtaToGG"
alg_gg = Algorithm.new(alg_gg_name)
alg_gg.set_header(["#{alg_gg_name}Alg/#{alg_gg_name}.h"])
      .set_constant({"ECMS" => [:double, 3.097]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      # The inclusive J/psi -> gamma eta channel with a converted radiative
      # photon (gamma -> e+e- via the Photon Conversion Finder) is an extra
      # channel whose reconstruction is not expressible in the DSL.
      .note(:photon_conversion_finder,
            "Inclusive J/psi -> gamma eta channel in which the radiative photon " \
            "converts, gamma -> e+e-; the conversion is reconstructed with the " \
            "Photon Conversion Finder. This conversion-reconstruction selection " \
            "is not expressible in the DSL and must be implemented directly in BOSS.")

sel_gg = Selection.new
sel_gg.select_track {                 # require no charged tracks
         nChrp "==0"
         nChrn "==0"
         nNet  "==0"
       }
       .select_photon {               # >= 3 barrel photons with E > 0.07 GeV
         tdc_emc_start 0
         tdc_emc_end 14
         angle_to_track 10.0
         energyThreshold_b 0.07
         nGam ">=3"
       }
       .kinematic_fit([:gamma, :gamma, :gamma]) {  # 4C fit to gamma gamma gamma
         nominal
         constrain_four_momentum                   # 4C energy-momentum constraint
         chi2_cut 80
       }

alg_gg.with_decay_card(decay_card_eta_gg).apply(sel_gg)
alg_gg.execute_on([jpsi_data, jpsi_incMC, exMC_eta_gg])

# ======================= Mode II : eta -> pi0 pi0 pi0 =======================
alg_3pi0_name = "JpsiEtaTo3Pi0"
alg_3pi0 = Algorithm.new(alg_3pi0_name)
alg_3pi0.set_header(["#{alg_3pi0_name}Alg/#{alg_3pi0_name}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_3pi0 = Selection.new
sel_3pi0.select_track {                # require no charged tracks
          nChrp "==0"
          nChrn "==0"
          nNet  "==0"
        }
        .select_photon {               # >= 6 photons, E > 0.025 (barrel) / 0.050 (endcap)
          tdc_emc_start 0
          tdc_emc_end 14
          angle_to_track 10.0
          energyThreshold_b 0.025
          energyThreshold_e 0.050
          nGam ">=6"
        }
        .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C mass constraint -> pi0
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          chi2_cut 25
          npi0 ">=3"
        }
        .kinematic_fit([:gamma, :pi0, :pi0, :pi0]) {  # 7C: 4C + three pi0 mass constraints
          nominal
          constrain_four_momentum                     # 4C energy-momentum constraint
          chi2_cut 100                                # best combination = smallest chi2
        }

alg_3pi0.with_decay_card(decay_card_eta_3pi0).apply(sel_3pi0)
alg_3pi0.execute_on([jpsi_data, jpsi_incMC, exMC_eta_3pi0])

# ======================= Mode III : eta -> pi+ pi- pi0 ======================
alg_pipimpi0_name = "JpsiEtaToPiPiPi0"
alg_pipimpi0 = Algorithm.new(alg_pipimpi0_name)
alg_pipimpi0.set_header(["#{alg_pipimpi0_name}Alg/#{alg_pipimpi0_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

sel_pipimpi0 = Selection.new
sel_pipimpi0.select_track {            # exactly one pi+ and one pi-
              cos_theta 0.93
              Vz 10.0
              Vr 1.0
              nChrp "==1"
              nChrn "==1"
              nNet  "==0"
            }
            .select_photon {           # >= 3 photons, E > 0.025 (barrel) / 0.050 (endcap)
              tdc_emc_start 0
              tdc_emc_end 14
              angle_to_track 10.0
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              nGam ">=3"
            }
            .pid(method: :probability) {     # pi/K separation, probability method
              prob_cut 0.001
              identify :pion, against: [:kaon, :proton]
              npip "==1"
              npim "==1"
            }
            .kinematic_fit([:pip, :pim, :gamma, :gamma]) {  # 5C: 4C + pi0 mass constraint
              nominal
              constrain_four_momentum                       # 4C energy-momentum constraint
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
              chi2_cut 100
            }

alg_pipimpi0.with_decay_card(decay_card_eta_pipimpi0).apply(sel_pipimpi0)
alg_pipimpi0.execute_on([jpsi_data, jpsi_incMC, exMC_eta_pipimpi0])

# ======================= Mode IV : eta -> pi+ pi- gamma =====================
alg_pipimg_name = "JpsiEtaToPiPiG"
alg_pipimg = Algorithm.new(alg_pipimg_name)
alg_pipimg.set_header(["#{alg_pipimg_name}Alg/#{alg_pipimg_name}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_pipimg = Selection.new
sel_pipimg.select_track {              # one pi+ and one pi-
            cos_theta 0.93
            Vz 10.0
            Vr 1.0
            nChrp "==1"
            nChrn "==1"
            nNet  "==0"
          }
          .select_photon {             # >= 2 photons, E > 0.025 (barrel) / 0.050 (endcap)
            tdc_emc_start 0
            tdc_emc_end 14
            angle_to_track 10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam ">=2"
          }
          .pid(method: :probability) {   # pi/K separation, probability method
            prob_cut 0.001
            identify :pion, against: [:kaon, :proton]
            npip "==1"
            npim "==1"
          }
          .kinematic_fit([:pip, :pim, :gamma, :gamma]) {  # nominal 4C fit (gamma pi+ pi- gamma)
            nominal
            constrain_four_momentum
            chi2_cut 60
          }
          # Competing non-nominal 5C fit under the pi+ pi- pi0 hypothesis,
          # re-using the same tracks; its chi2 is only stored so that the veto
          # (prob(gamma pi+ pi- gamma) > prob(gamma pi+ pi- pi0)) is applied in ROOT.
          .kinematic_fit([:pip, :pim, :gamma, :gamma]) {
            use_track_index_from_nominal_kmfit
            constrain_four_momentum
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          }

alg_pipimg.with_decay_card(decay_card_eta_pipimg).apply(sel_pipimg)
alg_pipimg.execute_on([jpsi_data, jpsi_incMC, exMC_eta_pipimg])