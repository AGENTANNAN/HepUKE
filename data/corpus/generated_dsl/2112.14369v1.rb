# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi real data at sqrt(s) = 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive MC sample

# Decay card for Mode I: J/psi -> e+ e- pi+ pi- eta', eta' -> gamma pi+ pi-
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000 e+ e- pi+ pi- eta'  PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi-  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for Mode II: J/psi -> e+ e- pi+ pi- eta', eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.0000 e+ e- pi+ pi- eta'  PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- eta  PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC for Mode I (500k events)
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_ee_4pi_etap_gam"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

# Exclusive signal MC for Mode II (500k events)
exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_ee_4pi_etap_pipim_eta"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ------------------------------------------------------------------
# Mode I : eta' -> gamma pi+ pi-
#          4C kinematic fit to e+ e- pi+ pi- pi+ pi- gamma
# ------------------------------------------------------------------
alg_name_modeI = "JpsiEE4piEtapGam"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})   # sqrt(s) = 3.097 GeV

sel_modeI = Selection.new
sel_modeI.select_track {                    # exactly six charged tracks, |cos(theta)|, Vz, Vr, net charge
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp     "==3"                  # 1 e+ + 2 pi+
           nChrn     "==3"                  # 1 e- + 2 pi-
           nNet      "==0"
         }
         .select_photon {                    # >=1 photon (eta' -> gamma pi+ pi-)
           tdc_emc_start    0
           tdc_emc_end      14
           angle_to_track   10.0
           energyThreshold_b 0.025          # 25 MeV in the barrel
           energyThreshold_e 0.050          # 50 MeV in the endcap
           nGam ">=1"
         }
         .pid(method: :probability) {         # lepton + pion identification
           identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                          treat_as_electron_if_energy_above: 0.6  # p>1.0 -> lepton; EMC eraw>0.6 -> e, else mu
           identify :pion, against: [:kaon, :proton]  # pi+ and pi- vs K and p
           nlp   "==1"                       # 1 e+
           nlm   "==1"                       # 1 e-
           npip  "==2"                       # 2 pi+
           npim  "==2"                       # 2 pi-
         }
         .kinematic_fit([:ep, :em, :pip, :pim, :pip, :pim, :gamma]) {  # 4C fit to e+e-4pi gamma
           nominal
           constrain_four_momentum
           chi2_cut 200                      # loose BOSS-level chi2 cut (tight cut applied in ROOT)
         }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

# ------------------------------------------------------------------
# Mode II : eta' -> pi+ pi- eta, eta -> gamma gamma
#           eta reconstructed via mass-constrained Kalman fit,
#           then 5C fit (4C + eta mass constraint on the eta token)
# ------------------------------------------------------------------
alg_name_modeII = "JpsiEE4piEtapPiPiEta"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})

sel_modeII = Selection.new
sel_modeII.select_track {                   # exactly six charged tracks, |cos(theta)|, Vz, Vr, net charge
            cos_theta 0.93
            Vz        10.0
            Vr        1.0
            nChrp     "==3"
            nChrn     "==3"
            nNet      "==0"
          }
          .select_photon {                    # >=2 photons (eta -> gamma gamma)
            tdc_emc_start    0
            tdc_emc_end      14
            angle_to_track   10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam ">=2"
          }
          .pid(method: :probability) {         # same lepton + pion identification as Mode I
            identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                           treat_as_electron_if_energy_above: 0.6
            identify :pion, against: [:kaon, :proton]
            nlp   "==1"
            nlm   "==1"
            npip  "==2"
            npim  "==2"
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct eta -> gamma gamma (1-C mass constraint)
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
            chi2_cut 25
            neta ">=1"
          }
          .kinematic_fit([:ep, :em, :pip, :pim, :pip, :pim, :eta]) {  # 5C fit: 4C + eta mass constraint
            nominal
            constrain_four_momentum
            chi2_cut 200                    # loose BOSS-level chi2 cut
          }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

### Execute both algorithms on data, inclusive MC and the corresponding signal MC ###
root_files_modeI  = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])
root_files_modeII = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])