# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
# J/psi (3.097 GeV) real data (2.253e8 events) and the corresponding inclusive MC.
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ---- Decay cards (EvtGen format) for the three eta' decay modes ----
# Mode I: eta' -> pi+ pi- e+ e-
decay_card_modeI = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'   PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- e+ e-   PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: eta' -> pi+ pi- mu+ mu-
decay_card_modeII = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'   PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- mu+ mu-   PHSP;
  Enddecay

  End
DECAYCARD

# Mode III: eta' -> gamma pi+ pi-
decay_card_modeIII = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'   PHSP;
  Enddecay

  Decay eta'
  1.0000 gamma pi+ pi-   PHSP;
  Enddecay

  End
DECAYCARD

# ---- Exclusive MC: 100k events for each of the three eta' decay modes ----
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etap_pipi_ee"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etap_pipi_mumu"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

exMC_modeIII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etap_gammapipi"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_modeIII
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Independent final states / PID hypotheses -> one Algorithm + Selection per mode (Rule T1).

# ===== Mode I: J/psi -> gamma eta', eta' -> pi+ pi- e+ e- =====
alg_name_modeI = "GammaEtapToPipPimee"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_modeI = Selection.new
sel_modeI.select_track {                 # charged track selection
           cos_theta 0.93               # |cos(theta)| < 0.93
           Vz        10.0               # |Vz| < 10 cm
           Vr        1.0                # Vr < 1 cm
           nChrp     "==2"              # exactly 2 positive tracks
           nChrn     "==2"              # exactly 2 negative tracks
           nNet      "==0"              # net charge zero
         }
         .select_photon {                # photon selection
           tdc_emc_start     0          # EMC timing window
           tdc_emc_end       14
           energyThreshold_b 0.025      # 25 MeV (barrel)
           energyThreshold_e 0.050      # 50 MeV (endcap)
           angle_to_track    15.0       # >= 15 deg from any good track
           nGam              ">=1"      # at least 1 photon
         }
         .pid(method: :chi2_sum) {       # combinatorial chi2-sum PID: pi vs e per charge
           chi_min_cut 4                # chi2_min < 4
           identify :pion, :electron    # two positive {pi+,e+}, two negative {pi-,e-}
         }
         .kinematic_fit([:gamma, :pip, :pim, :ep, :em]) {   # 4C fit under gamma pi+ pi- e+ e-
           nominal                                          # nominal fit; best photon/track combination via min chi2
           constrain_four_momentum
           chi2_cut 75
         }
         .kinematic_fit([:gamma, :pip, :pim, :ep, :em]) {   # competing-hypothesis 4C fit (Rule T2): store chi2 only
           constrain_four_momentum                          # no chi2_cut, not nominal -> chi2 kept for ROOT-level veto
         }
alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

# ===== Mode II: J/psi -> gamma eta', eta' -> pi+ pi- mu+ mu- =====
alg_name_modeII = "GammaEtapToPipPimumu"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})

sel_modeII = Selection.new
sel_modeII.select_track {
            cos_theta 0.93
            Vz        10.0
            Vr        1.0
            nChrp     "==2"
            nChrn     "==2"
            nNet      "==0"
          }
          .select_photon {
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            angle_to_track    15.0
            nGam              ">=1"
          }
          .pid(method: :chi2_sum) {       # combinatorial chi2-sum PID: pi vs mu per charge
            chi_min_cut 4
            identify :pion, :muon         # two positive {pi+,mu+}, two negative {pi-,mu-}
          }
          .kinematic_fit([:gamma, :pip, :pim, :mup, :mum]) {   # 4C fit under gamma pi+ pi- mu+ mu-
            nominal
            constrain_four_momentum
            chi2_cut 75
          }
          .kinematic_fit([:gamma, :pip, :pim, :mup, :mum]) {   # competing-hypothesis 4C fit: store chi2 only
            constrain_four_momentum
          }
alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

# ===== Mode III: J/psi -> gamma eta', eta' -> gamma pi+ pi- =====
alg_name_modeIII = "GammaEtapToGammaPipi"
alg_modeIII = Algorithm.new(alg_name_modeIII)
alg_modeIII.set_header(["#{alg_name_modeIII}Alg/#{alg_name_modeIII}.h"])
           .set_constant({"ECMS" => [:double, 3.097]})

sel_modeIII = Selection.new
sel_modeIII.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     "==2"
             nChrn     "==2"
             nNet      "==0"
           }
           .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track    15.0
             nGam              ">=2"      # at least 2 photons (radiative gamma + eta' gamma)
           }
           .assign({:chrgp => :pip, :chrgn => :pim})   # no PID in this mode: treat tracks as pions
           .kinematic_fit([:gamma, :gamma, :pip, :pim]) {   # 4C fit under pi+ pi- gamma gamma
             nominal
             constrain_four_momentum
             chi2_cut 75
             invariant_mass_of(:gamma, :gamma).larger_than(0.16)  # M(gamma gamma) > 0.16 GeV/c^2 to veto pi0
           }
# Combination choice by eta' mass proximity (rather than fit chi2) is not expressible in the BOSS DSL.
alg_modeIII.note(:combination_selection, "among the available gamma pi+ pi- combinations the one whose invariant mass is closest to the nominal eta' mass (0.9578 GeV/c^2) is retained; this mass-based combination choice cannot be expressed in the BOSS selection (the 4C fit selects by chi2) and is applied downstream in ROOT")
alg_modeIII.with_decay_card(decay_card_modeIII).apply(sel_modeIII)

### Execute on real data, inclusive MC, and the three exclusive MC samples ###
root_files_modeI   = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])
root_files_modeII  = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])
root_files_modeIII = alg_modeIII.execute_on([jpsi_data, jpsi_incMC, exMC_modeIII])