# Core DSL classes and dependencies are loaded automatically at execution time.

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/ψ(3097) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # matching inclusive MC

# ------------------------------------------------------------------
# Decay cards (EvtGen format)
# ------------------------------------------------------------------
# Mode I: J/psi -> p pbar eta, eta -> gamma gamma
decay_card_gg = <<~DECAYCARD
  Decay J/psi
  1.000 p+ anti-p- eta PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: J/psi -> p pbar eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_3pi = <<~DECAYCARD
  Decay J/psi
  1.000 p+ anti-p- eta PHSP;
  Enddecay

  Decay eta
  1.000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ------------------------------------------------------------------
# Exclusive MC: 1,000,000 events for each eta decay mode
# ------------------------------------------------------------------
exMC_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_ppbar_eta_gg"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_gg
  config.cross_section   = :default
end

exMC_3pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_ppbar_eta_3pipi0"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_3pi
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ==================================================================
# Mode I: J/psi -> p pbar eta, eta -> gamma gamma
# ==================================================================
alg_name_gg = "JpsiPPbarEtaGG"
alg_gg = Algorithm.new(alg_name_gg)
alg_gg.set_header(["#{alg_name_gg}Alg/#{alg_name_gg}.h"])
      .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy (GeV)

selection_gg = Selection.new
selection_gg.select_track {                 # charged-track quality + multiplicity
                cos_theta 0.93              # |cos(theta)| < 0.93
                Vz        10.0              # |Vz| < 10 cm
                Vr        1.0               # Vr < 1 cm
                nChrp     ">=1"             # at least one positive track
                nChrn     ">=1"             # at least one negative track
                nNet      "==0"             # net charge zero
              }
            .select_photon {                # photon selection
                tdc_emc_start     0         # TDC window 0-700 ns
                tdc_emc_end       14
                angle_to_track    20.0      # > 20 deg from nearest charged track
                energyThreshold_b 0.025     # 25 MeV (barrel)
                energyThreshold_e 0.050     # 50 MeV (endcap)
                nGam              ">=2"     # at least two photons
              }
            .pid(method: :probability) {
                prob_cut 0.001
                identify :proton, against: [:kaon]   # p / pbar identified against K
                nprp "==1"                  # exactly one proton
                nprm "==1"                  # exactly one anti-proton
              }
            # 4C kinematic fit to p pbar eta (eta == gamma gamma, eta mass unconstrained)
            .kinematic_fit([:prp, :prm, :gamma, :gamma]) {
                nominal                     # mark as nominal fit
                constrain_four_momentum     # 4C energy-momentum constraint
                invariant_mass_of(:gamma, :gamma).within(0.200, 0.900)  # eta -> gamma gamma mass window
                chi2_cut 200                # loose chi2 cut; smallest-chi2 combination kept by default
              }

alg_gg.with_decay_card(decay_card_gg).apply(selection_gg)
root_files_gg = alg_gg.execute_on([jpsi_data, jpsi_incMC, exMC_gg])

# ==================================================================
# Mode II: J/psi -> p pbar eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma
# ==================================================================
alg_name_3pi = "JpsiPPbarEta3PiPi0"
alg_3pi = Algorithm.new(alg_name_3pi)
alg_3pi.set_header(["#{alg_name_3pi}Alg/#{alg_name_3pi}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})

selection_3pi = Selection.new
selection_3pi.select_track {
                cos_theta 0.93
                Vz        10.0
                Vr        1.0
                nChrp     ">=2"             # at least two positive tracks
                nChrn     ">=2"             # at least two negative tracks
                nNet      "==0"
              }
             .select_photon {
                tdc_emc_start     0
                tdc_emc_end       14
                angle_to_track    20.0
                energyThreshold_b 0.025
                energyThreshold_e 0.050
                nGam              ">=2"
              }
             .pid(method: :probability) {
                prob_cut 0.001
                identify :proton, against: [:pion]    # p / pbar identified against pi (no kaon requirement)
                identify :pion,   against: [:proton]  # pi / pi identified against p
                nprp "==1"                  # exactly one proton
                nprm "==1"                  # exactly one anti-proton
                npip ">=1"                  # at least one pi+
                npim ">=1"                  # at least one pi-
              }
             .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma (1-C mass constraint)
                invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                chi2_cut 25
                npi0 ">=1"
              }
             # 4C kinematic fit to p pbar pi+ pi- pi0 (eta mass unconstrained)
             .kinematic_fit([:prp, :prm, :pip, :pim, :pi0]) {
                nominal
                constrain_four_momentum
                chi2_cut 200
              }

alg_3pi.with_decay_card(decay_card_3pi).apply(selection_3pi)
root_files_3pi = alg_3pi.execute_on([jpsi_data, jpsi_incMC, exMC_3pi])