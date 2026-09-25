# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi real data at sqrt(s) = 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive J/psi MC sample

# Decay card for Mode 1: J/psi -> Sigma+ anti-Sigma-, Sigma+ -> p pi0, tag anti-Sigma- -> anti-p- pi0
decay_card_mode1 = <<~DECAYCARD
    Decay J/psi
    1.0000 Sigma+ anti-Sigma-      PHSP;
    Enddecay

    Decay Sigma+
    1.0000 p+ pi0                  PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000 anti-p- pi0             PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma             PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for Mode 2: J/psi -> Sigma+ anti-Sigma-, Sigma+ -> n pi+, tag anti-Sigma- -> anti-p- pi0
decay_card_mode2 = <<~DECAYCARD
    Decay J/psi
    1.0000 Sigma+ anti-Sigma-      PHSP;
    Enddecay

    Decay Sigma+
    1.0000 n0 pi+                  PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000 anti-p- pi0             PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma             PHSP;
    Enddecay

    End
DECAYCARD

# Two 200k-event exclusive MC samples (one per Sigma+ decay mode), both with the Sigma- tag
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "jpsi_sigmaP_to_p_pi0_sigmaMtag_exmc"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_mode1
  config.cross_section   = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "jpsi_sigmaP_to_n_pip_sigmaMtag_exmc"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_mode2
  config.cross_section   = :default
end

exMC_mode1.save_to_config(format: :yaml, file_path: 'temp_for_test')
exMC_mode2.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
# ==================== Mode 1: Sigma+ -> p pi0 ====================
alg_name_mode1 = "SigmaPToPPi0SigmaTag"
alg_mode1 = Algorithm.new(alg_name_mode1)
alg_mode1.set_header(["#{alg_name_mode1}Alg/#{alg_name_mode1}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})   # sqrt(s) = 3.097 GeV (J/psi)
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_mode1 = Selection.new
sel_mode1.select_track {              # charged track selection
            cos_theta 0.93            # |cos(theta)| < 0.93
            Vz        10.0            # |Vz| < 10 cm
            Vr        2.0             # Vr < 2 cm
            nChrp     ">=2"           # at least two positive tracks
            nChrn     ">=1"           # at least one negative track
            nNet      ">=0"           # net charge >= 0
          }
         .select_photon {             # photon selection
            tdc_emc_start     0       # TDC 0-14
            tdc_emc_end       14
            angle_to_track    10.0    # > 10 deg from any charged track
            energyThreshold_b 0.025   # >= 25 MeV (barrel)
            energyThreshold_e 0.050   # >= 50 MeV (endcap)
            nGam              ">=2"   # at least two photons
          }
         .pid(method: :probability) { # probability-method PID
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]  # p / pbar vs K and pi
            identify :pion,   against: [:kaon]         # pi vs K
            nprp ">=1"                # at least one proton (signal side of Sigma+ -> p pi0)
            nprm ">=1"                # at least one antiproton (tag side of anti-Sigma- -> anti-p- pi0)
          }
         .select_isolated_photon {    # reject photon showers from the antiproton
            angle_to_prm_track 20.0   # > 20 deg from the antiproton candidate
            nGam ">=2"
          }
         .kalman_kinematic_fit([:gamma, :gamma]) {  # reconstruct pi0 from two photons
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"                # at least one pi0
          }
         .kinematic_fit([:prp, :prm, :pi0]) {  # nominal fit: form the Sigma- tag (pbar pi0) + extra proton
            nominal
            constrain_four_momentum
            invariant_mass_of(:prm, :pi0).within(1.169, 1.205)  # Sigma- tag mass window; best combination kept by fit
            chi2_cut 200
          }

# The signal pi0 (from Sigma+ -> p pi0) is not reconstructed: its recoil/missing mass
# window cannot be expressed with the current DSL (recoil_mass_of not implemented).
alg_mode1.note(:missing_mass_window,
  "Sigma+ -> p pi0 (double tag): after removing the tag (anti-p- pi0) and the extra proton, " \
  "the unreconstructed signal pi0 recoil/missing mass M(J/psi) - (p anti-p- pi0) must lie in " \
  "[0.034, 0.231] GeV (3-sigma window around the pi0 mass). Not expressible in the BOSS DSL.")

alg_mode1.with_decay_card(decay_card_mode1).apply(sel_mode1)

# ==================== Mode 2: Sigma+ -> n pi+ ====================
alg_name_mode2 = "SigmaPToNPiPSigmaTag"
alg_mode2 = Algorithm.new(alg_name_mode2)
alg_mode2.set_header(["#{alg_name_mode2}Alg/#{alg_name_mode2}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})   # sqrt(s) = 3.097 GeV (J/psi)
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_mode2 = Selection.new
sel_mode2.select_track {              # charged track selection
            cos_theta 0.93            # |cos(theta)| < 0.93
            Vz        10.0            # |Vz| < 10 cm
            Vr        2.0             # Vr < 2 cm
            nChrp     ">=1"           # at least one positive track
            nChrn     ">=1"           # at least one negative track
            nNet      "==0"           # net charge zero
          }
         .select_photon {             # photon selection
            tdc_emc_start     0       # TDC 0-14
            tdc_emc_end       14
            angle_to_track    10.0    # > 10 deg from any charged track
            energyThreshold_b 0.025   # >= 25 MeV (barrel)
            energyThreshold_e 0.050   # >= 50 MeV (endcap)
            nGam              ">=2"   # at least two photons
          }
         .pid(method: :probability) { # probability-method PID
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]  # p / pbar vs K and pi
            identify :pion,   against: [:kaon]         # pi vs K
            nprm ">=1"                # at least one antiproton (tag side of anti-Sigma- -> anti-p- pi0)
            npip ">=1"                # at least one pi+ (signal side of Sigma+ -> n pi+)
          }
         .select_isolated_photon {    # reject photon showers from the antiproton
            angle_to_prm_track 20.0   # > 20 deg from the antiproton candidate
            nGam ">=2"
          }
         .kalman_kinematic_fit([:gamma, :gamma]) {  # reconstruct pi0 from two photons
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"                # at least one pi0
          }
         .kinematic_fit([:pip, :prm, :pi0]) {  # nominal fit: form the Sigma- tag (pbar pi0) + extra pi+
            nominal
            constrain_four_momentum
            invariant_mass_of(:prm, :pi0).within(1.169, 1.205)  # Sigma- tag mass window; best combination kept by fit
            chi2_cut 200
          }

# The neutron (from Sigma+ -> n pi+) is not detected: its recoil/missing mass window
# cannot be expressed with the current DSL (recoil_mass_of not implemented).
alg_mode2.note(:missing_mass_window,
  "Sigma+ -> n pi+ (double tag): after removing the tag (anti-p- pi0) and the extra pi+, " \
  "the undetected neutron recoil/missing mass M(J/psi) - (pi+ anti-p- pi0) must lie in " \
  "[0.881, 0.997] GeV (3-sigma window around the neutron mass). Not expressible in the BOSS DSL.")

alg_mode2.with_decay_card(decay_card_mode2).apply(sel_mode2)

# ==================== Execution ====================
root_files_mode1 = alg_mode1.execute_on([jpsi_data, jpsi_incMC, exMC_mode1])
root_files_mode2 = alg_mode2.execute_on([jpsi_data, jpsi_incMC, exMC_mode2])