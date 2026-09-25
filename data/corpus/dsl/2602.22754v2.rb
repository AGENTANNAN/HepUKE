# Analysis: Measurements of branching fractions of
#   Lambda_c+ -> Sigma0 K_S0 pi+  and  Lambda_c+ -> Sigma0 K_S0 K+
# using 6.4 fb-1 of e+e- data collected with BESIII at 13 CMS energies
# ranging from 4.600 to 4.950 GeV.  Sigma0 -> gamma Lambda, Lambda -> p pi-,
# K_S0 -> pi+ pi-.  Signal reconstructed as single-tag of Lambda_c+ candidates,
# the recoiling Lambda_c-bar being constrained via a 4C kinematic fit
# (recoil mass fixed to the Lambda_c-bar nominal mass along with K_S0,
#  Lambda and Sigma0 intermediate masses).

### Dataset description ###
# All 13 high-energy XYZ points from the BESIII XYZ program
xyz_datasets = [
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("706_4610"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4740"),
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
  DatasetManager.real_data.find("707_4840"),
  DatasetManager.real_data.find("707_4914"),
  DatasetManager.real_data.find("707_4946"),
]

xyz_incMCs = [
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
  DatasetManager.inclusive_mc.find("707_4740"),
  DatasetManager.inclusive_mc.find("707_4750"),
  DatasetManager.inclusive_mc.find("707_4780"),
  DatasetManager.inclusive_mc.find("707_4840"),
  DatasetManager.inclusive_mc.find("707_4914"),
  DatasetManager.inclusive_mc.find("707_4946"),
]

# ---------------------------------------------------------------------------
# Signal decay cards (PHSP model for the Lambda_c+ signal three-body decay)
# ---------------------------------------------------------------------------
decay_card_SigKspi = <<~DECAYCARD
    Decay Lambda_c+
    1.0000 Sigma0 K_S0 pi+                    PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-Sigma0 K_S0 pi-               PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0                      PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 gamma anti-Lambda0                 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                             PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                        PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                            PHSP;
    Enddecay

    End
DECAYCARD

decay_card_SigKsK = <<~DECAYCARD
    Decay Lambda_c+
    1.0000 Sigma0 K_S0 K+                     PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-Sigma0 K_S0 K-                PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0                      PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 gamma anti-Lambda0                 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                             PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                        PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                            PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for each mode is generated per energy point via
# create_exclusive_mc_for; the e+e- -> Lambda_c+ Lambda_c-bar production
# uses the measured cross section (default).
exMC_SigKspi = DatasetManager.create_exclusive_mc_for(xyz_datasets) do |config|
  config.sample_name   = "LcpToSigma0KsPi_signalMC"
  config.events        = 200000
  config.decay_card    = decay_card_SigKspi
  config.cross_section = :default
end
exMC_SigKspi.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

exMC_SigKsK = DatasetManager.create_exclusive_mc_for(xyz_datasets) do |config|
  config.sample_name   = "LcpToSigma0KsK_signalMC"
  config.events        = 200000
  config.decay_card    = decay_card_SigKsK
  config.cross_section = :default
end
exMC_SigKsK.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

# ===========================================================================
# Algorithm 1 : Lambda_c+ -> Sigma0 K_S0 pi+
# ===========================================================================
alg_name1 = "LcpSigma0KsPi"
alg_SigKspi = Algorithm.new(alg_name1)
alg_SigKspi.set_header(["#{alg_name1}Alg/#{alg_name1}.h"])
           .set_constant({ "ECMS" => [:double, 4.681] })
           .set_alias({ "std::vector<double>" => "Vdouble" })

sel_SigKspi = Selection.new
sel_SigKspi.select_track {
              cos_theta 0.93        # |cos(theta)| < 0.93 (MDC acceptance)
              Vz        10.0        # |Vz|  < 10 cm  (non-V0 tracks only)
              Vr        1.0         # |Vxy| < 1  cm  (non-V0 tracks only)
              nChrp    ">=2"
              nChrn    ">=2"
              nNet     "==0"
            }
            .select_photon {
              tdc_emc_start     0
              tdc_emc_end       700 / 50   # keep within [0,700] ns spec
              angle_to_track    10.0
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              nGam ">=1"
            }
            .pid(method: :probability) {
              prob_cut 0.001
              identify :proton, against: [:kaon, :pion]
              identify :pion,   against: [:kaon, :proton]
              nprp ">=1"
              npip ">=1"
              npim ">=1"
            }
            # Lambda -> p pi- (secondary vertex, mass window 1.111-1.121 GeV/c^2)
            .secondary_vertex_fit([:prp, :pim]) {
              build_virtual_particle(:Lambda).by_minimizing_mass_difference
              remove_used_particle_from_candidate_list
            }
            # K_S0 -> pi+ pi- (secondary vertex, mass window 0.487-0.511 GeV/c^2)
            .secondary_vertex_fit([:pip, :pim]) {
              build_virtual_particle(:K_S0).by_minimizing_mass_difference
              remove_used_particle_from_candidate_list
            }
            # Sigma0 -> gamma Lambda (mass window 1.179-1.203 GeV/c^2).
            # Reconstruct as a virtual particle from gamma + Lambda.
            .kalman_kinematic_fit([:gamma, :Lambda]) {
              invariant_mass_of(:gamma, :Lambda).constrain_to_nominal_mass_of(:Sigma0)
              chi2_cut 200
              nSigma0 ">=1"
            }
            # Main 4C kinematic fit: constrain the total 4-momentum to ECMS,
            # constrain m(pi+pi-)/m(p pi-)/m(gamma p pi-) to K_S0/Lambda/Sigma0
            # nominal masses, and constrain the recoil (missing) invariant
            # mass of Lambda_c+ candidates to the Lambda_c-bar nominal mass.
            .kinematic_fit([:Sigma0, :K_S0, :pip]) {
              nominal
              constrain_four_momentum
              miss_track_of(:Lambda_c_bar)
              invariant_mass_of(:Lambda_c_bar).constrain_to_nominal_mass_of(:Lambda_c_bar)
              chi2_cut 200
            }

alg_SigKspi
  .note(:ks_lambda_track_vz,
        "Charged tracks assigned to K_S0/Lambda candidates use a looser |Vz| < 20 cm " \
        "cut and are not required to satisfy the primary-vertex IP cuts.")
  .note(:ks_lambda_vertex_chi2,
        "K_S0 and Lambda vertex fits are required to satisfy chi^2 < 100 " \
        "and the fitted decay length must exceed twice the vertex resolution.")
  .note(:mass_windows,
        "M(p pi-) in [1.111,1.121] GeV/c^2 for Lambda; " \
        "M(pi+ pi-) in [0.487,0.511] GeV/c^2 for K_S0; " \
        "M(gamma Lambda) in [1.179,1.203] GeV/c^2 for Sigma0. " \
        "3-sigma windows around the PDG masses.")
  .note(:chi2_4c_tight_cut,
        "Optimised tight requirement chi2_4C < 29 applied in ROOT after " \
        "S/sqrt(S+B) optimisation on the M_BC signal region " \
        "[2.282, 2.291] GeV/c^2.")
  .note(:mBC_definition,
        "Signal identified via beam-constrained mass " \
        "M_BC = sqrt(E_beam^2 - p_Lc^2) evaluated at each energy point.")

alg_SigKspi.with_decay_card(decay_card_SigKspi).apply(sel_SigKspi)
alg_SigKspi.execute_on(xyz_datasets + xyz_incMCs + exMC_SigKspi)

# ===========================================================================
# Algorithm 2 : Lambda_c+ -> Sigma0 K_S0 K+
# ===========================================================================
alg_name2 = "LcpSigma0KsK"
alg_SigKsK = Algorithm.new(alg_name2)
alg_SigKsK.set_header(["#{alg_name2}Alg/#{alg_name2}.h"])
          .set_constant({ "ECMS" => [:double, 4.681] })
          .set_alias({ "std::vector<double>" => "Vdouble" })

sel_SigKsK = Selection.new
sel_SigKsK.select_track {
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp    ">=2"
             nChrn    ">=2"
             nNet     "==0"
           }
           .select_photon {
             tdc_emc_start     0
             tdc_emc_end       700 / 50
             angle_to_track    10.0
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             nGam ">=1"
           }
           .pid(method: :probability) {
             prob_cut 0.001
             identify :proton, against: [:kaon, :pion]
             identify :kaon,   against: [:pion, :proton]
             identify :pion,   against: [:kaon, :proton]
             nprp ">=1"
             nkp  ">=1"
             npip ">=1"
             npim ">=1"
           }
           # Lambda -> p pi-
           .secondary_vertex_fit([:prp, :pim]) {
             build_virtual_particle(:Lambda).by_minimizing_mass_difference
             remove_used_particle_from_candidate_list
           }
           # K_S0 -> pi+ pi-
           .secondary_vertex_fit([:pip, :pim]) {
             build_virtual_particle(:K_S0).by_minimizing_mass_difference
             remove_used_particle_from_candidate_list
           }
           # Sigma0 -> gamma Lambda
           .kalman_kinematic_fit([:gamma, :Lambda]) {
             invariant_mass_of(:gamma, :Lambda).constrain_to_nominal_mass_of(:Sigma0)
             chi2_cut 200
             nSigma0 ">=1"
           }
           # Main 4C kinematic fit for the Sigma0 K_S0 K+ mode
           .kinematic_fit([:Sigma0, :K_S0, :kp]) {
             nominal
             constrain_four_momentum
             miss_track_of(:Lambda_c_bar)
             invariant_mass_of(:Lambda_c_bar).constrain_to_nominal_mass_of(:Lambda_c_bar)
             chi2_cut 200
           }

alg_SigKsK
  .note(:ks_lambda_track_vz,
        "Charged tracks assigned to K_S0/Lambda candidates use a looser " \
        "|Vz| < 20 cm cut and are not required to satisfy the primary-vertex " \
        "IP cuts.")
  .note(:ks_lambda_vertex_chi2,
        "K_S0 and Lambda vertex fits are required to satisfy chi^2 < 100 " \
        "and the fitted decay length must exceed twice the vertex resolution.")
  .note(:mass_windows,
        "M(p pi-) in [1.111,1.121] GeV/c^2 for Lambda; " \
        "M(pi+ pi-) in [0.487,0.511] GeV/c^2 for K_S0; " \
        "M(gamma Lambda) in [1.179,1.203] GeV/c^2 for Sigma0.")
  .note(:chi2_4c_tight_cut,
        "Optimised tight requirement chi2_4C < 171 applied in ROOT after " \
        "S/sqrt(S+B) optimisation on the M_BC signal region " \
        "[2.282, 2.291] GeV/c^2.")
  .note(:mBC_definition,
        "Signal identified via beam-constrained mass " \
        "M_BC = sqrt(E_beam^2 - p_Lc^2) evaluated at each energy point.")

alg_SigKsK.with_decay_card(decay_card_SigKsK).apply(sel_SigKsK)
alg_SigKsK.execute_on(xyz_datasets + xyz_incMCs + exMC_SigKsK)
