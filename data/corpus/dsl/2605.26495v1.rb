### Dataset description ###
jpsi_data   = DatasetManager.real_data.find("708_3097")
jpsi_incMC  = DatasetManager.inclusive_mc.find("708_3097")

# ---------- Decay card for mode I: J/psi -> gamma K_S0 K_S0 pi0 ----------
decay_card_ksks_pi0 = <<~DECAYCARD
    Decay J/psi
    1.0000  gamma K_S0 K_S0 pi0                PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+ pi-                            PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma gamma                        PHSP;
    Enddecay

    End
DECAYCARD

# ---------- Decay card for mode II: J/psi -> gamma pi0 pi0 eta ----------
decay_card_pipi_eta = <<~DECAYCARD
    Decay J/psi
    1.0000  gamma pi0 pi0 eta                  PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma gamma                        PHSP;
    Enddecay

    Decay eta
    1.0000  gamma gamma                        PHSP;
    Enddecay

    End
DECAYCARD

exMC_ksks_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_ksks_pi0"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_ksks_pi0
  config.cross_section   = :default
end
exMC_ksks_pi0.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_pipi_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_pipi_eta"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_pipi_eta
  config.cross_section   = :default
end
exMC_pipi_eta.save_to_config(format: :yaml, file_path: 'temp_for_test')

########################################################################
# Algorithm I: J/psi -> gamma K_S0 K_S0 pi0
########################################################################
alg_ksks = Algorithm.new("X2370KsKsPi0")
alg_ksks.set_header(["X2370KsKsPi0Alg/X2370KsKsPi0.h"])
        .set_constant({"ECMS" => [:double, 3.097]})

sel_ksks = Selection.new
sel_ksks.select_track {                    # Charged tracks come from K_S0 decays only
            cos_theta 0.93                 # |cos(theta)| < 0.93
            nNet      "==0"                # net charge zero
          }
        .select_photon {
            tdc_emc_start      0
            tdc_emc_end        14          # EMC time within 700 ns
            angle_to_track     10.0        # opening angle to nearest track > 10 deg
            energyThreshold_b  0.025       # barrel:   E > 25 MeV (|cos(theta)| < 0.80)
            energyThreshold_e  0.050       # end cap:  E > 50 MeV (0.86 < |cos(theta)| < 0.92)
            nGam               ">=3"       # at least three photons
          }
        .assign({:chrgp => :pip, :chrgn => :pim})   # all charged tracks assumed to be pions
        # Reconstruct first K_S0 -> pi+ pi- via secondary vertex fit
        .secondary_vertex_fit([:pip, :pim]) {
            build_virtual_particle(:K_S0).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
        # Reconstruct second K_S0 -> pi+ pi-
        .secondary_vertex_fit([:pip, :pim]) {
            build_virtual_particle(:K_S0).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
        # Kalman fit: reconstruct pi0 -> gamma gamma  (5C = 4C+1C step will use it)
        .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 200
            npi0 ">=1"
          }
        # Main kinematic fit: J/psi -> gamma_rad K_S0 K_S0 pi0  (4C + pi0 mass constraint = 5C)
        .kinematic_fit([:gamma, :K_S0, :K_S0, :pi0]) {
            nominal
            constrain_four_momentum
            chi2_cut 200                   # loose cut; tight 40 applied in ROOT
          }

alg_ksks.note(:ks_selection,
              "|M(pi+ pi-) - m_KS0| < 11 MeV/c^2 and decay length > 2 * vertex resolution; " \
              "applied inside the secondary-vertex-fit stage.")
        .note(:pi0_mass_window,
              "|M(gamma gamma) - m_pi0| < 22 MeV/c^2 and each pi0 photon E > 100 MeV.")
        .note(:chi2_4c_cut,
              "chi2 of 4C kinematic fit J/psi -> 3gamma 2K_S0 required < 40.")
        .note(:chi2_5c_cut,
              "chi2 of 5C kinematic fit (adds pi0 mass constraint) used to pick best pi0 combination.")
        .note(:background_veto,
              "reject events with |M(gamma_rad gamma_pi0) - m_pi0| < 22 MeV/c^2, " \
              "|M(gamma_rad gamma_pi0) - m_eta| < 25 MeV/c^2, and " \
              "|M(gamma_rad pi0) - m_omega| < 40 MeV/c^2 to suppress mis-combinations and eta/omega backgrounds.")

alg_ksks.with_decay_card(decay_card_ksks_pi0).apply(sel_ksks)
alg_ksks.execute_on([jpsi_data, jpsi_incMC, exMC_ksks_pi0])

########################################################################
# Algorithm II: J/psi -> gamma pi0 pi0 eta
########################################################################
alg_pipi_eta = Algorithm.new("X2370PiPiEta")
alg_pipi_eta.set_header(["X2370PiPiEtaAlg/X2370PiPiEta.h"])
            .set_constant({"ECMS" => [:double, 3.097]})

sel_pipi_eta = Selection.new
sel_pipi_eta.select_track {                # no charged tracks in the final state
              cos_theta 0.93
              nNet      "==0"
              nChrp     "==0"
              nChrn     "==0"
            }
            .select_photon {
              tdc_emc_start      0
              tdc_emc_end        14
              angle_to_track     10.0
              energyThreshold_b  0.025
              energyThreshold_e  0.050
              nGam               ">=7"     # at least seven photons
            }
            # 1C fit: reconstruct pi0 candidates (chi2_1C < 10, at least 2 pi0)
            .kalman_kinematic_fit([:gamma, :gamma]) {
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
              chi2_cut 10
              npi0 ">=2"
            }
            # 7C fit: reconstruct eta candidate from gamma pair
            .kalman_kinematic_fit([:gamma, :gamma]) {
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
              chi2_cut 200
              neta ">=1"
            }
            # Main 6C -> 7C kinematic fit: J/psi -> gamma_rad pi0 pi0 eta
            .kinematic_fit([:gamma, :pi0, :pi0, :eta]) {
              nominal
              constrain_four_momentum
              chi2_cut 200                 # loose cut; tight 40 applied in ROOT
            }

alg_pipi_eta.note(:chi2_1c_cut,
                  "For pi0 candidates the 1C fit chi^2 must be < 10; at least two pi0 pairs required.")
            .note(:chi2_6c_cut,
                  "chi2 of 6C kinematic fit J/psi -> 3gamma 2pi0 required < 40.")
            .note(:eta_mass_window,
                  "|M(gamma gamma) - m_eta| < 27 MeV/c^2 for the eta candidate.")
            .note(:photon_hypothesis_veto,
                  "chi2_4C(J/psi -> 7 gamma) must be smaller than chi2_4C(J/psi -> 8 gamma) " \
                  "and chi2_4C(J/psi -> 9 gamma) to suppress multi-photon backgrounds.")
            .note(:background_veto,
                  "reject events with |M(gamma_rad gamma_pi0) - m_pi0| < 20 MeV/c^2, " \
                  "|M(gamma_rad gamma_eta) - m_pi0| < 20 MeV/c^2, " \
                  "|M(gamma_rad gamma_pi0) - m_eta| < 30 MeV/c^2, " \
                  "|M(gamma_rad gamma_eta) - m_eta| < 50 MeV/c^2, and " \
                  "|M(gamma_rad pi0) - m_omega| < 40 MeV/c^2 to suppress photon " \
                  "mis-combinations and eta/omega backgrounds.")

alg_pipi_eta.with_decay_card(decay_card_pipi_eta).apply(sel_pipi_eta)
alg_pipi_eta.execute_on([jpsi_data, jpsi_incMC, exMC_pipi_eta])
