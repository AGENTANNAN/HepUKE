# Search for a massless BSM particle in Xi0 -> Lambda + invisible
# J/psi data ~ 1.0087e10 events; Xi0 anti-Xi0 pairs from J/psi decays.
# ST side: anti-Xi0 -> anti-Lambda pi0 (anti-Lambda -> anti-p pi+; pi0 -> gamma gamma)
# Signal (DT) side: Xi0 -> Lambda + invisible; Lambda -> p pi-
# This is a hand-built single-tag/double-tag analysis, NOT a DTag-based one.
# We express it as an ordinary Algorithm with selection reproducing DT reconstruction.

### Dataset ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Signal MC: J/psi -> Xi0(-> Lambda + invisible) anti-Xi0(-> anti-Lambda pi0)
decay_card_signal = <<~DECAYCARD
    Alias  Xi0inv     Xi0
    Alias  anti-Xi0tag  anti-Xi0

    Decay J/psi
    1.000  Xi0inv  anti-Xi0tag                    PHSP;
    Enddecay

    Decay Xi0inv
    1.000  Lambda0  gamma                         PHSP;
    Enddecay

    Decay anti-Xi0tag
    1.000  anti-Lambda0  pi0                      PHSP;
    Enddecay

    Decay Lambda0
    1.000  p+  pi-                                PHSP;
    Enddecay

    Decay anti-Lambda0
    1.000  anti-p-  pi+                           PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                           PHSP;
    Enddecay

    End
DECAYCARD

# ST MC for anti-Xi0 -> anti-Lambda pi0 tag efficiency (Xi0 -> Lambda anything)
decay_card_st = <<~DECAYCARD
    Decay J/psi
    1.000  Xi0  anti-Xi0                          PHSP;
    Enddecay

    Decay anti-Xi0
    1.000  anti-Lambda0  pi0                      PHSP;
    Enddecay

    Decay anti-Lambda0
    1.000  anti-p-  pi+                           PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                           PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "Xi0_to_Lambda_invisible"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'exMC_Xi0_Lambda_inv_config')

exMC_st = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "Jpsi_Xi0_antiXi0_ST_MC"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_st
  config.cross_section   = :default
end
exMC_st.save_to_config(format: :yaml, file_path: 'exMC_Xi0_ST_config')

### Event Selection ###
alg = Algorithm.new("Xi0LambdaInv")
alg.set_header(["Xi0LambdaInvAlg/Xi0LambdaInv.h"])
   .set_constant({"ECMS" => [:double, 3.097]})

sel = Selection.new
sel.select_track {
       cos_theta 0.93     # |cos(theta)| < 0.93
       Vz        20.0     # |Vz| < 20 cm (loose along beam)
       nTot      ">=2"    # At least 2 charged tracks
       nNet      "==0"    # zero net charge
     }
    .select_photon {
       tdc_emc_start     0
       tdc_emc_end       14
       angle_to_track    10.0
       angle_to_prm_track 20.0   # 20 deg from anti-proton track
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam              ">=2"
     }
    .pid(method: :probability) {
       prob_cut 0.0
       identify :proton, against: [:pion, :kaon]
       identify :pion,   against: [:kaon, :proton]
       nprp ">=1"  # proton on DT side
       nprm ">=1"  # anti-proton on ST side
       npip ">=1"  # pi+ on ST side (from anti-Lambda)
       npim ">=1"  # pi- on DT side (from Lambda)
     }
    # Reconstruct pi0 -> gamma gamma with 1C mass constraint (chi2 < 25)
    .kalman_kinematic_fit([:gamma, :gamma]) {
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
       chi2_cut 25
       npi0 ">=1"
     }
    # Secondary vertex fit for anti-Lambda -> anti-p pi+ (ST side)
    .secondary_vertex_fit([:prm, :pip]) {
       build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
    # Secondary vertex fit for Lambda -> p pi- (DT side)
    .secondary_vertex_fit([:prp, :pim]) {
       build_virtual_particle(:Lambda).by_minimizing_mass_difference
       remove_used_particle_from_candidate_list
     }
    # Main kinematic fit: J/psi -> p anti-p pi+ pi- pi0 + invisible (massless)
    # 4-momentum conservation + missing particle inferred by kinematic fit
    .kinematic_fit([:Lambda, :Lambda_bar, :pi0]) {
       nominal
       constrain_four_momentum
       chi2_cut 200
     }

alg
  .note(:st_mBC_definition,
        "Single-tag mBC = sqrt(Ecms^2/4 - |P(anti-Lambda + pi0)|^2); anti-Xi0 tag via anti-Lambda pi0.")
  .note(:st_deltaE_best,
        "Best ST anti-Xi0 candidate chosen by minimum |E(anti-Lambda + pi0) - Ecms/2|.")
  .note(:st_deltaM_cut,
        "|M(anti-Lambda pi0) - M(anti-Xi0)_PDG| < 0.02 GeV/c^2 for anti-Xi0 candidate.")
  .note(:pi0_mass_window,
        "pi0 photon-pair invariant mass window (0.115, 0.150) GeV/c^2 (~2sigma) before 1C fit.")
  .note(:Lambda_mass_window,
        "Lambda (anti-Lambda) M(p pi) window (1.111, 1.120) GeV/c^2 (~3sigma) after secondary vertex fit.")
  .note(:Lambda_vertex_cuts,
        "Primary and secondary vertex fit chi^2 < 100 each; secondary L/sigma_L > 2 for good Lambda candidate.")
  .note(:dt_2C_kinematic_fit,
        "Two 2C kinematic fits under J/psi -> p anti-p pi+ pi- pi0 + invisible: (a) invisible massless, (b) invisible mass = M(pi0). Require chi2_2C < 8.5 AND chi2_2C_pi0 > chi2_2C to suppress Xi0 -> Lambda pi0.")
  .note(:Minv_window,
        "M(p pi- + invisible) in (1.310, 1.322) GeV/c^2.")
  .note(:cos_theta_inv_cut,
        "|cos(theta_inv)| < 0.8 to suppress Xi0 -> Lambda gamma / Lambda pi0 with missing photon.")
  .note(:p_inv_cut,
        "|P_inv| > 0.15 GeV/c to suppress J/psi -> Sigma0 anti-Sigma0.")
  .note(:extra_photon_5C_veto,
        "If N(good photons) > 2, perform 5C fit J/psi -> p anti-p pi+ pi- pi0 gamma over each extra photon; require min chi2_3gamma > 1000.")
  .note(:extra_photon_6C_veto,
        "If N(good photons) > 3, perform 6C fit J/psi -> p anti-p pi+ pi- pi0 gamma gamma over each extra photon pair; require min chi2_4gamma > 1000.")
  .note(:Eextra_shape_correction,
        "Data-driven shape correction of E_extra^other using J/psi -> anti-Xi0(anti-Lambda pi0) Xi0(Lambda pi0) control sample, binned in anti-proton momentum and polar angle.")
  .with_decay_card(decay_card_signal)
  .apply(sel)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_st])
