# Multi-channel joint analysis of T_ccbar(4020)-
# Three channels at sqrt(s) = 4.395 and 4.416 GeV:
#   1) e+e- -> D*0 D*- pi+   (partial reconstruction of D*-)
#   2) e+e- -> pi+ pi- J/psi ; J/psi -> l+ l-
#   3) e+e- -> pi+ pi- h_c   ; h_c -> gamma eta_c (16 eta_c hadronic decay modes)
# Each channel -> its own Algorithm (Rule T1).

### Datasets ###
data_4390 = DatasetManager.real_data.find("703_4390")   # 4.395 GeV (55.57 pb^-1)
data_4420 = DatasetManager.real_data.find("703_4420")   # 4.416 GeV (~1090 pb^-1)
incMC_4390 = DatasetManager.inclusive_mc.find("703_4390")
incMC_4420 = DatasetManager.inclusive_mc.find("703_4420")

######################################################################
### Channel 1: e+e- -> D*0 D*- pi+                                  ###
### D*0 -> D0 pi0 or D0 gamma ; D*- -> Dbar0 pi- (pi- not reconstructed)
### D0/Dbar0 via K- pi+ , K- pi+ pi0 , K- pi+ pi+ pi-               ###
######################################################################

# For partial reconstruction, we use decay card with the primary process; the pi- from D*- is missed
decay_card_DstarDstarPi = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*0 D*- pi+                       PHSP;
    Enddecay

    Decay D*0
    0.5   D0 pi0                            VSS;
    0.5   D0 gamma                          VSP_PWAVE;
    Enddecay

    Decay D*-
    1.0   anti-D0 pi-                       VSS;
    Enddecay

    Decay D0
    0.333 K- pi+                            PHSP;
    0.333 K- pi+ pi0                        PHSP;
    0.334 K- pi+ pi+ pi-                    PHSP;
    Enddecay

    Decay anti-D0
    0.333 K+ pi-                            PHSP;
    0.333 K+ pi- pi0                        PHSP;
    0.334 K+ pi- pi- pi+                    PHSP;
    Enddecay

    Decay pi0
    1.0   gamma gamma                       PHSP;
    Enddecay

    End
DECAYCARD

exMC_ch1 = DatasetManager.create_exclusive_mc_for([data_4390, data_4420]) do |config|
  config.sample_name    = "TccBar_DstarDstarPi_signal"
  config.events          = 500000
  config.decay_card      = decay_card_DstarDstarPi
  config.cross_section   = :default
end
exMC_ch1.each { |m| m.save_to_config(format: :yaml, file_path: "exMC_ch1_#{m.name}") }

alg_ch1 = Algorithm.new("TccDstarDstarPi")
alg_ch1.set_header(["TccDstarDstarPiAlg/TccDstarDstarPi.h"])
       .set_constant({"ECMS" => [:double, 4.416]})

sel_ch1 = Selection.new
sel_ch1.select_track {
           cos_theta 0.93
           Vr        1.0
           Vz        10.0
           nTot      ">=4"
           nNet      "==0"
         }
        .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14
           angle_to_track    10.0
           energyThreshold_b 0.025
           energyThreshold_e 0.050
         }
        .pid(method: :probability) {
           prob_cut 0.0
           identify :kaon, against: [:pion]
           identify :pion, against: [:kaon]
           nkp ">=1"
           nkm ">=1"
           npip ">=1"
           npim ">=1"
         }
        .kalman_kinematic_fit([:gamma, :gamma]) {
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 200
           npi0 ">=1"
         }
        # Partial reconstruction: pi- from D*- is missed
        .kinematic_fit([:kp, :km, :pip, :pim, :pi0, :gamma]) {
           nominal
           miss_track_of :pim
           constrain_four_momentum
           invariant_mass_of(:kp, :pim).constrain_to_nominal_mass_of(:D0)
           chi2_cut 200
         }

alg_ch1
  .note(:D0_mass_window,
        "D0 / Dbar0 candidate mass window (1.850, 1.880) GeV/c^2 for K-pi+, K-pi+pi0 and K-pi+pi+pi- modes.")
  .note(:pi0_mass_window,
        "pi0 photon-pair mass window (0.120, 0.145) GeV/c^2 for D*0 -> D0 pi0 mode; 1C mass-constraint fit to nominal pi0 mass.")
  .note(:MQ_bachelor_pi_veto,
        "MQ(D0 pi+) > 2.03 GeV/c^2 to reject pi+ from D*+ decay.")
  .note(:Dstar0_mass_windows,
        "MQ(pi0 D0) in (2.004, 2.009) GeV/c^2 for D*0 -> D0 pi0; MQ(gamma D0) in (1.995, 2.015) GeV/c^2 for D*0 -> D0 gamma.")
  .note(:recoil_mass_windows,
        "RQ(D*0 Dbar0 pi+) in (0.120, 0.160) GeV/c^2; RQ(D0 pi0 pi+) in (1.990, 2.030) GeV/c^2 or RQ(D0 gamma pi+) in (1.990, 2.040) GeV/c^2 depending on D*0 mode.")
  .note(:global_kinematic_fit,
        "Global kinematic fit constrains D0, D*0, Dbar0 to nominal masses; recoil of D*0 Dbar0 pi+ constrained to pi- mass; recoil of D*0 pi+ constrained to D*- mass. Best candidate has minimum chi^2.")
  .with_decay_card(decay_card_DstarDstarPi)
  .apply(sel_ch1)
alg_ch1.execute_on([data_4390, data_4420, incMC_4390, incMC_4420] + exMC_ch1)

######################################################################
### Channel 2: e+e- -> pi+ pi- J/psi ; J/psi -> l+ l-               ###
######################################################################
decay_card_pipiJpsi = <<~DECAYCARD
    Decay psi(4260)
    1.000 pi+ pi- J/psi                      PHSP;
    Enddecay

    Decay J/psi
    0.5   e+  e-                             PHOTOS VLL;
    0.5   mu+ mu-                            PHOTOS VLL;
    Enddecay

    End
DECAYCARD

exMC_ch2 = DatasetManager.create_exclusive_mc_for([data_4390, data_4420]) do |config|
  config.sample_name    = "TccBar_pipiJpsi_signal"
  config.events          = 500000
  config.decay_card      = decay_card_pipiJpsi
  config.cross_section   = :default
end
exMC_ch2.each { |m| m.save_to_config(format: :yaml, file_path: "exMC_ch2_#{m.name}") }

alg_ch2 = Algorithm.new("TccPipiJpsi")
alg_ch2.set_header(["TccPipiJpsiAlg/TccPipiJpsi.h"])
       .set_constant({"ECMS" => [:double, 4.416]})

sel_ch2 = Selection.new
sel_ch2.select_track {
           cos_theta 0.93
           Vr        1.0
           Vz        10.0
           nTot      "==4"
           nNet      "==0"
         }
        .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14
           angle_to_track    10.0
           energyThreshold_b 0.025
           energyThreshold_e 0.050
         }
        .pid(method: :probability) {
           prob_cut 0.0
           identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.06,
                                           treat_as_electron_if_energy_above: 1.1
           identify :pion, against: [:kaon]
           npip "==1"
           npim "==1"
           nlp  "==1"
           nlm  "==1"
         }
        .kinematic_fit([:pip, :pim, :lp, :lm]) {
           nominal
           constrain_four_momentum
           invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
           chi2_cut 40
         }

alg_ch2
  .note(:lepton_classification,
        "Charged tracks with p < 1.06 GeV/c treated as pions; others as leptons. Muons: E_EMC < 0.35 GeV. Electrons: E_EMC > 1.1 GeV.")
  .note(:pipi_angle_veto,
        "cos(theta_pi+pi-) < 0.98 to remove radiative Bhabha / dimuon with photon conversion mimicking pi+pi-.")
  .note(:pi_e_angle_veto,
        "For J/psi -> e+e-: cos(theta_pi_e) < 0.98 to reject gamma-conversion background.")
  .note(:muon_hit_requirement,
        "For J/psi -> mu+ mu-: at least one muon candidate must have >5 MUC hit layers.")
  .note(:Jpsi_mass_window,
        "M(l+ l-) in (3.090, 3.105) GeV/c^2 for signal region; sidebands (3.030, 3.060) U (3.140, 3.170) GeV/c^2.")
  .with_decay_card(decay_card_pipiJpsi)
  .apply(sel_ch2)
alg_ch2.execute_on([data_4390, data_4420, incMC_4390, incMC_4420] + exMC_ch2)

######################################################################
### Channel 3: e+e- -> pi+ pi- h_c ; h_c -> gamma eta_c             ###
### eta_c reconstructed via 16 hadronic decay modes                  ###
######################################################################
decay_card_pipiHc = <<~DECAYCARD
    Decay psi(4260)
    1.000 pi+ pi- h_c                        PHSP;
    Enddecay

    Decay h_c
    1.000 gamma eta_c                        HELAMP 1.0 0.0 1.0 0.0;
    Enddecay

    Decay eta_c
    1.000 pi+ pi- pi+ pi-                    PHSP;
    Enddecay

    End
DECAYCARD

exMC_ch3 = DatasetManager.create_exclusive_mc_for([data_4390, data_4420]) do |config|
  config.sample_name    = "TccBar_pipiHc_signal"
  config.events          = 500000
  config.decay_card      = decay_card_pipiHc
  config.cross_section   = :default
end
exMC_ch3.each { |m| m.save_to_config(format: :yaml, file_path: "exMC_ch3_#{m.name}") }

alg_ch3 = Algorithm.new("TccPipiHc")
alg_ch3.set_header(["TccPipiHcAlg/TccPipiHc.h"])
       .set_constant({"ECMS" => [:double, 4.416]})

sel_ch3 = Selection.new
sel_ch3.select_track {
           cos_theta 0.93
           Vr        1.0
           Vz        10.0
           nTot      ">=4"     # channels have 4, 6, or 8 charged tracks
           nNet      "==0"
         }
        .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14
           angle_to_track    10.0
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           nGam              ">=1"
         }
        .pid(method: :probability) {
           prob_cut 0.0
           identify :kaon,   against: [:pion, :proton]
           identify :pion,   against: [:kaon, :proton]
           identify :proton, against: [:kaon, :pion]
           npip ">=1"
           npim ">=1"
         }
        .kinematic_fit([:pip, :pim, :pip, :pim, :gamma]) {
           nominal
           constrain_four_momentum
           chi2_cut 35
         }

alg_ch3
  .note(:eta_c_decay_modes,
        "eta_c reconstructed via 16 hadronic modes: pp-bar, 2(pi+pi-), 2(K+K-), pi+pi-K+K-, pi+pi-pp-bar, 3(pi+pi-), 2(pi+pi-)K+K-, KS0 K+/- pi-/+, KS0 K+/- pi-/+ pi+pi-, K+K-pi0, pp-bar pi0, K+K-eta, pi+pi-eta, 2(pi+pi-)eta, pi+pi-pi0pi0, 2(pi+pi-)pi0. The number of charged tracks per mode is 4, 6, or 8.")
  .note(:KS0_reconstruction,
        "For eta_c modes with K_S0: two charged tracks constrained to common vertex; |M(pi+pi-) - m(K_S0)| < 20 MeV/c^2; decay length > 2 sigma.")
  .note(:pi0_eta_mass_windows,
        "pi0 from gamma gamma: M(gamma gamma) in (0.110, 0.150) GeV/c^2. eta: M(gamma gamma) in (0.500, 0.570) GeV/c^2.")
  .note(:best_candidate_selection,
        "For multi-channel matches with same number of tracks, best candidate chosen by minimum chi^2 = chi^2_4C + chi^2_PID + chi^2_(pi0/eta).")
  .note(:hc_signal_window,
        "For final states with only charged tracks: RM(pi+pi- gamma) in (2.934, 3.034) GeV/c^2 and chi^2_4C < 35. For final states with pi0/eta: RM(pi+pi- gamma) in (2.939, 3.029) GeV/c^2 and chi^2_4C < 20. h_c signal window: RM(pi+pi-) in (3.515, 3.535) GeV/c^2; sidebands (3.475,3.495) U (3.555, 3.575) GeV/c^2.")
  .note(:eta_c_mass_selection,
        "Multiple pi+ pi- (from e+e-) and gamma (from h_c) combinations: select candidate with minimum |M(eta_c) - m(eta_c)|.")
  .with_decay_card(decay_card_pipiHc)
  .apply(sel_ch3)
alg_ch3.execute_on([data_4390, data_4420, incMC_4390, incMC_4420] + exMC_ch3)
