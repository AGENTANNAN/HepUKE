# BESIII analysis: Observation of Lambda_c+ -> Lambda a0(980)+ and evidence for Sigma(1380)+ in Lambda_c+ -> Lambda pi+ eta
# arXiv: 2407.12270v2
# Single-tag method at threshold: e+e- -> Lambda_c+ Lambda_c- at sqrt(s) = 4.600-4.843 GeV, 6.1 fb^-1
# PWA of Lambda_c+ -> Lambda pi+ eta with Lambda -> p pi-, eta -> gamma gamma / pi+ pi- pi0, pi0 -> gamma gamma

### Datasets (11 energy points) ###
data_4600 = DatasetManager.real_data.find("703_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")
data_4740 = DatasetManager.real_data.find("707_4740")
data_4750 = DatasetManager.real_data.find("707_4750")
data_4780 = DatasetManager.real_data.find("707_4780")
data_4840 = DatasetManager.real_data.find("707_4840")

all_energy_points = [data_4600, data_4610, data_4620, data_4640, data_4660,
                     data_4680, data_4700, data_4740, data_4750, data_4780,
                     data_4840]

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_4740 = DatasetManager.inclusive_mc.find("707_4740")
incMC_4750 = DatasetManager.inclusive_mc.find("707_4750")
incMC_4780 = DatasetManager.inclusive_mc.find("707_4780")
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")

all_incMC = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660,
             incMC_4680, incMC_4700, incMC_4740, incMC_4750, incMC_4780,
             incMC_4840]

### Decay cards ###
# Mode I: eta -> gamma gamma
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.000 Lambda pi+ eta PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.000 anti-p- pi+ anti-Lambda PHSP;
    Enddecay

    Decay anti-Lambda
    1.000 anti-p- pi+ PHSP;
    Enddecay

    Decay Lambda
    1.000 p+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: eta -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.000 Lambda pi+ eta PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.000 anti-p- pi+ anti-Lambda PHSP;
    Enddecay

    Decay anti-Lambda
    1.000 anti-p- pi+ PHSP;
    Enddecay

    Decay Lambda
    1.000 p+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC ###
exMC_modeI = DatasetManager.create_exclusive_mc_for(all_energy_points) do |config|
  config.sample_name   = "lcp_lapi_eta_gg"
  config.events        = 1_000_000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc_for(all_energy_points) do |config|
  config.sample_name   = "lcp_lapi_eta_pipipi0"
  config.events        = 1_000_000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

### Algorithm I: Lambda_c+ -> Lambda pi+ eta, eta -> gamma gamma ###
alg_modeI = Algorithm.new("LcLambdaPiEtaGG")
alg_modeI.set_header(["LcLambdaPiEtaGGAlg/LcLambdaPiEtaGG.h"])
  .set_constant({"ECMS" => [:double, 4.682]})
  .note(:single_tag_method, "Single-tag reconstruction of Lambda_c+ from Lambda_c+ Lambda_c- pair
    production at threshold. Lambda_c+ reconstructed from Lambda pi+ eta; recoiling anti-Lambda_c-
    inferred from beam constraint.")
  .note(:deltae_window, "DeltaE window: -0.1 < DeltaE < 0.1 GeV. Keep candidate with minimum |DeltaE|
    per event. DeltaE = E_tag - E_beam in e+e- rest frame.")
  .note(:mbc_fit, "Signal yield extracted from extended unbinned maximum likelihood fit to M_BC
    distribution. Components: signal MC shape, mismatched background, Lambda_c+ decay backgrounds,
    combinatorial background (ARGUS function).")
  .note(:bdtg, "BDTG applied post-selection with score > 0.95 for eta->gamma gamma channel.
    Input variables: DeltaE, M(p pi-), L/sigma_L, M(gamma gamma), cos_theta_eta,
    Lat(gamma_high), Lat(gamma_low). Trained on inclusive MC.")
  .note(:post_kinematic_fit, "Three-constraint kinematic fit applied post-BOSS-level: Lambda pi+ eta
    invariant mass constrained to Lambda_c+ mass, p pi- constrained to Lambda mass, recoil
    anti-Lambda_c- constrained to Lambda_c+ mass. Updated momenta used in PWA.")
  .note(:pwa, "Partial wave analysis performed using helicity amplitude formalism via TF-PWA.
    Components: Lambda a0(980)+, Sigma(1385)+ eta, Lambda(1670) pi+, and non-resonant S-wave.")

sel_modeI = Selection.new
sel_modeI.select_track {
            cos_theta 0.93
            Vz 10.0
            Vr 1.0
            nChrp ">=2"
            nChrn ">=1"
          }
         .select_photon {
            tdc_emc_start 0
            tdc_emc_end 14
            angle_to_track 10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam ">=2"
         }
         .pid(method: :probability) {
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]
            identify :pion, against: [:kaon]
            nprp ">=1"
            npip ">=1"
         }
         .remove([:prp <= :chrgp])
         .remove([:prm <= :chrgn])
         .secondary_vertex_fit([:prp, :pim]) {
            build_virtual_particle(:Lambda).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
         }
         .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
            chi2_cut 20
            neta ">=1"
         }
         .assign({:chrgp => :pip, :chrgn => :pim})
         .kinematic_fit([:Lambda, :pip, :eta, :pip]) {
            nominal
            constrain_four_momentum
            chi2_cut 200
         }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on(all_energy_points + all_incMC + exMC_modeI)

### Algorithm II: Lambda_c+ -> Lambda pi+ eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma ###
alg_modeII = Algorithm.new("LcLambdaPiEtaPiPiPi0")
alg_modeII.set_header(["LcLambdaPiEtaPiPiPi0Alg/LcLambdaPiEtaPiPiPi0.h"])
  .set_constant({"ECMS" => [:double, 4.682]})
  .note(:single_tag_method, "Single-tag reconstruction of Lambda_c+. Same method as Mode I.")
  .note(:deltae_window, "DeltaE window: -0.1 < DeltaE < 0.1 GeV. Best candidate per event
    selected by minimum |DeltaE|.")
  .note(:mbc_fit, "Signal yield from extended unbinned maximum likelihood fit to M_BC distribution.")
  .note(:bdtg, "BDTG score > 0.97 for eta -> pi+ pi- pi0 channel.
    Additional input variable: M(pi+ pi- pi0).")
  .note(:post_kinematic_fit, "3C kinematic fit: Lambda pi+ eta -> Lambda_c+ mass,
    p pi- -> Lambda mass, recoil anti-Lambda_c- -> Lambda_c+ mass.")
  .note(:pwa, "PWA via TF-PWA helicity amplitude formalism.")

sel_modeII = Selection.new
sel_modeII.select_track {
             cos_theta 0.93
             Vz 10.0
             Vr 1.0
             nChrp ">=2"
             nChrn ">=2"
           }
          .select_photon {
             tdc_emc_start 0
             tdc_emc_end 14
             angle_to_track 10.0
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             nGam ">=2"
          }
          .pid(method: :probability) {
             prob_cut 0.001
             identify :proton, against: [:kaon, :pion]
             identify :pion, against: [:kaon]
             nprp ">=1"
             npip ">=1"
             npim ">=1"
          }
          .remove([:prp <= :chrgp])
          .remove([:prm <= :chrgn])
          .secondary_vertex_fit([:prp, :pim]) {
             build_virtual_particle(:Lambda).by_minimizing_mass_difference
             remove_used_particle_from_candidate_list
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 200
             npi0 ">=1"
          }
          .assign({:chrgp => :pip, :chrgn => :pim})
          .kinematic_fit([:Lambda, :pip, :pip, :pim, :pi0]) {
             nominal
             constrain_four_momentum
             chi2_cut 200
          }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on(all_energy_points + all_incMC + exMC_modeII)