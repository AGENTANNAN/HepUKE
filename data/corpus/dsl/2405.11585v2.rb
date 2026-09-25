# BESIII DSL: 2405.11585v2 — Improved measurement of h_c -> gamma eta'/eta and search for h_c -> gamma pi0
# psi(3686) decays, (27.12e8) psi(3686) events, BOSS 709

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# =============================================
# Shared decay card
# =============================================
common_card = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# ===========================================================================
# Mode I: h_c -> gamma eta', eta' -> pi+ pi- eta, eta -> gamma gamma
# Final state: pi0 (6C: pi0 mass + eta mass), 5 photons, 2 charged
# ===========================================================================
decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- eta DIY;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_gamma_etap_pipimeta_gg"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_modeI
  config.cross_section = :default
end

alg_modeI = Algorithm.new("HcGamEtapPipimEta")
alg_modeI.set_header(["HcGamEtapPipimEtaAlg/HcGamEtapPipimEta.h"])
         .set_constant({ "ECMS" => [:double, 3.686] })

sel_modeI = Selection.new
sel_modeI.select_track {
           cos_theta 0.93
           Vz 10.0
           Vr 1.0
           nChrp "==1"
           nChrn "==1"
           nNet "==0"
         }
         .select_photon {
           tdc_emc_start 0
           tdc_emc_end 14
           angle_to_track 10.0
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           nGam ">=5"
         }
         .pid(method: :probability) {
           prob_cut 0.001
           identify :pion, against: [:kaon, :proton]
           npip "==1"
           npim "==1"
         }
         .kalman_kinematic_fit([:gamma, :gamma]) {
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 25
           npi0 ">=2"
         }
         .kalman_kinematic_fit([:gamma, :gamma]) {
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
           chi2_cut 25
           neta ">=1"
         }
         # 6C: 4C + pi0 mass (from psi(3686)) + eta mass (from eta'->pip pim eta chain)
         .kinematic_fit([:pip, :pim, :gamma, :pi0, :eta]) {
           nominal
           constrain_four_momentum
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
           chi2_cut 70
         }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on([psip_data, psip_incMC, exMC_modeI])

# ===========================================================================
# Mode II: h_c -> gamma eta', eta' -> gamma pi+ pi-
# Final state: pi0 + gamma + gamma + pip + pim (4 photons + 2 charged), 5C (4C + pi0 mass)
# ===========================================================================
decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.000 gamma pi+ pi- DIY;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_gamma_etap_gampipim"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_modeII
  config.cross_section = :default
end

alg_modeII = Algorithm.new("HcGamEtapGamPipim")
alg_modeII.set_header(["HcGamEtapGamPipimAlg/HcGamEtapGamPipim.h"])
          .set_constant({ "ECMS" => [:double, 3.686] })

sel_modeII = Selection.new
sel_modeII.select_track {
            cos_theta 0.93
            Vz 10.0
            Vr 1.0
            nChrp "==1"
            nChrn "==1"
            nNet "==0"
          }
          .select_photon {
            tdc_emc_start 0
            tdc_emc_end 14
            angle_to_track 10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam ">=4"
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :pion, against: [:kaon, :proton]
            npip "==1"
            npim "==1"
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"
          }
          # 5C: 4C + pi0 mass (from psi(3686) decay)
          .kinematic_fit([:pip, :pim, :gamma, :gamma, :pi0]) {
            nominal
            constrain_four_momentum
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 40
          }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on([psip_data, psip_incMC, exMC_modeII])

# ===========================================================================
# Mode III: h_c -> gamma eta, eta -> gamma gamma  (full neutral)
# Final state: pi0 + gamma + gamma + gamma (5 photons, no charged), 5C (4C + pi0 mass)
# ===========================================================================
decay_card_modeIII = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 gamma eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_modeIII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_gamma_eta_gg"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_modeIII
  config.cross_section = :default
end

alg_modeIII = Algorithm.new("HcGamEtaGG")
alg_modeIII.set_header(["HcGamEtaGGAlg/HcGamEtaGG.h"])
           .set_constant({ "ECMS" => [:double, 3.686] })
           .note(:full_neutral_timing, "If no charged particle in final state, require |T - T_max| <= 500 ns for EMC timing; no DSL method for this.")

sel_modeIII = Selection.new
sel_modeIII.select_photon {
             tdc_emc_start 0
             tdc_emc_end 14
             angle_to_track 10.0
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             nGam ">=5"
           }
           .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 25
             npi0 ">=1"
           }
           # 5C: 4C + pi0 mass (from psi(3686))
           .kinematic_fit([:gamma, :gamma, :gamma, :pi0]) {
             nominal
             constrain_four_momentum
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 40
           }

alg_modeIII.with_decay_card(decay_card_modeIII).apply(sel_modeIII)
alg_modeIII.execute_on([psip_data, psip_incMC, exMC_modeIII])

# ===========================================================================
# Mode IV: h_c -> gamma eta, eta -> pi+ pi- pi0
# Final state: 5 photons + pip + pim, 6C (4C + two pi0 masses), chi2 < 40
# Competing hypothesis veto: chi2_6C(gamma pip pim pi0 pi0) < chi2_6C(gamma pip pim eta pi0)
# ===========================================================================
decay_card_modeIV = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
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

exMC_modeIV = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_gamma_eta_pipimpi0"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_modeIV
  config.cross_section = :default
end

alg_modeIV = Algorithm.new("HcGamEtaPipimPi0")
alg_modeIV.set_header(["HcGamEtaPipimPi0Alg/HcGamEtaPipimPi0.h"])
          .set_constant({ "ECMS" => [:double, 3.686] })

sel_modeIV = Selection.new
sel_modeIV.select_track {
            cos_theta 0.93
            Vz 10.0
            Vr 1.0
            nChrp "==1"
            nChrn "==1"
            nNet "==0"
          }
          .select_photon {
            tdc_emc_start 0
            tdc_emc_end 14
            angle_to_track 10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam ">=5"
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify :pion, against: [:kaon, :proton]
            npip "==1"
            npim "==1"
          }
          .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=2"
          }
          # Nominal 6C: 4C + two pi0 mass constraints
          .kinematic_fit([:pip, :pim, :gamma, :pi0, :pi0]) {
            nominal
            constrain_four_momentum
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 40
          }
          # Competing hypothesis: gamma pip pim pi0 pi0 (no eta intermediate)
          # Used for chi2 veto: chi2_sig < chi2_bkg
          .kinematic_fit([:pip, :pim, :gamma, :pi0, :pi0]) {
            constrain_four_momentum
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          }

alg_modeIV.with_decay_card(decay_card_modeIV).apply(sel_modeIV)
alg_modeIV.execute_on([psip_data, psip_incMC, exMC_modeIV])

# ===========================================================================
# Mode V: h_c -> gamma pi0, pi0 -> gamma gamma  (full neutral, search)
# Final state: 5 photons, 6C (4C + two pi0 mass constraints), chi2 < 20
# ===========================================================================
decay_card_modeV = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.000 gamma pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_modeV = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_gamma_pi0"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_modeV
  config.cross_section = :default
end

alg_modeV = Algorithm.new("HcGamPi0")
alg_modeV.set_header(["HcGamPi0Alg/HcGamPi0.h"])
         .set_constant({ "ECMS" => [:double, 3.686] })
         .note(:full_neutral_timing, "If no charged particle, require |T - T_max| <= 500 ns for EMC timing.")

sel_modeV = Selection.new
sel_modeV.select_photon {
           tdc_emc_start 0
           tdc_emc_end 14
           angle_to_track 10.0
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           nGam ">=5"
         }
         .kalman_kinematic_fit([:gamma, :gamma]) {
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 25
           npi0 ">=2"
         }
         # 6C: 4C + two pi0 mass constraints
         .kinematic_fit([:gamma, :pi0, :pi0]) {
           nominal
           constrain_four_momentum
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 20
         }

alg_modeV.with_decay_card(decay_card_modeV).apply(sel_modeV)
alg_modeV.execute_on([psip_data, psip_incMC, exMC_modeV])