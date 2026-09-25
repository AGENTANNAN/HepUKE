# 1902.03823v2: J/psi -> gamma eta' and eta' absolute BFs
# Data: J/psi (3.097 GeV), 1310.6 x 10^6 events
# Inclusive tag via photon conversion + 5 exclusive eta' decay modes

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_jpsi = DatasetManager.real_data.find("708_3097")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")

# ---- Inclusive: J/psi -> gamma eta' (photon conversion method) ----
# Photon conversion finder used; hard to express in standard DSL.
# Monitored via note; the exclusive modes below are the main analysis.

# ============================================================
# Mode I: eta' -> gamma pi+ pi-  (J/psi -> gamma gamma pi+ pi-)
# ============================================================
decay_card_I = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.000 gamma pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

exMC_I = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_JpsitoGammaEtap_ModeI"
  config.related_dataset = data_jpsi
  config.events = 100000
  config.decay_card = decay_card_I
  config.cross_section = :default
end

alg_I = Algorithm.new("EtapGammaPiPi")
alg_I.set_header(["EtapGammaPiPiAlg/EtapGammaPiPi.h"])
     .set_constant({"ECMS" => [:double, 3.097]})

sel_I = Selection.new
sel_I.select_track {
       cos_theta 0.93
       Vz 100.0
       Vr 10.0
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
       nGam ">=2"
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :pion, against: [:kaon]
       npip "==1"
       npim "==1"
     }
     .remove([:pip <= :chrgp, :pim <= :chrgn])
     .kinematic_fit([:gamma, :gamma, :pip, :pim]) {
       nominal
       constrain_four_momentum
       chi2_cut 200
     }
alg_I.with_decay_card(decay_card_I).apply(sel_I)
alg_I.execute_on([data_jpsi, incMC_jpsi, exMC_I])

# ============================================================
# Mode II: eta' -> eta pi+ pi-  (eta -> gamma gamma)
# ============================================================
decay_card_II = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.000 eta pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_II = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_JpsitoGammaEtap_ModeII"
  config.related_dataset = data_jpsi
  config.events = 100000
  config.decay_card = decay_card_II
  config.cross_section = :default
end

alg_II = Algorithm.new("EtapEtaPiPi")
alg_II.set_header(["EtapEtaPiPiAlg/EtapEtaPiPi.h"])
      .set_constant({"ECMS" => [:double, 3.097]})

sel_II = Selection.new
sel_II.select_track {
        cos_theta 0.93
        Vz 100.0
        Vr 10.0
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
        nGam ">=3"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon]
        npip "==1"
        npim "==1"
      }
      .remove([:pip <= :chrgp, :pim <= :chrgn])
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=1"
      }
      .kinematic_fit([:gamma, :eta, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg_II.with_decay_card(decay_card_II).apply(sel_II)
alg_II.execute_on([data_jpsi, incMC_jpsi, exMC_II])

# ============================================================
# Mode III: eta' -> eta pi0 pi0  (eta -> gamma gamma, pi0 -> gamma gamma)
# ============================================================
decay_card_III = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.000 eta pi0 pi0 PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_III = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_JpsitoGammaEtap_ModeIII"
  config.related_dataset = data_jpsi
  config.events = 100000
  config.decay_card = decay_card_III
  config.cross_section = :default
end

alg_III = Algorithm.new("EtapEtaPi0Pi0")
alg_III.set_header(["EtapEtaPi0Pi0Alg/EtapEtaPi0Pi0.h"])
       .set_constant({"ECMS" => [:double, 3.097]})

sel_III = Selection.new
sel_III.select_track {
         cos_theta 0.93
         Vz 100.0
         Vr 10.0
         nChrp "==0"
         nChrn "==0"
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
       .kalman_kinematic_fit([:gamma, :gamma]) {
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
         chi2_cut 25
         neta ">=1"
       }
       .kalman_kinematic_fit([:gamma, :gamma]) {
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
         chi2_cut 25
         npi0 ">=2"
       }
       .kinematic_fit([:gamma, :eta, :pi0, :pi0]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }
alg_III.with_decay_card(decay_card_III).apply(sel_III)
alg_III.execute_on([data_jpsi, incMC_jpsi, exMC_III])

# ============================================================
# Mode IV: eta' -> gamma omega  (omega -> pi+ pi- pi0, pi0 -> gamma gamma)
# ============================================================
decay_card_IV = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.000 gamma omega PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_IV = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_JpsitoGammaEtap_ModeIV"
  config.related_dataset = data_jpsi
  config.events = 100000
  config.decay_card = decay_card_IV
  config.cross_section = :default
end

alg_IV = Algorithm.new("EtapGammaOmega")
alg_IV.set_header(["EtapGammaOmegaAlg/EtapGammaOmega.h"])
      .set_constant({"ECMS" => [:double, 3.097]})

sel_IV = Selection.new
sel_IV.select_track {
        cos_theta 0.93
        Vz 100.0
        Vr 10.0
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
        nGam ">=3"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon]
        npip "==1"
        npim "==1"
      }
      .remove([:pip <= :chrgp, :pim <= :chrgn])
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      }
      .kinematic_fit([:gamma, :gamma, :pip, :pim, :pi0]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg_IV.with_decay_card(decay_card_IV).apply(sel_IV)
alg_IV.execute_on([data_jpsi, incMC_jpsi, exMC_IV])

# ============================================================
# Mode V: eta' -> gamma gamma  (J/psi -> gamma gamma gamma)
# ============================================================
decay_card_V = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_V = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_JpsitoGammaEtap_ModeV"
  config.related_dataset = data_jpsi
  config.events = 100000
  config.decay_card = decay_card_V
  config.cross_section = :default
end

alg_V = Algorithm.new("EtapGammaGamma")
alg_V.set_header(["EtapGammaGammaAlg/EtapGammaGamma.h"])
     .set_constant({"ECMS" => [:double, 3.097]})

sel_V = Selection.new
sel_V.select_track {
       cos_theta 0.93
       Vz 100.0
       Vr 10.0
       nChrp "==0"
       nChrn "==0"
       nNet "==0"
     }
     .select_photon {
       tdc_emc_start 0
       tdc_emc_end 14
       angle_to_track 10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam ">=3"
     }
     .kinematic_fit([:gamma, :gamma, :gamma]) {
       nominal
       constrain_four_momentum
       chi2_cut 200
     }
alg_V.with_decay_card(decay_card_V).apply(sel_V)
alg_V.execute_on([data_jpsi, incMC_jpsi, exMC_V])

# Photon conversion method for inclusive eta' tag not expressed in DSL;
# captured as notes for each algorithm because it provides the absolute BF
# normalization in ROOT analysis.
[alg_I, alg_II, alg_III, alg_IV, alg_V].each do |a|
  a.note(:photon_conversion_inclusive,
    "Inclusive eta' yield obtained from photon conversion method: " \
    "radiative photon from J/psi -> gamma eta' converts to e+e- pair; " \
    "photon conversion finder used; pi0 conversion veto applied; " \
    "BF(J/psi->gamma eta') = (5.27+/-0.03+/-0.05) x 10^-3 used for normalization")
  a.note(:radiative_photon_selection,
    "radiative photon selected as highest-energy photon in event (E~1.4 GeV at J/psi); " \
    "kinematic fit applied per mode; photon conversion efficiency correction factor " \
    "f = 1.0085+/-0.0050 applied in ROOT")
end