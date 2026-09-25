# =====================================================================
# BESIII @ 3.686 GeV: search for psi(3686) -> pi0 hc
#   Mode I   : hc -> p pbar pi+ pi- pi0                (11 pi's? no: ppbar pi+pi- pi0)
#   Mode IIa : hc -> p pbar eta, eta -> gamma gamma
#   Mode IIb : hc -> p pbar eta, eta -> pi+ pi- pi0
#   Mode III : hc -> p pbar pi0
# BOSS (dataset preparation + event selection) part only.
# =====================================================================

### -------------------- Datasets -------------------- ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(3686) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # matching inclusive MC

### -------------------- Decay cards (EvtGen) -------------------- ###

# Mode I : psi(2S) -> pi0 hc , hc -> p+ anti-p- pi+ pi- pi0
decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 hc PHSP;
    Enddecay
    Decay hc
    1.0000 p+ anti-p- pi+ pi- pi0 PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode IIa : psi(2S) -> pi0 hc , hc -> p+ anti-p- eta , eta -> gamma gamma
decay_card_modeIIa = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 hc PHSP;
    Enddecay
    Decay hc
    1.0000 p+ anti-p- eta PHSP;
    Enddecay
    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode IIb : psi(2S) -> pi0 hc , hc -> p+ anti-p- eta , eta -> pi+ pi- pi0
decay_card_modeIIb = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 hc PHSP;
    Enddecay
    Decay hc
    1.0000 p+ anti-p- eta PHSP;
    Enddecay
    Decay eta
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Mode III : psi(2S) -> pi0 hc , hc -> p+ anti-p- pi0
decay_card_modeIII = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 hc PHSP;
    Enddecay
    Decay hc
    1.0000 p+ anti-p- pi0 PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

### -------------------- Exclusive MC (500k events each) -------------------- ###
# All four samples share the psi(2S) -> pi0 hc production.
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_pi0hc_hc_pppipimpip0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end
exMC_modeI.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_modeIIa = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_pi0hc_hc_ppeta_etagammagamma"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_modeIIa
  config.cross_section   = :default
end
exMC_modeIIa.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_modeIIb = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_pi0hc_hc_ppeta_etapipimpip0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_modeIIb
  config.cross_section   = :default
end
exMC_modeIIb.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_modeIII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_pi0hc_hc_pppip0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_modeIII
  config.cross_section   = :default
end
exMC_modeIII.save_to_config(format: :yaml, file_path: 'temp_for_test')

### =====================================================================
### Mode I : hc -> p pbar pi+ pi- pi0   (fit: p pbar pi+ pi- pi0 pi0 pi0)
### =====================================================================
alg_I = Algorithm.new("HcToPpPipPimPi0")
alg_I.set_header(["HcToPpPipPimPi0Alg/HcToPpPipPimPi0.h"])
     .set_constant({"ECMS" => [:double, 3.686]})
     .set_alias({"std::vector<double>" => "Vdouble"})

sel_I = Selection.new
  .select_track {              # >=2 positive and >=2 negative tracks
    cos_theta 0.93
    Vz        10.0             # |Vz| < 10 cm
    Vr        1.0              # Vr   < 1 cm
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {             # >=6 photons (3 pi0)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=6"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:pion]   # p / pbar vs pions
    identify :pion,   against: [:kaon]   # pi+ / pi- vs kaons (where present)
    nprp ">=1"
    nprm ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma mass-constrained fit
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0     ">=3"
  }
  .kinematic_fit([:prp, :prm, :pip, :pim, :pi0, :pi0, :pi0]) {  # 6C fit
    nominal
    constrain_four_momentum
    chi2_cut 200                                              # loose; tight cut in ROOT
  }

alg_I.note(:chi2_6c_final_cut, "final chi2_6C < 45 applied at ROOT level (BOSS runs the loose chi2_cut 200)")
     .with_decay_card(decay_card_modeI).apply(sel_I)

alg_I.execute_on([psip_data, psip_incMC, exMC_modeI])

### =====================================================================
### Mode IIa : hc -> p pbar eta, eta -> gamma gamma   (fit: p pbar pi0 eta)
### =====================================================================
alg_IIa = Algorithm.new("HcToPpEtaGamGam")
alg_IIa.set_header(["HcToPpEtaGamGamAlg/HcToPpEtaGamGam.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_IIa = Selection.new
  .select_track {              # >=1 positive and >=1 negative track
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=1"
    nChrn     ">=1"
    nNet      "==0"
  }
  .select_photon {             # >=4 photons (eta -> gg + pi0 -> gg)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:pion]
    nprp ">=1"
    nprm ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0     ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # eta -> gamma gamma (chi2 < 200)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta     ">=1"
  }
  .kinematic_fit([:prp, :prm, :pi0, :eta]) {  # 6C fit
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_IIa.note(:chi2_6c_final_cut, "final chi2_6C < 45 for mode IIa applied at ROOT level (BOSS runs chi2_cut 200)")
       .with_decay_card(decay_card_modeIIa).apply(sel_IIa)

alg_IIa.execute_on([psip_data, psip_incMC, exMC_modeIIa])

### =====================================================================
### Mode IIb : hc -> p pbar eta, eta -> pi+ pi- pi0   (fit: p pbar pi0 eta)
### =====================================================================
alg_IIb = Algorithm.new("HcToPpEtaPipPimPi0")
alg_IIb.set_header(["HcToPpEtaPipPimPi0Alg/HcToPpEtaPipPimPi0.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_IIb = Selection.new
  .select_track {              # >=2 positive and >=2 negative tracks (p, pbar, pi+, pi-)
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {             # >=4 photons (2 pi0)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:pion]
    identify :pion,   against: [:kaon]
    nprp ">=1"
    nprm ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma (>=2 pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0     ">=2"
  }
  .kinematic_fit([:prp, :prm, :pip, :pim, :pi0, :pi0]) {   # 6C fit, eta built from pi+pi-pi0
    nominal
    constrain_four_momentum
    invariant_mass_of(:pip, :pim, :pi0).within(0.532, 0.562)  # eta -> pi+ pi- pi0 window
    chi2_cut 200
  }

alg_IIb.with_decay_card(decay_card_modeIIb).apply(sel_IIb)

alg_IIb.execute_on([psip_data, psip_incMC, exMC_modeIIb])

### =====================================================================
### Mode III : hc -> p pbar pi0   (fit: p pbar pi0 pi0)
### =====================================================================
alg_III = Algorithm.new("HcToPpPi0")
alg_III.set_header(["HcToPpPi0Alg/HcToPpPi0.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_III = Selection.new
  .select_track {              # >=1 positive and >=1 negative track
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=1"
    nChrn     ">=1"
    nNet      "==0"
  }
  .select_photon {             # >=4 photons (2 pi0)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:pion]
    nprp ">=1"
    nprm ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma (>=2 pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0     ">=2"
  }
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) {  # 6C fit
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_III.note(:chi2_6c_final_cut, "final chi2_6C < 64 for mode III applied at ROOT level (BOSS runs chi2_cut 200)")
       .with_decay_card(decay_card_modeIII).apply(sel_III)

alg_III.execute_on([psip_data, psip_incMC, exMC_modeIII])