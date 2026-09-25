# ============================================================
# Paper: arXiv:2410.13368v1
# First observation of Lambda_c+ → p pi0 (SCS decay)
# DNN-based single-tag strategy, 5.4 sigma significance
# Also measures Lambda_c+ → p eta as reference channel
# Data: 4.5 fb^-1 at 7 energy points (4.600-4.699 GeV)
# ============================================================

### ============================================================
### Datasets — 7 energy points
### ============================================================

data_4600 = DatasetManager.real_data.find("703_4600")   # 4.5995 GeV
data_4610 = DatasetManager.real_data.find("706_4610")   # 4.6119 GeV
data_4620 = DatasetManager.real_data.find("706_4620")   # 4.6280 GeV
data_4640 = DatasetManager.real_data.find("706_4640")   # 4.6409 GeV
data_4660 = DatasetManager.real_data.find("706_4660")   # 4.6612 GeV
data_4680 = DatasetManager.real_data.find("706_4680")   # 4.6819 GeV
data_4700 = DatasetManager.real_data.find("706_4700")   # 4.6988 GeV

datasets_all = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")

incMCs_all = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

### ============================================================
### Decay Cards
### ============================================================

# Lambda_c+ → p pi0, pi0 → gamma gamma; anti-Lambda_c- decays inclusively
decay_Lc_ppi0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay
  Decay Lambda_c+
  1.000 p+ pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Lambda_c+ → p eta, eta → gamma gamma; anti-Lambda_c- decays inclusively
decay_Lc_peta = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay
  Decay Lambda_c+
  1.000 p+ eta PHSP;
  Enddecay
  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

### ============================================================
### Exclusive MC
### ============================================================

exMC_ppi0 = DatasetManager.create_exclusive_mc_for(datasets_all) do |config|
  config.sample_name   = "sig_Lc_ppi0"
  config.events        = 50000
  config.decay_card    = decay_Lc_ppi0
  config.cross_section = :default
end

exMC_peta = DatasetManager.create_exclusive_mc_for(datasets_all) do |config|
  config.sample_name   = "sig_Lc_peta"
  config.events        = 50000
  config.decay_card    = decay_Lc_peta
  config.cross_section = :default
end

### ============================================================
### Algorithm 1: Lambda_c+ → p pi0
### Final state: proton + pi0 → proton + gamma gamma
### ============================================================

alg_Lc_ppi0 = Algorithm.new("LcP2pPi0")
alg_Lc_ppi0.set_header(["LcP2pPi0Alg/LcP2pPi0.h"])
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before the kinematic fit")
  .note(:dnn_selection, "Particle Transformer (ParT) DNN with iterative weighting applied to classify signal vs background; DNN output score > 0.95 required; DNN uses all charged tracks and isolated showers in the event not associated with the Lambda_c+ candidate; 20-model ensemble; model uncertainty estimated via ParticleNet alternative architecture; domain shift uncertainty estimated via control samples")

sel_Lc_ppi0 = Selection.new
sel_Lc_ppi0.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=1"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=2"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprp ">=1"
}
.remove([:prp <= :chrgp])
.assign({chrgp: :pip, chrgn: :pim})
.select_isolated_photon {
  angle_to_prm_track 20.0
  nGam ">=2"
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0 ">=1"
}
.kinematic_fit([:prp, :pi0]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_Lc_ppi0.with_decay_card(decay_Lc_ppi0).apply(sel_Lc_ppi0)
alg_Lc_ppi0.execute_on(datasets_all + incMCs_all + exMC_ppi0)

### ============================================================
### Algorithm 2: Lambda_c+ → p eta
### Final state: proton + eta → proton + gamma gamma
### ============================================================

alg_Lc_peta = Algorithm.new("LcP2pEta")
alg_Lc_peta.set_header(["LcP2pEtaAlg/LcP2pEta.h"])
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before the kinematic fit")
  .note(:dnn_selection, "Particle Transformer (ParT) DNN with iterative weighting applied to classify signal vs background; DNN output score > 0.95 required; DNN uses all charged tracks and isolated showers in the event not associated with the Lambda_c+ candidate; 20-model ensemble; model uncertainty estimated via ParticleNet alternative architecture; domain shift uncertainty estimated via control samples; combined DNN training with Lambda_c+ → p pi0")

sel_Lc_peta = Selection.new
sel_Lc_peta.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=1"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=2"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprp ">=1"
}
.remove([:prp <= :chrgp])
.assign({chrgp: :pip, chrgn: :pim})
.select_isolated_photon {
  angle_to_prm_track 20.0
  nGam ">=2"
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 200
  neta ">=1"
}
.kinematic_fit([:prp, :eta]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_Lc_peta.with_decay_card(decay_Lc_peta).apply(sel_Lc_peta)
alg_Lc_peta.execute_on(datasets_all + incMCs_all + exMC_peta)