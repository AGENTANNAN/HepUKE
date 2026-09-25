# ============================================================
# Paper: arXiv:2410.13515v2
# First observation of Lambda_c+ → n e+ nu_e (Cabibbo-suppressed beta decay)
# Double-tag technique: anti-Lambda_c- tagged via 10 hadronic modes,
# signal Lambda_c+ → neutron + positron + neutrino (semileptonic ST + missing)
# GNN-based (ParticleNet) signal/background classification
# Data: 4.5 fb^-1 at 7 energy points (4.600-4.699 GeV)
# Significance: >10 sigma
# ============================================================

### ============================================================
### Datasets — 7 energy points
### ============================================================

data_4600 = DatasetManager.real_data.find("703_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")

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
### Decay Card — e+e- → Lambda_c+ anti-Lambda_c-
### Signal: Lambda_c+ → n e+ nu_e
### Tag: anti-Lambda_c- → hadronic modes (reconstructed by DTagAlg)
### ============================================================

decay_semilep = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay
  Decay Lambda_c+
  1.000 n e+ nu_e PHSP;
  Enddecay
  End
DECAYCARD

### ============================================================
### Exclusive MC — batch for energy scan
### ============================================================

exMC_semilep = DatasetManager.create_exclusive_mc_for(datasets_all) do |config|
  config.sample_name   = "sig_Lc_ne_nu"
  config.events        = 100000
  config.decay_card    = decay_semilep
  config.cross_section = :default
end

### ============================================================
### Tag Analysis: anti-Lambda_c- tag (10 hadronic modes) + signal Lambda_c+ → n e+ nu_e
### ST + missing pattern (semileptonic)
### Signal side: 1 positron + EMC showers (neutron) + missing nu_e
### GNN classifies signal vs background using EMC shower patterns
### ============================================================

alg_Lc_ne_nu = TagAnalysis.new("LcNeNu")
alg_Lc_ne_nu.set_header(["LcNeNuAlg/LcNeNu.h"])
  .with_decay_card(decay_semilep)
  .note(:gnn_selection, "Graph Neural Network (ParticleNet) used for binary classification of signal vs background using EMC shower patterns; GNN trained on J/psi neutron and Lambda control samples with data-driven calibration; 100-model ensemble; domain shift uncertainty estimated via J/psi → Sigma+(→n pi+) anti-Sigma-(→anti-p pi0) and J/psi → Xi+(→Lambda pi+) anti-Xi-(→anti-Lambda pi-) control samples; iterative weighting method for mass sculpting mitigation")
  .note(:neutron_detection, "neutron detected via EMC shower patterns only; neutron cannot be reconstructed in MDC; EMC showers from neutron interactions analyzed by GNN; neutron control samples from J/psi → anti-p pi+ n used for GNN training and calibration")
  .note(:positron_pid, "positron ID: likelihood L(e) > 0.001 and L(e)/(L(e)+L(pi)+L(K)) > 0.8; EMC deposited energy / MDC momentum > 0.5")
  .note(:dt_yield_fit, "simultaneous binned maximum-likelihood fits to GNN output distributions; Lambda_c+ → n e+ nu_e and charge conjugate fitted separately; signal BF = 0.357% ± 0.034(stat) ± 0.014(syst); |V_cd| = 0.208 ± 0.011(exp) ± 0.007(LQCD) ± 0.001(tau)")

# Tag side: anti-Lambda_c- with 10 hadronic decay modes
alg_Lc_ne_nu.tag_side(:Lambdac) do |t|
  t.mode_group :hadronic
  t.charm -1
end

# Signal side: Lambda_c+ → n e+ nu_e
# neutron → EMC showers (at least 1), positron (1), neutrino (missing, massless)
alg_Lc_ne_nu.signal_side do |s|
  s.charged(ep: 1)
  s.photons 1
  s.require_charge 1
  s.missing :nu_e
end

# Kinematic fit: tag + e+ + nu = ecms_lab (4C)
alg_Lc_ne_nu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_Lc_ne_nu.apply
alg_Lc_ne_nu.execute_on(datasets_all + incMCs_all + exMC_semilep)