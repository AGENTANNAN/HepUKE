# ============================================================
# Paper: arXiv:2410.16912v1
# Measurement of Lambda_c+ → Lambda K_S0 K+, Lambda_c+ → Lambda K_S0 pi+, and Lambda_c+ → Lambda K*+
# Lambda → p pi-, K_S0 → pi+ pi-, K*+ → K_S0 pi+
# Data: 4.5 fb^-1 at 7 energy points (4.600-4.699 GeV)
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
### Decay Cards
### ============================================================

# Lambda_c+ → Lambda K_S0 K+, Lambda → p pi-, K_S0 → pi+ pi-
decay_Lc_LKsK = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay
  Decay Lambda_c+
  1.000 Lambda K_S0 K+ PHSP;
  Enddecay
  Decay Lambda
  1.000 p+ pi- PHSP;
  Enddecay
  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Lambda_c+ → Lambda K_S0 pi+, Lambda → p pi-, K_S0 → pi+ pi-
decay_Lc_LKspi = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay
  Decay Lambda_c+
  1.000 Lambda K_S0 pi+ PHSP;
  Enddecay
  Decay Lambda
  1.000 p+ pi- PHSP;
  Enddecay
  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Lambda_c+ → Lambda K*+, Lambda → p pi-, K*+ → K_S0 pi+, K_S0 → pi+ pi-
decay_Lc_LKstar = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay
  Decay Lambda_c+
  1.000 Lambda K*+ PHSP;
  Enddecay
  Decay Lambda
  1.000 p+ pi- PHSP;
  Enddecay
  Decay K*+
  1.000 K_S0 pi+ VSS;
  Enddecay
  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

### ============================================================
### Exclusive MC
### ============================================================

exMC_Lc_LKsK = DatasetManager.create_exclusive_mc_for(datasets_all) do |config|
  config.sample_name   = "sig_Lc_LKsK"
  config.events        = 100000
  config.decay_card    = decay_Lc_LKsK
  config.cross_section = :default
end

exMC_Lc_LKspi = DatasetManager.create_exclusive_mc_for(datasets_all) do |config|
  config.sample_name   = "sig_Lc_LKspi"
  config.events        = 100000
  config.decay_card    = decay_Lc_LKspi
  config.cross_section = :default
end

exMC_Lc_LKstar = DatasetManager.create_exclusive_mc_for(datasets_all) do |config|
  config.sample_name   = "sig_Lc_LKstar"
  config.events        = 100000
  config.decay_card    = decay_Lc_LKstar
  config.cross_section = :default
end

### ============================================================
### Algorithm 1: Lambda_c+ → Lambda K_S0 K+
### Final state: Lambda(p pi-) + K_S0(pi+ pi-) + K+
### 5 charged tracks; Lambda and K_S0 via secondary vertex fits
### Signal extraction: 2-D M_BC vs M(p pi-) fit (ROOT-level)
### ============================================================

alg_Lc_LKsK = Algorithm.new("LcToLKsK")
alg_Lc_LKsK.set_header(["LcToLKsKAlg/LcToLKsK.h"])
  .set_constant({"ECMS" => [:double, 4.600]})
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before secondary vertex and kinematic fits")
  .note(:lambda_selection, "Lambda candidates selected via secondary vertex fit (chi2 < 100) with decay vertex separated from IP by >= 2 sigma; mass window [1.090, 1.140] GeV/c^2 applied; proton PID L(p)>L(K) and L(p)>L(pi); pion daughter from Lambda has no PID requirement")
  .note(:ks_selection, "K_S0 candidates selected via secondary vertex fit (chi2 < 100) with decay vertex separated from IP by >= 2 sigma; mass window [0.450, 0.540] GeV/c^2 applied; K_S0 daughter pions have no PID requirement; tracks from K_S0/Lambda decays have Vz < 20 cm, no Vr requirement")
  .note(:best_candidate, "if multiple Lambda_c+ candidates per event, the one with smallest |DeltaE| (within [-0.02, 0.02] GeV) is retained")
  .note(:signal_yield_fit, "2-D extended unbinned maximum likelihood fit on M_BC and M(p pi-) distributions, simultaneous in K_S0 signal region [0.487,0.511] and sideband [0.450,0.470]U[0.520,0.540] GeV/c^2; fake-K_S0 background ratio f_Ks = 0.56; signal yield = 128.9 +/- 12.7; BF = (3.04 +/- 0.30 +/- 0.16) x 10^-3")

sel_Lc_LKsK = Selection.new
sel_Lc_LKsK.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=3"
  nChrn ">=2"
  nNet "==1"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam "<4"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :kaon, against: [:pion]
  nprp ">=1"
  nkp ">=1"
}
.remove([:prp <= :chrgp, :kp <= :chrgp])
.assign({chrgp: :pip, chrgn: :pim})
.secondary_vertex_fit([:prp, :pim]) {
  invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
  chi2_cut 100
  nLambda ">=1"
}
.secondary_vertex_fit([:pip, :pim]) {
  invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  chi2_cut 100
  nKs ">=1"
}
.kinematic_fit([:Lambda, :K_S0, :kp]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_Lc_LKsK.with_decay_card(decay_Lc_LKsK).apply(sel_Lc_LKsK)
alg_Lc_LKsK.execute_on(datasets_all + incMCs_all + exMC_Lc_LKsK)

### ============================================================
### Algorithm 2: Lambda_c+ → Lambda K_S0 pi+
### Final state: Lambda(p pi-) + K_S0(pi+ pi-) + pi+
### 5 charged tracks; same secondary vertex fit pattern as mode 1
### Signal extraction: 3-D M_BC, M(pi+pi-), M(K_S0 pi+) fits (ROOT)
### Also used for Lambda_c+ → Lambda K*+ measurement
### ============================================================

alg_Lc_LKspi = Algorithm.new("LcToLKspi")
alg_Lc_LKspi.set_header(["LcToLKspiAlg/LcToLKspi.h"])
  .set_constant({"ECMS" => [:double, 4.600]})
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before secondary vertex and kinematic fits")
  .note(:lambda_selection, "Lambda candidates selected via secondary vertex fit (chi2 < 100) with decay vertex separated from IP by >= 2 sigma; mass window [1.090, 1.140] GeV/c^2 applied; proton PID L(p)>L(K) and L(p)>L(pi); pion daughter from Lambda has no PID requirement")
  .note(:ks_selection, "K_S0 candidates selected via secondary vertex fit (chi2 < 100) with decay vertex separated from IP by >= 2 sigma; mass window [0.450, 0.540] GeV/c^2 applied; K_S0 daughter pions have no PID requirement")
  .note(:best_candidate, "if multiple Lambda_c+ candidates per event, the one with smallest |DeltaE| (within [-0.02, 0.02] GeV) is retained")
  .note(:signal_yield_fit, "3-D extended unbinned maximum likelihood fit on M_BC, M(pi+pi-) and M(K_S0 pi+) distributions, simultaneous in Lambda signal region [1.111,1.121] and sideband [1.090,1.100]U[1.130,1.140] GeV/c^2; fake-Lambda background ratio f_Lambda = 0.50; non-resonant and K*+ components fitted with interference via phase angle theta_0; Lambda_c+ → Lambda K_S0 pi+ observed at 8.9 sigma; BF = (1.73 +/- 0.26 +/- 0.10) x 10^-3")
  .note(:kstar_analysis, "Lambda_c+ → Lambda K*+ studied as intermediate resonance; K*+ → K_S0 pi+ reconstructed via M(K_S0 pi+) distribution; 3-D fits under different interference assumptions (theta_0 = 0° to 360°); evidence at 4.7 sigma; BF under three scenarios: (2.40 +/- 0.58 +/- 0.11)x10^-3 (no interference), (5.21 +/- 0.71 +/- 0.25)x10^-3 (theta_0=109°), (1.29 +/- 0.44 +/- 0.06)x10^-3 (theta_0=221°)")

sel_Lc_LKspi = Selection.new
sel_Lc_LKspi.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=3"
  nChrn ">=2"
  nNet "==1"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam "<4"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon]
  nprp ">=1"
}
.remove([:prp <= :chrgp])
.assign({chrgp: :pip, chrgn: :pim})
.secondary_vertex_fit([:prp, :pim]) {
  invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
  chi2_cut 100
  nLambda ">=1"
}
.secondary_vertex_fit([:pip, :pim]) {
  invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  chi2_cut 100
  nKs ">=1"
}
.kinematic_fit([:Lambda, :K_S0, :pip]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_Lc_LKspi.with_decay_card(decay_Lc_LKspi).apply(sel_Lc_LKspi)
alg_Lc_LKspi.execute_on(datasets_all + incMCs_all + exMC_Lc_LKspi)

### ============================================================
### Algorithm 3: Lambda_c+ → Lambda K*+ (K*+ → K_S0 pi+)
### Same final state and selection as mode 2; K*+ resonance in decay card
### K*+ component extracted via M(K_S0 pi+) fit in ROOT
### ============================================================

alg_Lc_LKstar = Algorithm.new("LcToLKstar")
alg_Lc_LKstar.set_header(["LcToLKstarAlg/LcToLKstar.h"])
  .set_constant({"ECMS" => [:double, 4.600]})
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before secondary vertex and kinematic fits")
  .note(:lambda_selection, "same Lambda selection as mode 2 (secondary vertex fit, chi2 < 100, flight significance > 2 sigma, mass window [1.090,1.140] GeV/c^2)")
  .note(:ks_selection, "same K_S0 selection as mode 2 (secondary vertex fit, chi2 < 100, flight significance > 2 sigma, mass window [0.450,0.540] GeV/c^2)")
  .note(:best_candidate, "if multiple Lambda_c+ candidates per event, the one with smallest |DeltaE| (within [-0.02, 0.02] GeV) is retained")
  .note(:kstar_yield_fit, "K*+ signal extracted from 3-D simultaneous M_BC/M(pi+pi-)/M(K_S0 pi+) fit; K*+ described by Breit-Wigner (mass 891.66 MeV, width 50.8 MeV from PDG) convolved with Gaussian resolution; interference with non-resonant Lambda_c+ → Lambda K_S0 pi+ modeled via phase angle theta_0; evidence at 4.7 sigma significance")

sel_Lc_LKstar = Selection.new
sel_Lc_LKstar.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=3"
  nChrn ">=2"
  nNet "==1"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam "<4"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon]
  nprp ">=1"
}
.remove([:prp <= :chrgp])
.assign({chrgp: :pip, chrgn: :pim})
.secondary_vertex_fit([:prp, :pim]) {
  invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
  chi2_cut 100
  nLambda ">=1"
}
.secondary_vertex_fit([:pip, :pim]) {
  invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  chi2_cut 100
  nKs ">=1"
}
.kinematic_fit([:Lambda, :K_S0, :pip]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_Lc_LKstar.with_decay_card(decay_Lc_LKstar).apply(sel_Lc_LKstar)
alg_Lc_LKstar.execute_on(datasets_all + incMCs_all + exMC_Lc_LKstar)