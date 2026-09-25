# ============================================================
# Paper: arXiv:2410.12620v1
# Search for e+e- → phi chi_c0 and phi eta_c2(1D) at 4.47-4.95 GeV
# Data: 6.7 fb^-1, 16 energy points, BOSS 703/706/707
# ============================================================

### ============================================================
### Datasets — phi chi_c0 (all 16 energy points, 4.47-4.95 GeV)
### ============================================================

# BOSS 703
data_4470 = DatasetManager.real_data.find("703_4470")
data_4530 = DatasetManager.real_data.find("703_4530")
data_4575 = DatasetManager.real_data.find("703_4575")
data_4600 = DatasetManager.real_data.find("703_4600")
# BOSS 706
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")
# BOSS 707
data_4740 = DatasetManager.real_data.find("707_4740")
data_4750 = DatasetManager.real_data.find("707_4750")
data_4780 = DatasetManager.real_data.find("707_4780")
data_4840 = DatasetManager.real_data.find("707_4840")
data_4914 = DatasetManager.real_data.find("707_4914")
data_4946 = DatasetManager.real_data.find("707_4946")

datasets_all = [
  data_4470, data_4530, data_4575, data_4600,
  data_4610, data_4620, data_4640, data_4660,
  data_4680, data_4700, data_4740, data_4750,
  data_4780, data_4840, data_4914, data_4946
]

# Inclusive MC at each energy point
incMC_4470 = DatasetManager.inclusive_mc.find("703_4470")
incMC_4530 = DatasetManager.inclusive_mc.find("703_4530")
incMC_4575 = DatasetManager.inclusive_mc.find("703_4575")
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
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

incMCs_all = [
  incMC_4470, incMC_4530, incMC_4575, incMC_4600,
  incMC_4610, incMC_4620, incMC_4640, incMC_4660,
  incMC_4680, incMC_4700, incMC_4740, incMC_4750,
  incMC_4780, incMC_4840, incMC_4914, incMC_4946
]

# eta_c2(1D) search uses only datasets with sqrt(s) >= 4.84 GeV
datasets_eta_c2 = [data_4840, data_4914, data_4946]
incMCs_eta_c2  = [incMC_4840, incMC_4914, incMC_4946]

### ============================================================
### Decay Cards
### ============================================================

# Channel 1: e+e- → phi chi_c0, phi → K+K-, chi_c0 → pi+pi-
decay_chi_c0_pipi = <<~DECAYCARD
  Decay psi(4260)
  1.000 phi chi_c0 PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay chi_c0
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Channel 2: e+e- → phi chi_c0, phi → K+K-, chi_c0 → pi+pi-pi0pi0
decay_chi_c0_pipipi0pi0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 phi chi_c0 PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay chi_c0
  1.000 pi+ pi- pi0 pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Channel 3: e+e- → phi chi_c0, phi → K+K-, chi_c0 → K+K-pi+pi-
decay_chi_c0_KKpipi = <<~DECAYCARD
  Decay psi(4260)
  1.000 phi chi_c0 PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay chi_c0
  1.000 K+ K- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Channel 4: e+e- → phi chi_c0, phi → K+K-, chi_c0 → 2(pi+pi-)
decay_chi_c0_4pi = <<~DECAYCARD
  Decay psi(4260)
  1.000 phi chi_c0 PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay chi_c0
  1.000 pi+ pi- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Channel 5: e+e- → phi chi_c0, phi → K+K-, chi_c0 → 3(pi+pi-)
decay_chi_c0_6pi = <<~DECAYCARD
  Decay psi(4260)
  1.000 phi chi_c0 PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay chi_c0
  1.000 pi+ pi- pi+ pi- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# eta_c2(1D) decay card (same final states as chi_c0 but with eta_c2)
decay_eta_c2_pipi = <<~DECAYCARD
  Decay psi(4260)
  1.000 phi eta_c2 PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay eta_c2
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

### ============================================================
### Exclusive MC — batch-create for energy scan
### ============================================================

exMC_chi_c0_pipi = DatasetManager.create_exclusive_mc_for(datasets_all) do |config|
  config.sample_name   = "sig_phi_chi_c0_pipi"
  config.events        = 50000
  config.decay_card    = decay_chi_c0_pipi
  config.cross_section = :default
end

exMC_chi_c0_pipipi0pi0 = DatasetManager.create_exclusive_mc_for(datasets_all) do |config|
  config.sample_name   = "sig_phi_chi_c0_pipipi0pi0"
  config.events        = 50000
  config.decay_card    = decay_chi_c0_pipipi0pi0
  config.cross_section = :default
end

exMC_chi_c0_KKpipi = DatasetManager.create_exclusive_mc_for(datasets_all) do |config|
  config.sample_name   = "sig_phi_chi_c0_KKpipi"
  config.events        = 50000
  config.decay_card    = decay_chi_c0_KKpipi
  config.cross_section = :default
end

exMC_chi_c0_4pi = DatasetManager.create_exclusive_mc_for(datasets_all) do |config|
  config.sample_name   = "sig_phi_chi_c0_4pi"
  config.events        = 50000
  config.decay_card    = decay_chi_c0_4pi
  config.cross_section = :default
end

exMC_chi_c0_6pi = DatasetManager.create_exclusive_mc_for(datasets_all) do |config|
  config.sample_name   = "sig_phi_chi_c0_6pi"
  config.events        = 50000
  config.decay_card    = decay_chi_c0_6pi
  config.cross_section = :default
end

exMC_eta_c2_pipi = DatasetManager.create_exclusive_mc_for(datasets_eta_c2) do |config|
  config.sample_name   = "sig_phi_eta_c2_pipi"
  config.events        = 50000
  config.decay_card    = decay_eta_c2_pipi
  config.cross_section = :default
end

### ============================================================
### Algorithm 1: phi chi_c0 → K+K- pi+pi- (FULL reconstruction)
### Channel 1: chi_c0 → pi+pi-
### Final state: K+ K- pi+ pi- (4 charged tracks, no photons)
### ============================================================

alg_chi_c0_pipi_full = Algorithm.new("PhiChiC0PiPiFull")
alg_chi_c0_pipi_full.set_header(["PhiChiC0PiPiFullAlg/PhiChiC0PiPiFull.h"])
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 4C kinematic fit; efficiency difference between with/without estimated by re-running BOSS selection")

sel_chi_c0_pipi_full = Selection.new
sel_chi_c0_pipi_full.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==2"
  nNet "==0"
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
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
  nkp "==1"
  nkm "==1"
  npip "==1"
  npim "==1"
}
.kinematic_fit([:kp, :km, :pip, :pim]) {
  nominal
  constrain_four_momentum
  chi2_cut 40
}

alg_chi_c0_pipi_full.with_decay_card(decay_chi_c0_pipi).apply(sel_chi_c0_pipi_full)
alg_chi_c0_pipi_full.execute_on(datasets_all + incMCs_all + exMC_chi_c0_pipi)

### ============================================================
### Algorithm 2: phi chi_c0 → K+K- pi+pi-pi0pi0 (FULL reconstruction)
### Channel 2: chi_c0 → pi+pi-pi0pi0
### Final state: K+ K- pi+ pi- pi0 pi0 (4 charged + 2 pi0)
### ============================================================

alg_chi_c0_pipipi0pi0_full = Algorithm.new("PhiChiC0PiPiPi0Pi0Full")
alg_chi_c0_pipipi0pi0_full.set_header(["PhiChiC0PiPiPi0Pi0FullAlg/PhiChiC0PiPiPi0Pi0Full.h"])
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 4C kinematic fit")

sel_chi_c0_pipipi0pi0_full = Selection.new
sel_chi_c0_pipipi0pi0_full.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==2"
  nNet "==0"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam "<8"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
  nkp "==1"
  nkm "==1"
  npip "==1"
  npim "==1"
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0 ">=2"
}
.kinematic_fit([:kp, :km, :pip, :pim, :pi0, :pi0]) {
  nominal
  constrain_four_momentum
  chi2_cut 40
}

alg_chi_c0_pipipi0pi0_full.with_decay_card(decay_chi_c0_pipipi0pi0).apply(sel_chi_c0_pipipi0pi0_full)
alg_chi_c0_pipipi0pi0_full.execute_on(datasets_all + incMCs_all + exMC_chi_c0_pipipi0pi0)

### ============================================================
### Algorithm 3: phi chi_c0 → K+K- K+K-pi+pi- (FULL reconstruction)
### Channel 3: chi_c0 → K+K-pi+pi-
### Final state: 4 kaons + 2 pions (6 charged tracks)
### Note: higher-momentum kaon pair assigned to chi_c0, lower to phi
### ============================================================

alg_chi_c0_KKpipi_full = Algorithm.new("PhiChiC0KKPiPiFull")
alg_chi_c0_KKpipi_full.set_header(["PhiChiC0KKPiPiFullAlg/PhiChiC0KKPiPiFull.h"])
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 4C kinematic fit")

sel_chi_c0_KKpipi_full = Selection.new
sel_chi_c0_KKpipi_full.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==3"
  nChrn "==3"
  nNet "==0"
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
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
  nkp "==2"
  nkm "==2"
  npip "==1"
  npim "==1"
}
.kinematic_fit([:kp, :kp, :km, :km, :pip, :pim]) {
  nominal
  constrain_four_momentum
  chi2_cut 40
}

alg_chi_c0_KKpipi_full.with_decay_card(decay_chi_c0_KKpipi).apply(sel_chi_c0_KKpipi_full)
alg_chi_c0_KKpipi_full.execute_on(datasets_all + incMCs_all + exMC_chi_c0_KKpipi)

### ============================================================
### Algorithm 4: phi chi_c0 → K+K- 2(pi+pi-) (FULL reconstruction)
### Channel 4: chi_c0 → 2(pi+pi-)
### Final state: 2 kaons + 4 pions (6 charged tracks)
### ============================================================

alg_chi_c0_4pi_full = Algorithm.new("PhiChiC04PiFull")
alg_chi_c0_4pi_full.set_header(["PhiChiC04PiFullAlg/PhiChiC04PiFull.h"])
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 4C kinematic fit")

sel_chi_c0_4pi_full = Selection.new
sel_chi_c0_4pi_full.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==3"
  nChrn "==3"
  nNet "==0"
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
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
  nkp "==1"
  nkm "==1"
  npip "==2"
  npim "==2"
}
.kinematic_fit([:kp, :km, :pip, :pip, :pim, :pim]) {
  nominal
  constrain_four_momentum
  chi2_cut 40
}

alg_chi_c0_4pi_full.with_decay_card(decay_chi_c0_4pi).apply(sel_chi_c0_4pi_full)
alg_chi_c0_4pi_full.execute_on(datasets_all + incMCs_all + exMC_chi_c0_4pi)

### ============================================================
### Algorithm 5: phi chi_c0 → K+K- 3(pi+pi-) (FULL reconstruction)
### Channel 5: chi_c0 → 3(pi+pi-)
### Final state: 2 kaons + 6 pions (8 charged tracks)
### ============================================================

alg_chi_c0_6pi_full = Algorithm.new("PhiChiC06PiFull")
alg_chi_c0_6pi_full.set_header(["PhiChiC06PiFullAlg/PhiChiC06PiFull.h"])
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 4C kinematic fit")

sel_chi_c0_6pi_full = Selection.new
sel_chi_c0_6pi_full.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==4"
  nChrn "==4"
  nNet "==0"
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
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
  nkp "==1"
  nkm "==1"
  npip "==3"
  npim "==3"
}
.kinematic_fit([:kp, :km, :pip, :pip, :pip, :pim, :pim, :pim]) {
  nominal
  constrain_four_momentum
  chi2_cut 40
}

alg_chi_c0_6pi_full.with_decay_card(decay_chi_c0_6pi).apply(sel_chi_c0_6pi_full)
alg_chi_c0_6pi_full.execute_on(datasets_all + incMCs_all + exMC_chi_c0_6pi)

### ============================================================
### Algorithm 6: phi chi_c0 → K+K- pi+pi- (PARTIAL reco, missing K±)
### Channel 1 with missing kaon: 1C kinematic fit
### Net charge ∓1 depending on which kaon charge is missing
### chi2_cut < 5 if missing K from phi, < 45 if from chi_c0
### ============================================================

# --- 6a: Missing K- (K+ pi+ pi- in final state, nNet = +1) ---
alg_chi_c0_pipi_missKm = Algorithm.new("PhiChiC0PiPiMissKm")
alg_chi_c0_pipi_missKm.set_header(["PhiChiC0PiPiMissKmAlg/PhiChiC0PiPiMissKm.h"])
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 1C kinematic fit")

sel_chi_c0_pipi_missKm = Selection.new
sel_chi_c0_pipi_missKm.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==1"
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
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
  nkp "==1"
  npip "==1"
  npim "==1"
}
.kinematic_fit([:kp, :pip, :pim]) {
  nominal
  miss_track_of :km
  constrain_four_momentum
  chi2_cut 45
}

alg_chi_c0_pipi_missKm.with_decay_card(decay_chi_c0_pipi).apply(sel_chi_c0_pipi_missKm)
alg_chi_c0_pipi_missKm.execute_on(datasets_all + incMCs_all + exMC_chi_c0_pipi)

# --- 6b: Missing K+ (K- pi+ pi- in final state, nNet = -1) ---
alg_chi_c0_pipi_missKp = Algorithm.new("PhiChiC0PiPiMissKp")
alg_chi_c0_pipi_missKp.set_header(["PhiChiC0PiPiMissKpAlg/PhiChiC0PiPiMissKp.h"])
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 1C kinematic fit")

sel_chi_c0_pipi_missKp = Selection.new
sel_chi_c0_pipi_missKp.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==1"
  nChrn "==2"
  nNet "==-1"
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
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
  nkm "==1"
  npip "==1"
  npim "==1"
}
.kinematic_fit([:km, :pip, :pim]) {
  nominal
  miss_track_of :kp
  constrain_four_momentum
  chi2_cut 45
}

alg_chi_c0_pipi_missKp.with_decay_card(decay_chi_c0_pipi).apply(sel_chi_c0_pipi_missKp)
alg_chi_c0_pipi_missKp.execute_on(datasets_all + incMCs_all + exMC_chi_c0_pipi)

### ============================================================
### Algorithm 7: phi chi_c0 → K+K- pi+pi-pi0pi0 (PARTIAL reco, missing pi0)
### Channel 2 with missing pi0: 1C kinematic fit
### chi2_cut < 6
### ============================================================

alg_chi_c0_pipipi0pi0_missPi0 = Algorithm.new("PhiChiC0PiPiPi0Pi0MissPi0")
alg_chi_c0_pipipi0pi0_missPi0.set_header(["PhiChiC0PiPiPi0Pi0MissPi0Alg/PhiChiC0PiPiPi0Pi0MissPi0.h"])
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 1C kinematic fit")

sel_chi_c0_pipipi0pi0_missPi0 = Selection.new
sel_chi_c0_pipipi0pi0_missPi0.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==2"
  nNet "==0"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam "<8"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
  nkp "==1"
  nkm "==1"
  npip "==1"
  npim "==1"
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0 ">=1"
}
.kinematic_fit([:kp, :km, :pip, :pim, :pi0]) {
  nominal
  miss_track_of :pi0
  constrain_four_momentum
  chi2_cut 6
}

alg_chi_c0_pipipi0pi0_missPi0.with_decay_card(decay_chi_c0_pipipi0pi0).apply(sel_chi_c0_pipipi0pi0_missPi0)
alg_chi_c0_pipipi0pi0_missPi0.execute_on(datasets_all + incMCs_all + exMC_chi_c0_pipipi0pi0)

### ============================================================
### Algorithm 8: phi eta_c2(1D) → K+K- pi+pi- (FULL reconstruction)
### eta_c2(1D) search — only for sqrt(s) >= 4.84 GeV
### Same final state as channel 1 but chi2_cut < 60 (optimized for eta_c2)
### ============================================================

alg_eta_c2_pipi_full = Algorithm.new("PhiEtaC2PiPiFull")
alg_eta_c2_pipi_full.set_header(["PhiEtaC2PiPiFullAlg/PhiEtaC2PiPiFull.h"])
  .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 4C kinematic fit")

sel_eta_c2_pipi_full = Selection.new
sel_eta_c2_pipi_full.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==2"
  nNet "==0"
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
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
  nkp "==1"
  nkm "==1"
  npip "==1"
  npim "==1"
}
.kinematic_fit([:kp, :km, :pip, :pim]) {
  nominal
  constrain_four_momentum
  chi2_cut 60
}

alg_eta_c2_pipi_full.with_decay_card(decay_eta_c2_pipi).apply(sel_eta_c2_pipi_full)
alg_eta_c2_pipi_full.execute_on(datasets_eta_c2 + incMCs_eta_c2 + exMC_eta_c2_pipi)