# Paper: 2504.04096v3 — Observation of three resonant structures in e+e- -> pi+ pi- hc
# Cross section scan at 59 energy points (4.009-4.950 GeV)
# hc -> gamma eta_c, eta_c -> 16 hadronic final states
# Uses XYZ-I, XYZ-II, R-scan datasets

### Dataset preparation ###
# Energy scan: 59 points from 4.009 to 4.950 GeV across XYZ-I, XYZ-II, R-scan
# Common samples accessible via DatasetManager.real_data; exact sample names from BES3_dataset.md

# Representative energy points for XYZ-I/XYZ-II/R-scan
# XYZ-I: 703_4009, 703_4090, 703_4100, ... (multiple points)
# XYZ-II: 703_4200, 703_4210, ... (multiple points)
# R-scan: 713_Rscan_2125, ... (multiple points)
# Full 59-point enumeration deferred to execute_on via dataset table lookup

# Load all scan data points (representative structure)
data_psi4260 = DatasetManager.real_data.find("703_4260")  # representative
incMC_psi4260 = DatasetManager.inclusive_mc.find("703_4260")

# Note: full 59-point energy scan requires enumerating all energy points from
# BES3_dataset.md. The pattern for multi-point energy scans uses
# create_exclusive_mc_for(data_points) or individual create_exclusive_mc calls.

### Mode classification: 16 eta_c decay modes ###
# Each mode has different track/photon content; Rule T1 applies
# Grouped by final state composition for DSL expressibility:

# --- MODE GROUP A: 4-charged-track modes ---
# Mode 1a: K_S0 K+ pi- (2 charged tracks + K_S0 -> 2 more)
# Mode 1b: K_S0 K- pi+ (CC of 1a)
# Mode 4: K+ K- pi+ pi-
# Mode 14: 3(pi+ pi-)
# ... (other 4-track modes)

# --- Example: Mode 4: K+ K- pi+ pi- (4 charged tracks only, no photon) ---
decay_card_M4 = <<~DECAYCARD
    Decay psi(4260)
    1.000 pi+ pi- hc PHSP;
    Enddecay

    Decay hc
    1.000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.000 K+ K- pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_M4 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_pipihc_etac_KKpipi"
  config.related_dataset = data_psi4260
  config.events = 200_000
  config.decay_card = decay_card_M4
  config.cross_section = :default
end

alg_M4 = Algorithm.new("HcEtaC_KKpipi")
alg_M4.set_header(["HcEtaC_KKpipiAlg/HcEtaC_KKpipi.h"])
       .set_constant({"ECMS" => [:double, 4.260]})

sel_M4 = Selection.new
sel_M4
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"     # pi+ + K+
    nChrn ">=2"     # pi- + K-
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=1"      # gamma from hc decay
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"
    nkm ">=1"
    npip ">=1"
    npim ">=1"
  }
  .kinematic_fit([:pip, :pim, :kp, :km, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_M4.with_decay_card(decay_card_M4).apply(sel_M4)

# --- Example: Mode 5: pi+ pi- eta, eta -> gamma gamma (4 charged + 2 photon from eta) ---
decay_card_M5 = <<~DECAYCARD
    Decay psi(4260)
    1.000 pi+ pi- hc PHSP;
    Enddecay

    Decay hc
    1.000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_M5 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_pipihc_etac_pipieta_gg"
  config.related_dataset = data_psi4260
  config.events = 200_000
  config.decay_card = decay_card_M5
  config.cross_section = :default
end

alg_M5 = Algorithm.new("HcEtaC_pipieta_gg")
alg_M5.set_header(["HcEtaC_pipieta_ggAlg/HcEtaC_pipieta_gg.h"])
       .set_constant({"ECMS" => [:double, 4.260]})

sel_M5 = Selection.new
sel_M5
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=2"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=3"      # 1 (hc) + 2 (eta->gamma gamma)
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=2"
    npim ">=2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:pip, :pim, :pip, :pim, :gamma, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_M5.with_decay_card(decay_card_M5).apply(sel_M5)

# --- Example: Mode 13: 2(pi+ pi-) pi0 (eta_c -> 2(pi+ pi-) pi0) ---
decay_card_M13 = <<~DECAYCARD
    Decay psi(4260)
    1.000 pi+ pi- hc PHSP;
    Enddecay

    Decay hc
    1.000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.000 pi+ pi- pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_M13 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_pipihc_etac_4pi_pi0"
  config.related_dataset = data_psi4260
  config.events = 200_000
  config.decay_card = decay_card_M13
  config.cross_section = :default
end

alg_M13 = Algorithm.new("HcEtaC_4pi_pi0")
alg_M13.set_header(["HcEtaC_4pi_pi0Alg/HcEtaC_4pi_pi0.h"])
        .set_constant({"ECMS" => [:double, 4.260]})

sel_M13 = Selection.new
sel_M13
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=3"
    nChrn ">=3"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=3"      # 1 (hc) + 2 (pi0->gamma gamma)
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=3"
    npim ">=3"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pip, :pim, :pip, :pim, :pip, :pim, :gamma, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_M13.with_decay_card(decay_card_M13).apply(sel_M13)


# === NOTES: Inexpressible portions ===
# 59 energy points across XYZ-I, XYZ-II, R-scan require full BES3_dataset.md enumeration
alg_M4.note(:energy_scan,
  "59 energy points from 4.009-4.950 GeV across XYZ-I, XYZ-II, R-scan datasets. Full scan: load all data points via DatasetManager.real_data.find() with exact sample names from BES3_dataset.md; create_exclusive_mc_for(scan_points) for each mode.")

# Remaining 13 eta_c decay modes (beyond the 3 shown as examples):
# Mode 1a: K_S0 K+ pi- (requires secondary_vertex_fit for K_S0)
# Mode 1b: K_S0 K- pi+ (CC)
# Mode 2: K+ K- pi0 (requires kalman_kinematic_fit for pi0)
# Mode 3: K+ K- pi+ pi- pi0 (similar to M13)
# Mode 6: pi+ pi- eta, eta->pi+pi-pi0 (additional pi0 from eta)
# Mode 7: pi+ pi- eta', eta'->pi+pi-eta, eta->gamma gamma
# Mode 8: pi+ pi- eta', eta'->gamma rho0, rho0->pi+pi-
# Mode 9: K+ K- eta, eta->gamma gamma
# Mode 10: K+ K- eta, eta->pi+pi-pi0
# Mode 11: p anti-p
# Mode 12: 2(pi+ pi-)
# Mode 14: 3(pi+ pi-)
# Mode 15: K+ K- 2(pi+ pi-)
# Mode 16: K+ K- 3(pi+ pi-)
alg_M4.note(:additional_modes,
  "13 additional eta_c decay modes (K_S0 K pi, KK pi0, KK 3pi pi0, pipi eta (pi+pi-pi0), pipi etap (pipi eta), pipi etap (gamma rho0), KK eta (gamma gamma), KK eta (pi+pi-pi0), p pbar, 4pi, 6pi, KK 4pi, KK 6pi) each require separate Algorithm+Selection per Rule T1")

# Cross section measurement specifics: Born cross section, ISR correction,
# vacuum polarization, efficiency maps per energy point — all ROOT-level
alg_M4.note(:cross_section,
  "Born cross section measurement: ISR correction factor (1+delta), vacuum polarization factor, and efficiency per energy point computed in ROOT; three resonant structures (Y(4220), Y(4390), Y(4660)) fitted in ROOT analysis")

# Execute on representative dataset only; full 59-point scan deferred
alg_M4.execute_on([data_psi4260, incMC_psi4260, exMC_M4])
alg_M5.execute_on([data_psi4260, incMC_psi4260, exMC_M5])
alg_M13.execute_on([data_psi4260, incMC_psi4260, exMC_M13])