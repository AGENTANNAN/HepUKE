# =============================================================================
# BOSS / HepScript DSL
#   e+e- -> pi+ pi- h_c ,  h_c -> gamma eta_c
#   eta_c reconstructed in 16 exclusive hadronic final states
#   (observation of Zc(4020), search for Zc(3900))
# Scope: dataset preparation + event selection up to the nominal 4C fit.
# =============================================================================

### ------------------------------- Datasets ------------------------------- ###
# 13 c.m. energies from 3.900 to 4.420 GeV (XYZ scan, BOSS 703)
scan_names = %w[
  703_3900 703_4009 703_4090 703_4180 703_4190 703_4210
  703_4220 703_4230 703_4245 703_4260 703_4310 703_4360 703_4420
]
scan_data = scan_names.map { |n| DatasetManager.real_data.find(n) }

# Inclusive MC available at 4.009, 4.180, 4.260 and 4.360 GeV only
incmc_names = %w[703_4009 703_4180 703_4260 703_4360]
scan_incMC  = incmc_names.map { |n| DatasetManager.inclusive_mc.find(n) }

### ---------------- Signal decay cards (one per eta_c mode) --------------- ###
# psi(4260) is used as the KKMC top mother (no intermediate resonance in
# e+e- -> pi+ pi- h_c); the prompt photon of h_c -> gamma eta_c is kept.

# --- Mode 01: eta_c -> pi+ pi- pi+ pi- ---
card_m01 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 pi+ pi- pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 02: eta_c -> pi+ pi- pi+ pi- pi0 ---
card_m02 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 pi+ pi- pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 03: eta_c -> pi+ pi- pi+ pi- pi0 pi0 ---
card_m03 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 pi+ pi- pi+ pi- pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 04: eta_c -> pi+ pi- pi+ pi- pi+ pi- ---
card_m04 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 pi+ pi- pi+ pi- pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 05: eta_c -> pi+ pi- pi+ pi- pi+ pi- pi0 ---
card_m05 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 pi+ pi- pi+ pi- pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 06: eta_c -> pi+ pi- pi+ pi- pi+ pi- pi+ pi- ---
card_m06 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 pi+ pi- pi+ pi- pi+ pi- pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 07: eta_c -> pi+ pi- pi+ pi- eta ---
card_m07 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 pi+ pi- pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 08: eta_c -> K+ K- pi+ pi- ---
card_m08 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 K+ K- pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 09: eta_c -> K+ K- pi+ pi- pi0 ---
card_m09 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 K+ K- pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 10: eta_c -> K+ K- pi+ pi- pi+ pi- ---
card_m10 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 K+ K- pi+ pi- pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 11: eta_c -> K+ K- K+ K- ---
card_m11 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 K+ K- K+ K- PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 12: eta_c -> K+ K- pi+ pi- eta ---
card_m12 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 K+ K- pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 13: eta_c -> K_S0 K+ pi- ---
card_m13 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 K_S0 K+ pi- PHSP;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 14: eta_c -> K_S0 K- pi+ ---
card_m14 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 K_S0 K- pi+ PHSP;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 15: eta_c -> p+ anti-p- pi+ pi- ---
card_m15 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 p+ anti-p- pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# --- Mode 16: eta_c -> p+ anti-p- pi+ pi- pi+ pi- ---
card_m16 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0 p+ anti-p- pi+ pi- pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

### --------- Exclusive MC: 50k PHSP events per mode & per energy --------- ###
# create_exclusive_mc_for runs the same signal card on every scan point and
# returns one ExclusiveMC per dataset (events = 50k per energy point).
mc_m01 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_2pip2pim"; c.events = 50_000
  c.decay_card = card_m01; c.cross_section = :default
end
mc_m02 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_2pip2pim_pi0"; c.events = 50_000
  c.decay_card = card_m02; c.cross_section = :default
end
mc_m03 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_2pip2pim_2pi0"; c.events = 50_000
  c.decay_card = card_m03; c.cross_section = :default
end
mc_m04 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_3pip3pim"; c.events = 50_000
  c.decay_card = card_m04; c.cross_section = :default
end
mc_m05 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_3pip3pim_pi0"; c.events = 50_000
  c.decay_card = card_m05; c.cross_section = :default
end
mc_m06 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_4pip4pim"; c.events = 50_000
  c.decay_card = card_m06; c.cross_section = :default
end
mc_m07 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_2pip2pim_eta"; c.events = 50_000
  c.decay_card = card_m07; c.cross_section = :default
end
mc_m08 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_KpKm_pip_pim"; c.events = 50_000
  c.decay_card = card_m08; c.cross_section = :default
end
mc_m09 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_KpKm_pip_pim_pi0"; c.events = 50_000
  c.decay_card = card_m09; c.cross_section = :default
end
mc_m10 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_KpKm_2pip2pim"; c.events = 50_000
  c.decay_card = card_m10; c.cross_section = :default
end
mc_m11 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_2Kp2Km"; c.events = 50_000
  c.decay_card = card_m11; c.cross_section = :default
end
mc_m12 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_KpKm_pip_pim_eta"; c.events = 50_000
  c.decay_card = card_m12; c.cross_section = :default
end
mc_m13 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_KS_Kp_pim"; c.events = 50_000
  c.decay_card = card_m13; c.cross_section = :default
end
mc_m14 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_KS_Km_pip"; c.events = 50_000
  c.decay_card = card_m14; c.cross_section = :default
end
mc_m15 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_ppbar_pip_pim"; c.events = 50_000
  c.decay_card = card_m15; c.cross_section = :default
end
mc_m16 = DatasetManager.create_exclusive_mc_for(scan_data) do |c|
  c.sample_name = "hc_etac_ppbar_2pip2pim"; c.events = 50_000
  c.decay_card = card_m16; c.cross_section = :default
end

### ============================ Event selection =========================== ###
# Common constants for all mode algorithms (see per-algorithm :per_run_ecms note).
ECMS_NOMINAL = [:double, 4.260]

# ---------------------------------------------------------------------------
# Mode 01: eta_c -> pi+ pi- pi+ pi-        (2 pi+, 2 pi-, 1 prompt photon)
# ---------------------------------------------------------------------------
alg_m01 = Algorithm.new("HcEtaC_2pip2pim")
alg_m01.set_header(["HcEtaC_2pip2pimAlg/HcEtaC_2pip2pim.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); the CMS energy is
         read per run, the ECMS constant is only a simulation fallback")

sel_m01 = Selection.new
  .select_track {
    cos_theta 0.93          # |cos(theta)| < 0.93
    Vz        10.0          # |Vz| < 10 cm
    Vr        1.0           # Vr < 1 cm
    nChrp     ">=2"         # >=2 positive tracks
    nChrn     ">=2"         # >=2 negative tracks
    nNet      "==0"         # neutral final state
  }
  .select_photon {
    tdc_emc_start      0
    tdc_emc_end        14
    angle_to_track     10.0      # > 10 deg from any charged track
    energyThreshold_b  0.025     # > 25 MeV (barrel)
    energyThreshold_e  0.050     # > 50 MeV (endcap)
    nGam               ">=1"     # the prompt photon from h_c -> gamma eta_c
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]   # pi+ and pi- (charge conjugation)
    npip ">=2"
    npim ">=2"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) {  # 4C fit of the full final state
    nominal
    constrain_four_momentum
    chi2_cut 200          # loose BOSS cut; tight 35/20 applied later in ROOT
  }
alg_m01.with_decay_card(card_m01).apply(sel_m01)

# ---------------------------------------------------------------------------
# Mode 02: eta_c -> pi+ pi- pi+ pi- pi0     (pi0 from Kalman fit, >=3 photons)
# ---------------------------------------------------------------------------
alg_m02 = Algorithm.new("HcEtaC_2pip2pim_pi0")
alg_m02.set_header(["HcEtaC_2pip2pim_pi0Alg/HcEtaC_2pip2pim_pi0.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m02 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=3"                       # 1 prompt photon + 2 from pi0
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=2"; npim ">=2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # build pi0 (1-C mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25                               # ~ 15 MeV/c^2 mass window
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m02.with_decay_card(card_m02).apply(sel_m02)

# ---------------------------------------------------------------------------
# Mode 03: eta_c -> pi+ pi- pi+ pi- pi0 pi0  (>=5 photons)
# ---------------------------------------------------------------------------
alg_m03 = Algorithm.new("HcEtaC_2pip2pim_2pi0")
alg_m03.set_header(["HcEtaC_2pip2pim_2pi0Alg/HcEtaC_2pip2pim_2pi0.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m03 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=5"                       # 1 prompt + 2 x 2 from the two pi0
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=2"; npim ">=2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m03.with_decay_card(card_m03).apply(sel_m03)

# ---------------------------------------------------------------------------
# Mode 04: eta_c -> pi+ pi- pi+ pi- pi+ pi-   (3 pi+, 3 pi-)
# ---------------------------------------------------------------------------
alg_m04 = Algorithm.new("HcEtaC_3pip3pim")
alg_m04.set_header(["HcEtaC_3pip3pimAlg/HcEtaC_3pip3pim.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m04 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=3"
    nChrn     ">=3"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=3"; npim ">=3"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m04.with_decay_card(card_m04).apply(sel_m04)

# ---------------------------------------------------------------------------
# Mode 05: eta_c -> pi+ pi- pi+ pi- pi+ pi- pi0
# ---------------------------------------------------------------------------
alg_m05 = Algorithm.new("HcEtaC_3pip3pim_pi0")
alg_m05.set_header(["HcEtaC_3pip3pim_pi0Alg/HcEtaC_3pip3pim_pi0.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m05 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=3"
    nChrn     ">=3"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=3"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=3"; npim ">=3"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim, :pip, :pim, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m05.with_decay_card(card_m05).apply(sel_m05)

# ---------------------------------------------------------------------------
# Mode 06: eta_c -> pi+ pi- pi+ pi- pi+ pi- pi+ pi-   (4 pi+, 4 pi-)
# ---------------------------------------------------------------------------
alg_m06 = Algorithm.new("HcEtaC_4pip4pim")
alg_m06.set_header(["HcEtaC_4pip4pimAlg/HcEtaC_4pip4pim.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m06 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=4"
    nChrn     ">=4"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=4"; npim ">=4"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim, :pip, :pim, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m06.with_decay_card(card_m06).apply(sel_m06)

# ---------------------------------------------------------------------------
# Mode 07: eta_c -> pi+ pi- pi+ pi- eta      (eta from Kalman fit)
# ---------------------------------------------------------------------------
alg_m07 = Algorithm.new("HcEtaC_2pip2pim_eta")
alg_m07.set_header(["HcEtaC_2pip2pim_etaAlg/HcEtaC_2pip2pim_eta.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m07 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=3"                            # 1 prompt + 2 from eta
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=2"; npim ">=2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m07.with_decay_card(card_m07).apply(sel_m07)

# ---------------------------------------------------------------------------
# Mode 08: eta_c -> K+ K- pi+ pi-
# ---------------------------------------------------------------------------
alg_m08 = Algorithm.new("HcEtaC_KpKm_pip_pim")
alg_m08.set_header(["HcEtaC_KpKm_pip_pimAlg/HcEtaC_KpKm_pip_pim.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m08 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"; nkm ">=1"
    npip ">=1"; npim ">=1"
  }
  .kinematic_fit([:gamma, :kp, :km, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m08.with_decay_card(card_m08).apply(sel_m08)

# ---------------------------------------------------------------------------
# Mode 09: eta_c -> K+ K- pi+ pi- pi0
# ---------------------------------------------------------------------------
alg_m09 = Algorithm.new("HcEtaC_KpKm_pip_pim_pi0")
alg_m09.set_header(["HcEtaC_KpKm_pip_pim_pi0Alg/HcEtaC_KpKm_pip_pim_pi0.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m09 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=3"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"; nkm ">=1"
    npip ">=1"; npim ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :kp, :km, :pip, :pim, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m09.with_decay_card(card_m09).apply(sel_m09)

# ---------------------------------------------------------------------------
# Mode 10: eta_c -> K+ K- pi+ pi- pi+ pi-     (3 pi+, 3 pi-)
# ---------------------------------------------------------------------------
alg_m10 = Algorithm.new("HcEtaC_KpKm_2pip2pim")
alg_m10.set_header(["HcEtaC_KpKm_2pip2pimAlg/HcEtaC_KpKm_2pip2pim.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m10 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=3"
    nChrn     ">=3"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"; nkm ">=1"
    npip ">=2"; npim ">=2"
  }
  .kinematic_fit([:gamma, :kp, :km, :pip, :pim, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m10.with_decay_card(card_m10).apply(sel_m10)

# ---------------------------------------------------------------------------
# Mode 11: eta_c -> K+ K- K+ K-
# ---------------------------------------------------------------------------
alg_m11 = Algorithm.new("HcEtaC_2Kp2Km")
alg_m11.set_header(["HcEtaC_2Kp2KmAlg/HcEtaC_2Kp2Km.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m11 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp ">=2"; nkm ">=2"
  }
  .kinematic_fit([:gamma, :kp, :km, :kp, :km]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m11.with_decay_card(card_m11).apply(sel_m11)

# ---------------------------------------------------------------------------
# Mode 12: eta_c -> K+ K- pi+ pi- eta
# ---------------------------------------------------------------------------
alg_m12 = Algorithm.new("HcEtaC_KpKm_pip_pim_eta")
alg_m12.set_header(["HcEtaC_KpKm_pip_pim_etaAlg/HcEtaC_KpKm_pip_pim_eta.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m12 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=3"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"; nkm ">=1"
    npip ">=1"; npim ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:gamma, :kp, :km, :pip, :pim, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m12.with_decay_card(card_m12).apply(sel_m12)

# ---------------------------------------------------------------------------
# Mode 13: eta_c -> K_S0 K+ pi-   (K_S0 from secondary-vertex fit)
# ---------------------------------------------------------------------------
alg_m13 = Algorithm.new("HcEtaC_KS_Kp_pim")
alg_m13.set_header(["HcEtaC_KS_Kp_pimAlg/HcEtaC_KS_Kp_pim.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m13 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"           # pi+ (K_S0) + K+
    nChrn     ">=2"           # pi- (K_S0) + pi-
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    identify :kaon, against: [:pion, :proton]
    npip ">=1"; npim ">=2"    # two pions (one from K_S0) + K+
    nkp  ">=1"
  }
  .secondary_vertex_fit([:pip, :pim]) {          # K_S0 -> pi+ pi-
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :K_S0, :kp, :pim]) {   # 4C fit of prompt gamma, K_S0, K+, pi-
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m13.with_decay_card(card_m13).apply(sel_m13)

# ---------------------------------------------------------------------------
# Mode 14: eta_c -> K_S0 K- pi+
# ---------------------------------------------------------------------------
alg_m14 = Algorithm.new("HcEtaC_KS_Km_pip")
alg_m14.set_header(["HcEtaC_KS_Km_pipAlg/HcEtaC_KS_Km_pip.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m14 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"           # pi+ (K_S0) + pi+
    nChrn     ">=2"           # pi- (K_S0) + K-
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    identify :kaon, against: [:pion, :proton]
    npip ">=2"; npim ">=1"
    nkm  ">=1"
  }
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :K_S0, :km, :pip]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m14.with_decay_card(card_m14).apply(sel_m14)

# ---------------------------------------------------------------------------
# Mode 15: eta_c -> p+ anti-p- pi+ pi-
# ---------------------------------------------------------------------------
alg_m15 = Algorithm.new("HcEtaC_ppbar_pip_pim")
alg_m15.set_header(["HcEtaC_ppbar_pip_pimAlg/HcEtaC_ppbar_pip_pim.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m15 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion,   against: [:kaon, :proton]
    nprp ">=1"; nprm ">=1"
    npip ">=1"; npim ">=1"
  }
  .kinematic_fit([:gamma, :prp, :prm, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m15.with_decay_card(card_m15).apply(sel_m15)

# ---------------------------------------------------------------------------
# Mode 16: eta_c -> p+ anti-p- pi+ pi- pi+ pi-
# ---------------------------------------------------------------------------
alg_m16 = Algorithm.new("HcEtaC_ppbar_2pip2pim")
alg_m16.set_header(["HcEtaC_ppbar_2pip2pimAlg/HcEtaC_ppbar_2pip2pim.h"])
       .set_constant({"ECMS" => ECMS_NOMINAL})
       .note(:per_run_ecms, "13 scan points (3.900-4.420 GeV); per-run CMS energy")
sel_m16 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=3"
    nChrn     ">=3"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0; tdc_emc_end 14; angle_to_track 10.0
    energyThreshold_b 0.025; energyThreshold_e 0.050
    nGam ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion,   against: [:kaon, :proton]
    nprp ">=1"; nprm ">=1"
    npip ">=2"; npim ">=2"
  }
  .kinematic_fit([:gamma, :prp, :prm, :pip, :pim, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_m16.with_decay_card(card_m16).apply(sel_m16)

### ------------------------------ Execution ------------------------------ ###
# Real scan data + inclusive MC + the mode-specific signal MC.
root_m01 = alg_m01.execute_on(scan_data + scan_incMC + mc_m01)
root_m02 = alg_m02.execute_on(scan_data + scan_incMC + mc_m02)
root_m03 = alg_m03.execute_on(scan_data + scan_incMC + mc_m03)
root_m04 = alg_m04.execute_on(scan_data + scan_incMC + mc_m04)
root_m05 = alg_m05.execute_on(scan_data + scan_incMC + mc_m05)
root_m06 = alg_m06.execute_on(scan_data + scan_incMC + mc_m06)
root_m07 = alg_m07.execute_on(scan_data + scan_incMC + mc_m07)
root_m08 = alg_m08.execute_on(scan_data + scan_incMC + mc_m08)
root_m09 = alg_m09.execute_on(scan_data + scan_incMC + mc_m09)
root_m10 = alg_m10.execute_on(scan_data + scan_incMC + mc_m10)
root_m11 = alg_m11.execute_on(scan_data + scan_incMC + mc_m11)
root_m12 = alg_m12.execute_on(scan_data + scan_incMC + mc_m12)
root_m13 = alg_m13.execute_on(scan_data + scan_incMC + mc_m13)
root_m14 = alg_m14.execute_on(scan_data + scan_incMC + mc_m14)
root_m15 = alg_m15.execute_on(scan_data + scan_incMC + mc_m15)
root_m16 = alg_m16.execute_on(scan_data + scan_incMC + mc_m16)