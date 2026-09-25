# BESIII DSL for arXiv:2403.03500v1
# Search for h_c hadronic decays via ψ(3686) → π0 h_c
# Modes: 3(π+π-)π0, 2(π+π-)π0η, 2(π+π-)η, p pbar, 2(π+π-)ω

# ============================================================
# Datasets
# ============================================================
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# ============================================================
# Decay card for ψ(3686) → π0 h_c (shared by all modes via Aliases)
# ============================================================
decay_card_psip_to_pi0_hc = <<~DECAYCARD
    Alias pi0_bachelor pi0

    Decay psi(2S)
    1.0000 pi0_bachelor h_c PHSP;
    Enddecay

    Decay pi0_bachelor
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ============================================================
# Mode 1: h_c → 3(π+π-)π0
# ============================================================
decay_card_hc_3pipi_pi0 = <<~DECAYCARD
    Alias pi0_bachelor pi0

    Decay psi(2S)
    1.0000 pi0_bachelor h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 pi+ pi+ pi+ pi- pi- pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0_bachelor
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_hc_3pipi_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_to_3pipi_pi0"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = decay_card_hc_3pipi_pi0
  config.cross_section = :default
end

alg_hc_3pipi_pi0 = Algorithm.new("HcTo3PiPiPi0")
alg_hc_3pipi_pi0.set_header(["HcTo3PiPiPi0Alg/HcTo3PiPiPi0.h"])
                 .set_constant({"ECMS" => [:double, 3.686]})
                 .note(:helix_correction, "helix parameter correction applied to all charged tracks before the (4+N)C kinematic fit")
                 .note(:jpsi_veto, "J/psi-related background vetoed by recoil-mass windows on pi0pi0, pi+pi-, and eta; mass windows listed in paper Table 1")
                 .note(:competing_hypothesis_veto, "events with chi2_4C_nGamma >= chi2_4C_(n-1)Gamma are rejected (photon-count hypothesis test); applied at ROOT level")

sel_hc_3pipi_pi0 = Selection.new

# Charged tracks: pion hypothesis
sel_hc_3pipi_pi0.select_track do
  cos_theta  0.93
  Vz         10.0
  Vr         1.0
  nChrp      "==3"
  nChrn      "==3"
  nNet       "==0"
end

# Photons: E > 25 MeV barrel, > 50 MeV endcap, timing [0,700]ns
sel_hc_3pipi_pi0.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=4"   # 2 for pi0 signal + 2 for bachelor pi0 => 4
end

sel_hc_3pipi_pi0.pid(method: :probability) do
  prob_cut 0.001
  identify :pion, against: [:kaon, :proton]
  npip "==3"
  npim "==3"
end

sel_hc_3pipi_pi0.remove([:pip <= :chrgp, :pim <= :chrgn])

# Reconstruct pi0 candidates from gamma-gamma pairs (1C Kalman fit)
sel_hc_3pipi_pi0.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 20
  npi0 ">=2"   # one bachelor + one from h_c
end

# (4+N)C kinematic fit: 4C + 2 pi0 mass constraints
sel_hc_3pipi_pi0.kinematic_fit([:pip, :pip, :pip, :pim, :pim, :pim, :pi0, :pi0]) do
  nominal
  constrain_four_momentum
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200  # loose; paper cut < 35; optimal cut in ROOT
end

alg_hc_3pipi_pi0.with_decay_card(decay_card_hc_3pipi_pi0).apply(sel_hc_3pipi_pi0)

# ============================================================
# Mode 2: h_c → 2(π+π-)π0η
# ============================================================
decay_card_hc_2pipi_pi0_eta = <<~DECAYCARD
    Alias pi0_bachelor pi0

    Decay psi(2S)
    1.0000 pi0_bachelor h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 pi+ pi+ pi- pi- pi0 eta PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0_bachelor
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_hc_2pipi_pi0_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_to_2pipi_pi0_eta"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = decay_card_hc_2pipi_pi0_eta
  config.cross_section = :default
end

alg_hc_2pipi_pi0_eta = Algorithm.new("HcTo2PiPiPi0Eta")
alg_hc_2pipi_pi0_eta.set_header(["HcTo2PiPiPi0EtaAlg/HcTo2PiPiPi0Eta.h"])
                      .set_constant({"ECMS" => [:double, 3.686]})
                      .note(:helix_correction, "helix parameter correction applied to all charged tracks before (4+N)C kinematic fit")
                      .note(:jpsi_veto, "J/psi-related background vetoed by recoil-mass windows; see Table 1")
                      .note(:competing_hypothesis_veto, "chi2_4C_nGamma < chi2_4C_(n-1)Gamma required; chi2_7C(2pi0 eta) < chi2_7C(3pi0) for pi0 swap background veto; ROOT level")
                      .note(:peaking_background, "peaking background from h_c → gamma eta_c, eta_c → 2(π+π-)eta and eta_c → 2(π+π-π0) included with fixed yields in fit; ROOT level")

sel_hc_2pipi_pi0_eta = Selection.new

sel_hc_2pipi_pi0_eta.select_track do
  cos_theta  0.93
  Vz         10.0
  Vr         1.0
  nChrp      "==2"
  nChrn      "==2"
  nNet       "==0"
end

sel_hc_2pipi_pi0_eta.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=6"   # pi0(2γ) + eta(2γ) + bachelor pi0(2γ) = 6γ
end

sel_hc_2pipi_pi0_eta.pid(method: :probability) do
  prob_cut 0.001
  identify :pion, against: [:kaon, :proton]
  npip "==2"
  npim "==2"
end

sel_hc_2pipi_pi0_eta.remove([:pip <= :chrgp, :pim <= :chrgn])

# Reconstruct pi0 and eta from gamma-gamma pairs
sel_hc_2pipi_pi0_eta.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 20
  npi0 ">=2"
end

sel_hc_2pipi_pi0_eta.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 200
  neta ">=1"
end

# (4+N)C kinematic fit: 4C + pi0 + pi0 + eta mass constraints
sel_hc_2pipi_pi0_eta.kinematic_fit([:pip, :pip, :pim, :pim, :pi0, :pi0, :eta]) do
  nominal
  constrain_four_momentum
  chi2_cut 200  # paper cut < 20; optimal in ROOT
end

alg_hc_2pipi_pi0_eta.with_decay_card(decay_card_hc_2pipi_pi0_eta).apply(sel_hc_2pipi_pi0_eta)

# ============================================================
# Mode 3: h_c → 2(π+π-)η
# ============================================================
decay_card_hc_2pipi_eta = <<~DECAYCARD
    Alias pi0_bachelor pi0

    Decay psi(2S)
    1.0000 pi0_bachelor h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 pi+ pi+ pi- pi- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0_bachelor
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_hc_2pipi_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_to_2pipi_eta"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = decay_card_hc_2pipi_eta
  config.cross_section = :default
end

alg_hc_2pipi_eta = Algorithm.new("HcTo2PiPiEta")
alg_hc_2pipi_eta.set_header(["HcTo2PiPiEtaAlg/HcTo2PiPiEta.h"])
                 .set_constant({"ECMS" => [:double, 3.686]})
                 .note(:helix_correction, "helix parameter correction applied before (4+N)C kinematic fit")
                 .note(:jpsi_veto, "J/psi-related background vetoed by recoil-mass windows; see Table 1")
                 .note(:peaking_background, "peaking background from h_c → gamma eta_c, eta_c → pi+pi-eta' and eta_c → 2(π+π-)π0 with fixed yields in fit; ROOT level")

sel_hc_2pipi_eta = Selection.new

sel_hc_2pipi_eta.select_track do
  cos_theta  0.93
  Vz         10.0
  Vr         1.0
  nChrp      "==2"
  nChrn      "==2"
  nNet       "==0"
end

sel_hc_2pipi_eta.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=4"   # eta(2γ) + bachelor pi0(2γ) = 4γ
end

sel_hc_2pipi_eta.pid(method: :probability) do
  prob_cut 0.001
  identify :pion, against: [:kaon, :proton]
  npip "==2"
  npim "==2"
end

sel_hc_2pipi_eta.remove([:pip <= :chrgp, :pim <= :chrgn])

sel_hc_2pipi_eta.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 20
  npi0 ">=1"
end

sel_hc_2pipi_eta.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 200
  neta ">=1"
end

# (4+N)C kinematic fit: 4C + pi0 + eta mass constraints
sel_hc_2pipi_eta.kinematic_fit([:pip, :pip, :pim, :pim, :pi0, :eta]) do
  nominal
  constrain_four_momentum
  chi2_cut 200  # paper cut < 25; optimal in ROOT
end

alg_hc_2pipi_eta.with_decay_card(decay_card_hc_2pipi_eta).apply(sel_hc_2pipi_eta)

# ============================================================
# Mode 4: h_c → p pbar
# ============================================================
decay_card_hc_ppbar = <<~DECAYCARD
    Alias pi0_bachelor pi0

    Decay psi(2S)
    1.0000 pi0_bachelor h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 p+ anti-p- PHSP;
    Enddecay

    Decay pi0_bachelor
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_hc_ppbar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_to_ppbar"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = decay_card_hc_ppbar
  config.cross_section = :default
end

alg_hc_ppbar = Algorithm.new("HcToPPbar")
alg_hc_ppbar.set_header(["HcToPPbarAlg/HcToPPbar.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .note(:helix_correction, "helix parameter correction applied before 4C kinematic fit")
            .note(:cos_theta_proton, "|cos_theta_p(pbar)| < 0.8 to suppress J/psi lepton pair background")
            .note(:vertex_fit_common, "vertex fit constraining all charged particles to a common vertex before kinematic fit")

sel_hc_ppbar = Selection.new

sel_hc_ppbar.select_track do
  cos_theta  0.93
  Vz         10.0
  Vr         1.0
  nChrp      "==1"
  nChrn      "==1"
  nNet       "==0"
end

# Bachelor pi0 requires >=2 photons
sel_hc_ppbar.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=2"
end

# Proton identification
sel_hc_ppbar.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprp "==1"
  nprm "==1"
end

sel_hc_ppbar.remove([:prp <= :chrgp, :prm <= :chrgn])
            .assign({chrgp: :pip, chrgn: :pim})

# Reconstruct bachelor pi0
sel_hc_ppbar.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 20
  npi0 ">=1"
end

# 4C kinematic fit: p + pbar + pi0
sel_hc_ppbar.kinematic_fit([:prp, :prm, :pi0]) do
  nominal
  constrain_four_momentum
  vertex_fit([0, 1])
  chi2_cut 200  # paper cut < 15; optimal in ROOT
end

alg_hc_ppbar.with_decay_card(decay_card_hc_ppbar).apply(sel_hc_ppbar)

# ============================================================
# Mode 5: h_c → 2(π+π-)ω (ω → π+π-π0)
# ============================================================
decay_card_hc_2pipi_omega = <<~DECAYCARD
    Alias pi0_bachelor pi0

    Decay psi(2S)
    1.0000 pi0_bachelor h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 pi+ pi+ pi- pi- omega PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0_bachelor
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# This mode shares the 3(π+π-)π0 final state with a subset of events
# containing an ω intermediate state; it is extracted via a simultaneous
# fit to ω signal and sideband regions at ROOT level

exMC_hc_2pipi_omega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_to_2pipi_omega"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = decay_card_hc_2pipi_omega
  config.cross_section = :default
end

alg_hc_2pipi_omega = Algorithm.new("HcTo2PiPiOmega")
alg_hc_2pipi_omega.set_header(["HcTo2PiPiOmegaAlg/HcTo2PiPiOmega.h"])
                   .set_constant({"ECMS" => [:double, 3.686]})
                   .note(:helix_correction, "helix parameter correction applied before (4+N)C kinematic fit")
                   .note(:omega_extraction, "h_c → 2(π+π-)ω extracted via simultaneous fit to π+π-π0 invariant mass in ω signal and sideband regions; ROOT level")
                   .note(:jpsi_veto, "J/psi-related background vetoed by recoil-mass windows; see Table 1")

# Same selection as 3(π+π-)π0 mode (same final state, ω identified at ROOT level)
alg_hc_2pipi_omega.with_decay_card(decay_card_hc_2pipi_omega).apply(sel_hc_3pipi_pi0)

# ============================================================
# Execute all algorithms
# ============================================================
alg_hc_3pipi_pi0.execute_on([psip_data, psip_incMC, exMC_hc_3pipi_pi0])
alg_hc_2pipi_pi0_eta.execute_on([psip_data, psip_incMC, exMC_hc_2pipi_pi0_eta])
alg_hc_2pipi_eta.execute_on([psip_data, psip_incMC, exMC_hc_2pipi_eta])
alg_hc_ppbar.execute_on([psip_data, psip_incMC, exMC_hc_ppbar])
alg_hc_2pipi_omega.execute_on([psip_data, psip_incMC, exMC_hc_2pipi_omega])