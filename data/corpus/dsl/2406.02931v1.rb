# DSL for paper 2406.02931v1: Measurements of the branching fractions of the P-wave
# charmonium spin-singlet state h_c via ψ(3686) → π0h_c
# 4 modes: I(h_c→π+π-π0), II(h_c→K+K-π0), III(h_c→K+K-η), IV(h_c→π+π-η)
# 2712×10⁶ ψ(3686) events

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# ============================================================
# Decay card for ψ(3686) → π0 h_c (shared)
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
# Mode I: h_c → π+π-π0
# ============================================================
decay_card_hc_pipi_pi0 = <<~DECAYCARD
    Alias pi0_bachelor pi0
    Alias pi0_signal pi0

    Decay psi(2S)
    1.0000 pi0_bachelor h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 pi+ pi- pi0_signal PHSP;
    Enddecay

    Decay pi0_bachelor
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0_signal
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_hc_pipi_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_to_pipi_pi0"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = decay_card_hc_pipi_pi0
  config.cross_section = :default
end

alg_hc_pipi_pi0 = Algorithm.new("HcToPiPiPi0")
alg_hc_pipi_pi0.set_header(["HcToPiPiPi0Alg/HcToPiPiPi0.h"])
               .set_constant({"ECMS" => [:double, 3.686]})
               .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 6C kinematic fit")
               .note(:jpsi_veto, "ψ(3686) → π0π0J/ψ and ψ(3686) → ηJ/ψ background vetoed by recoil-mass windows on π0π0 and η; see paper Table 1")
               .note(:competing_hypothesis_veto, "χ²_4C,4γ < χ²_4C,3γ required for photon-count hypothesis test; applied at ROOT level")
               .note(:resonance_veto, "background from ω→π+π-π0, f0(980)→π0_Lπ0_H, and K*(892)→Kπ0_L vetoed; ROOT level")
               .note(:muon_veto, "ψ(3686)→π0π0J/ψ, J/ψ→μ+μ- background vetoed by MUC penetration depth requirement; ROOT level")
               .note(:dalitz_efficiency_correction, "detection efficiency corrected using Dalitz plot distribution from data; ROOT level")
               .note(:fit_procedure, "unbinned ML fit to M(π+π-π0) with signal shape from MC convolved with Gaussian, background by ARGUS function; ROOT level")

sel_hc_pipi_pi0 = Selection.new

sel_hc_pipi_pi0.select_track do
  cos_theta  0.93
  Vz         10.0
  Vr         1.0
  nChrp      ">=1"
  nChrn      ">=1"
  nChrg      "==2"
end

sel_hc_pipi_pi0.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=4"   # 2 for bachelor π0 + 2 for signal π0
end

# PID: pions for Mode I
sel_hc_pipi_pi0.pid(method: :probability) do
  prob_cut 0.001
  identify :pion, against: [:kaon, :proton]
  npip ">=1"
  npim ">=1"
end

sel_hc_pipi_pi0.remove([:pip <= :chrgp, :pim <= :chrgn])

# Reconstruct π0 from γγ pairs (1C Kalman fit)
sel_hc_pipi_pi0.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0 ">=2"   # bachelor π0 + signal π0
end

# 6C kinematic fit: 4C + bachelor π0 mass + signal π0 mass
sel_hc_pipi_pi0.kinematic_fit([:pip, :pim, :pi0, :pi0]) do
  nominal
  constrain_four_momentum
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200  # paper cut < 30; optimal in ROOT
end

# Competing hypothesis fit (4C + 3γ, one fewer photon — no nominal/chi2_cut)
sel_hc_pipi_pi0.kinematic_fit([:pip, :pim, :pi0]) do
  constrain_four_momentum
end

alg_hc_pipi_pi0.with_decay_card(decay_card_hc_pipi_pi0).apply(sel_hc_pipi_pi0)

# ============================================================
# Mode II: h_c → K+K-π0
# ============================================================
decay_card_hc_kk_pi0 = <<~DECAYCARD
    Alias pi0_bachelor pi0
    Alias pi0_signal pi0

    Decay psi(2S)
    1.0000 pi0_bachelor h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 K+ K- pi0_signal PHSP;
    Enddecay

    Decay pi0_bachelor
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0_signal
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_hc_kk_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_to_kk_pi0"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = decay_card_hc_kk_pi0
  config.cross_section = :default
end

alg_hc_kk_pi0 = Algorithm.new("HcToKKPi0")
alg_hc_kk_pi0.set_header(["HcToKKPi0Alg/HcToKKPi0.h"])
             .set_constant({"ECMS" => [:double, 3.686]})
             .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 6C kinematic fit")
             .note(:jpsi_veto, "ψ(3686) → π0π0J/ψ background vetoed by recoil-mass window on π0π0; see paper Table 1")
             .note(:competing_hypothesis_veto, "χ²_4C,3γ < χ²_4C,5γ and χ²_4C,4γ < χ²_4C,5γ required; ROOT level")
             .note(:resonance_veto, "background from f0(980)→π0_Lπ0_H and K*(892)→Kπ0_L vetoed by mass windows; ROOT level")
             .note(:dalitz_efficiency_correction, "detection efficiency corrected using Dalitz plot distribution from data; ROOT level")
             .note(:fit_procedure, "unbinned ML fit to M(K+K-π0) with signal shape from MC convolved with Gaussian, background by ARGUS function; ROOT level")

sel_hc_kk_pi0 = Selection.new

sel_hc_kk_pi0.select_track do
  cos_theta  0.93
  Vz         10.0
  Vr         1.0
  nChrp      ">=1"
  nChrn      ">=1"
  nChrg      "==2"
end

sel_hc_kk_pi0.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=4"
end

# PID: kaons for Mode II
sel_hc_kk_pi0.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp ">=1"
  nkm ">=1"
end

sel_hc_kk_pi0.remove([:kp <= :chrgp, :km <= :chrgn])

sel_hc_kk_pi0.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0 ">=2"
end

# 6C kinematic fit: 4C + bachelor π0 mass + signal π0 mass
sel_hc_kk_pi0.kinematic_fit([:kp, :km, :pi0, :pi0]) do
  nominal
  constrain_four_momentum
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200  # paper cut < 50; optimal in ROOT
end

# Competing hypothesis fits for photon-count veto (3γ, 5γ)
sel_hc_kk_pi0.kinematic_fit([:kp, :km, :pi0]) do
  constrain_four_momentum
end

alg_hc_kk_pi0.with_decay_card(decay_card_hc_kk_pi0).apply(sel_hc_kk_pi0)

# ============================================================
# Mode III: h_c → K+K-η
# ============================================================
decay_card_hc_kk_eta = <<~DECAYCARD
    Alias pi0_bachelor pi0

    Decay psi(2S)
    1.0000 pi0_bachelor h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 K+ K- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0_bachelor
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_hc_kk_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_to_kk_eta"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = decay_card_hc_kk_eta
  config.cross_section = :default
end

alg_hc_kk_eta = Algorithm.new("HcToKKEta")
alg_hc_kk_eta.set_header(["HcToKKEtaAlg/HcToKKEta.h"])
             .set_constant({"ECMS" => [:double, 3.686]})
             .note(:helix_correction, "helix parameter correction applied before 6C kinematic fit")
             .note(:jpsi_veto, "ψ(3686) → ηJ/ψ background vetoed by recoil-mass window on η; see paper Table 1")
             .note(:competing_hypothesis_veto, "χ²_4C,3γ < χ²_4C,5γ and χ²_4C,4γ < χ²_4C,5γ required; ROOT level")
             .note(:gamma_swap_veto, "fake γγ from π0 decay removed by invariant mass requirement; ROOT level")
             .note(:fit_procedure, "unbinned ML fit to M(K+K-η) with MC signal shape convolved with Gaussian, ARGUS background; ROOT level")

sel_hc_kk_eta = Selection.new

sel_hc_kk_eta.select_track do
  cos_theta  0.93
  Vz         10.0
  Vr         1.0
  nChrp      ">=1"
  nChrn      ">=1"
  nChrg      "==2"
end

sel_hc_kk_eta.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=4"   # bachelor π0(2γ) + η(2γ)
end

sel_hc_kk_eta.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp ">=1"
  nkm ">=1"
end

sel_hc_kk_eta.remove([:kp <= :chrgp, :km <= :chrgn])

# Reconstruct bachelor π0
sel_hc_kk_eta.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0 ">=1"
end

# Reconstruct η
sel_hc_kk_eta.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 200
  neta ">=1"
end

# 6C kinematic fit: 4C + bachelor π0 mass + η mass
sel_hc_kk_eta.kinematic_fit([:kp, :km, :pi0, :eta]) do
  nominal
  constrain_four_momentum
  chi2_cut 200  # paper cut < 35; optimal in ROOT
end

# Competing hypothesis fits
sel_hc_kk_eta.kinematic_fit([:kp, :km, :pi0]) do
  constrain_four_momentum
end

alg_hc_kk_eta.with_decay_card(decay_card_hc_kk_eta).apply(sel_hc_kk_eta)

# ============================================================
# Mode IV: h_c → π+π-η
# ============================================================
decay_card_hc_pipi_eta = <<~DECAYCARD
    Alias pi0_bachelor pi0

    Decay psi(2S)
    1.0000 pi0_bachelor h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0_bachelor
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_hc_pipi_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "hc_to_pipi_eta"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = decay_card_hc_pipi_eta
  config.cross_section = :default
end

alg_hc_pipi_eta = Algorithm.new("HcToPiPiEta")
alg_hc_pipi_eta.set_header(["HcToPiPiEtaAlg/HcToPiPiEta.h"])
               .set_constant({"ECMS" => [:double, 3.686]})
               .note(:helix_correction, "helix parameter correction applied before 6C kinematic fit")
               .note(:jpsi_veto, "ψ(3686) → ηJ/ψ background vetoed by recoil-mass window on η; see paper Table 1")
               .note(:competing_hypothesis_veto, "χ²_4C,3γ < χ²_4C,5γ and χ²_4C,4γ < χ²_4C,5γ required; ROOT level")
               .note(:resonance_veto, "K*(892)→Kπ_L vetoed; ROOT level")
               .note(:muon_veto, "ψ(3686)→ηJ/ψ, J/ψ→μ+μ- background vetoed by MUC depth requirement; ROOT level")
               .note(:fit_procedure, "unbinned ML fit to M(π+π-η); upper limit at 90% CL via Bayesian method; ROOT level")

sel_hc_pipi_eta = Selection.new

sel_hc_pipi_eta.select_track do
  cos_theta  0.93
  Vz         10.0
  Vr         1.0
  nChrp      ">=1"
  nChrn      ">=1"
  nChrg      "==2"
end

sel_hc_pipi_eta.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=4"
end

sel_hc_pipi_eta.pid(method: :probability) do
  prob_cut 0.001
  identify :pion, against: [:kaon, :proton]
  npip ">=1"
  npim ">=1"
end

sel_hc_pipi_eta.remove([:pip <= :chrgp, :pim <= :chrgn])

sel_hc_pipi_eta.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0 ">=1"
end

sel_hc_pipi_eta.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 200
  neta ">=1"
end

# 6C kinematic fit: 4C + bachelor π0 mass + η mass
sel_hc_pipi_eta.kinematic_fit([:pip, :pim, :pi0, :eta]) do
  nominal
  constrain_four_momentum
  chi2_cut 200  # paper cut < 24; optimal in ROOT
end

# Competing hypothesis fits
sel_hc_pipi_eta.kinematic_fit([:pip, :pim, :pi0]) do
  constrain_four_momentum
end

alg_hc_pipi_eta.with_decay_card(decay_card_hc_pipi_eta).apply(sel_hc_pipi_eta)

# ============================================================
# Execute all algorithms
# ============================================================
alg_hc_pipi_pi0.execute_on([psip_data, psip_incMC, exMC_hc_pipi_pi0])
alg_hc_kk_pi0.execute_on([psip_data, psip_incMC, exMC_hc_kk_pi0])
alg_hc_kk_eta.execute_on([psip_data, psip_incMC, exMC_hc_kk_eta])
alg_hc_pipi_eta.execute_on([psip_data, psip_incMC, exMC_hc_pipi_eta])