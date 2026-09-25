# DSL: ψ(3686) → γ χ_cJ , χ_cJ → η η η'  (J=0,1,2)
# Paper: 2504.19087v2
# Two independent algorithm chains — one per η' decay mode (Rule T1)
# Mode 1: η' → γ π+ π-  (6C kinematic fit)
# Mode 2: η' → η π+ π-  (7C kinematic fit)
# Each η → γγ

# ============================================================
# Datasets
# ============================================================
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_psip  = DatasetManager.real_data.find("709_3686")
incMC_psip = DatasetManager.inclusive_mc.find("709_3686")

# ============================================================
# Decay cards
# ============================================================
decay_card_mode1 = <<~DECAYCARD
  Decay psi(2S)
  1 gamma chi_cJ PHSP;
  Enddecay
  Decay chi_cJ
  1 eta eta eta' PHSP;
  Enddecay
  Decay eta'
  1 gamma pi+ pi- PHSP;
  Enddecay
  Decay eta
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

decay_card_mode2 = <<~DECAYCARD
  Decay psi(2S)
  1 gamma chi_cJ PHSP;
  Enddecay
  Decay chi_cJ
  1 eta eta eta' PHSP;
  Enddecay
  Decay eta'
  1 eta pi+ pi- PHSP;
  Enddecay
  Decay eta
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# ============================================================
# Signal MC
# ============================================================
sig_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_chicj_etaetaetap_mode1"
  config.related_dataset = data_psip
  config.events          = 100_000
  config.decay_card      = decay_card_mode1
  config.cross_section   = :default
end

sig_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_chicj_etaetaetap_mode2"
  config.related_dataset = data_psip
  config.events          = 100_000
  config.decay_card      = decay_card_mode2
  config.cross_section   = :default
end

# ============================================================
# Algorithm 1: ψ(3686) → γ χ_cJ , χ_cJ → ηηη' , η' → γ π+ π-
# Final state: 8γ + π+π-  (1 radiative γ + 1 γ from η' + 6γ from 3η)
# 6C fit = 4C + 2 η mass constraints
# ============================================================
alg_mode1 = Algorithm.new("PsipToGamChicJ_EtapToGamPiPi")
alg_mode1.set_header(["PsipToGamChicJAlg/PsipToGamChicJ.h"])
          .set_constant({ "ECMS" => [:double, 3.686] })

sel_mode1 = Selection.new
  # 2 charged pions
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrg "==2"
    NetCharge 0
  end
  # At least 8 photons (1 radiative + 1 from η' + 6 from 3η)
  .select_photon do
    nGam ">=8"
    min_energy 0.025
    min_angle 10.0
  end
  # Pion PID
  .pid do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  # π0 veto: reject events where any γγ pair forms π0
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  end
  # Reconstruct 3 η → γγ via Kalman 1C fits
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  end
  # Competing photon hypothesis veto: 5γ fit (drop one η)
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :eta, :eta]) do
    constrain_four_momentum
  end
  # Competing photon hypothesis veto: 6γ fit (drop two η)
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :eta]) do
    constrain_four_momentum
  end
  # Competing photon hypothesis veto: 7γ fit (drop a different η combination)
  .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim, :eta, :eta]) do
    constrain_four_momentum
  end
  # Nominal 6C fit: radiative γ + 3η + η'(→γπ+π-)
  # 4C (four-momentum) + 2 η mass constraints
  .kinematic_fit([:gamma, :eta, :eta, :eta, :pip, :pim, :gamma]) do
    constrain_four_momentum
    invariant_mass_of(:eta).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:eta).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    nominal
  end

alg_mode1.with_decay_card(decay_card_mode1).apply(sel_mode1)
alg_mode1.note(:signal_description, "ψ(3686) → γ χ_cJ , χ_cJ → ηηη' , η' → γ π+ π-. χ_cJ = χ_c0, χ_c1, χ_c2.")
alg_mode1.note(:kinematic_fit, "6C fit = 4C (four-momentum conservation) + 2 η mass constraints. χ² cut optimized in ROOT.")
alg_mode1.note(:competing_hypothesis, "χ²_5γ < χ²_nominal → reject; χ²_6γ < χ²_nominal → reject; χ²_7γ < χ²_nominal → reject. Evaluated in ROOT by comparing stored χ² values.")
alg_mode1.note(:pi0_veto, "Events where any γγ combination forms π0 (|M_γγ - m_π0| < 15 MeV) are rejected. Applied in ROOT via stored π0 mass.")
alg_mode1.note(:eta_reconstruction, "η → γγ reconstructed via 1C Kalman fit with |M_γγ - m_η| < 50 MeV. Three η candidates reconstructed per event.")
alg_mode1.note(:etap_reconstruction, "η' reconstructed from γ π+ π- with |M_γππ - m_η'| < 30 MeV.")
alg_mode1.execute_on([data_psip, incMC_psip, sig_mode1])

# ============================================================
# Algorithm 2: ψ(3686) → γ χ_cJ , χ_cJ → ηηη' , η' → η π+ π-
# Final state: 9γ + π+π- + π+π-  (1 radiative γ + 8γ from 4η + 2π+2π from η')
# 7C fit = 4C + 3 η mass constraints
# ============================================================
alg_mode2 = Algorithm.new("PsipToGamChicJ_EtapToEtaPiPi")
alg_mode2.set_header(["PsipToGamChicJAlg/PsipToGamChicJ.h"])
          .set_constant({ "ECMS" => [:double, 3.686] })

sel_mode2 = Selection.new
  # 4 charged pions (2 from η' → ηπ+π-, then η → γγ only adds photons)
  # Wait: final state is 2π from the η' decay chain
  # Actually: η' → η π+ π-, so we have π+ and π-
  # The other 3 η → γγ add no charged tracks
  # So only 2 charged tracks, not 4
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrg "==2"
    NetCharge 0
  end
  # At least 9 photons (1 radiative + 8 from 4η)
  .select_photon do
    nGam ">=9"
    min_energy 0.025
    min_angle 10.0
  end
  # Pion PID
  .pid do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  # π0 veto
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  end
  # Reconstruct 4 η → γγ via Kalman 1C fits
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  end
  # Competing photon hypothesis veto: 6γ fit (drop two η)
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :eta, :eta]) do
    constrain_four_momentum
  end
  # Competing photon hypothesis veto: 7γ fit (drop one η)
  .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim, :eta, :eta]) do
    constrain_four_momentum
  end
  # Competing photon hypothesis veto: 8γ fit (drop a different η combination)
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :eta, :eta]) do
    constrain_four_momentum
  end
  # Nominal 7C fit: radiative γ + 4η + π+π-
  # 4C (four-momentum) + 3 η mass constraints
  .kinematic_fit([:gamma, :eta, :eta, :eta, :eta, :pip, :pim]) do
    constrain_four_momentum
    invariant_mass_of(:eta).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:eta).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:eta).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    nominal
  end

alg_mode2.with_decay_card(decay_card_mode2).apply(sel_mode2)
alg_mode2.note(:signal_description, "ψ(3686) → γ χ_cJ , χ_cJ → ηηη' , η' → η π+ π-. χ_cJ = χ_c0, χ_c1, χ_c2.")
alg_mode2.note(:kinematic_fit, "7C fit = 4C (four-momentum conservation) + 3 η mass constraints. χ² cut optimized in ROOT.")
alg_mode2.note(:competing_hypothesis, "χ²_6γ < χ²_nominal → reject; χ²_7γ < χ²_nominal → reject; χ²_8γ < χ²_nominal → reject. Evaluated in ROOT by comparing stored χ² values.")
alg_mode2.note(:pi0_veto, "Events where any γγ combination forms π0 (|M_γγ - m_π0| < 15 MeV) are rejected. Applied in ROOT via stored π0 mass.")
alg_mode2.note(:eta_reconstruction, "η → γγ reconstructed via 1C Kalman fit with |M_γγ - m_η| < 50 MeV. Four η candidates reconstructed per event.")
alg_mode2.note(:etap_reconstruction, "η' reconstructed from η π+ π- with |M_ηππ - m_η'| < 30 MeV.")
alg_mode2.execute_on([data_psip, incMC_psip, sig_mode2])