# DSL: ω → π+π+e-e- LNV search via J/ψ → ω η
# Paper: 2504.21539v1
# Reference channel: ω → π+π-π0 (η → γγ)
# Signal channel: ω → π+π+e-e- (η → γγ)
# Two algorithms — Rule T1

# ============================================================
# Datasets
# ============================================================
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_jpsi  = DatasetManager.real_data.find("708_3097")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")

# ============================================================
# Decay cards
# ============================================================

# Reference: J/ψ → ωη, ω → π+π-π0, η → γγ
decay_card_ref = <<~DECAYCARD
  Decay J/psi
  1 omega eta PHSP;
  Enddecay
  Decay omega
  1 pi+ pi- pi0 PHSP;
  Enddecay
  Decay eta
  1 gamma gamma PHSP;
  Enddecay
  Decay pi0
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Signal: J/ψ → ωη, ω → e+e-π+π-, η → γγ
decay_card_sig = <<~DECAYCARD
  Decay J/psi
  1 omega eta PHSP;
  Enddecay
  Decay omega
  1 e+ e- pi+ pi- PHSP;
  Enddecay
  Decay eta
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# ============================================================
# Signal MC
# ============================================================
sig_ref = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_omega_eta_ref"
  config.related_dataset = data_jpsi
  config.events          = 100_000
  config.decay_card      = decay_card_ref
  config.cross_section   = :default
end

sig_sig = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_omega_eta_sig"
  config.related_dataset = data_jpsi
  config.events          = 100_000
  config.decay_card      = decay_card_sig
  config.cross_section   = :default
end

# ============================================================
# Algorithm 1: Reference channel — ω → π+π-π0 (η → γγ)
# Final state: π+π- + 4γ (2 from π0, 2 from η)
# 5C fit = 4C + π0 mass constraint, χ² < 20
# ============================================================
alg_ref = Algorithm.new("JpsiToOmegaEta_Ref")
alg_ref.set_header(["JpsiToOmegaEtaAlg/JpsiToOmegaEta.h"])
        .set_constant({ "ECMS" => [:double, 3.097] })

sel_ref = Selection.new
  # 2 charged tracks, net zero charge
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrg "==2"
    NetCharge 0
  end
  # At least 4 photons (π0→γγ + η→γγ)
  .select_photon do
    nGam ">=4"
    min_energy 0.025
    min_angle 10.0
  end
  # Pion PID
  .pid do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  # Reconstruct π0 → γγ via Kalman 1C fit
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
  end
  # Reconstruct η → γγ via Kalman 1C fit
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  end
  # Competing hypothesis: 4 pions + 2 photons (KKKK γγ → all π mass)
  .kinematic_fit([:pip, :pim, :pip, :pim, :gamma, :gamma]) do
    constrain_four_momentum
  end
  # Competing hypothesis: 4 kaons + 2 photons (KKKK γγ)
  .kinematic_fit([:kp, :km, :kp, :km, :gamma, :gamma]) do
    constrain_four_momentum
  end
  # Nominal 5C fit: π+π- + π0 + η
  # 4C (four-momentum) + π0 mass constraint
  .kinematic_fit([:pip, :pim, :pi0, :eta]) do
    constrain_four_momentum
    invariant_mass_of(:pi0).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    nominal
  end
  # η mass constraint applied in ROOT via stored mass

alg_ref.with_decay_card(decay_card_ref).apply(sel_ref)
alg_ref.note(:signal_description, "Reference channel: J/ψ → ωη, ω → π+π-π0, η → γγ.")
alg_ref.note(:kinematic_fit, "5C fit = 4C (four-momentum conservation) + π0 mass constraint. χ² < 20. η mass window applied in ROOT.")
alg_ref.note(:competing_hypothesis, "χ²(ππππγγ) < χ²_nominal → reject; χ²(KKKKγγ) < χ²_nominal → reject. Evaluated in ROOT by comparing stored χ² values.")
alg_ref.note(:mass_windows, "π0 mass window |M_γγ - m_π0| < 15 MeV; η mass window |M_γγ - m_η| < 40 MeV. Applied in ROOT via stored invariant masses.")
alg_ref.execute_on([data_jpsi, incMC_jpsi, sig_ref])

# ============================================================
# Algorithm 2: Signal channel — ω → π+π+e-e- (η → γγ)
# Final state: e+e-π+π- + 2γ (from η)
# 4C fit, χ² < 10
# ============================================================
alg_sig = Algorithm.new("JpsiToOmegaEta_Sig")
alg_sig.set_header(["JpsiToOmegaEtaAlg/JpsiToOmegaEta.h"])
        .set_constant({ "ECMS" => [:double, 3.097] })

sel_sig = Selection.new
  # 4 charged tracks, net zero charge
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrg "==4"
    NetCharge 0
  end
  # At least 2 photons (η → γγ)
  .select_photon do
    nGam ">=2"
    min_energy 0.025
    min_angle 10.0
  end
  # PID: identify electrons and pions
  # e: CL_e > 0.001 and CL_e/(CL_e+CL_K+CL_π) > 0.8
  .pid do
    prob_cut 0.001
    identify :electron, against: [:kaon, :pion]
    identify :pion, against: [:kaon, :proton]
  end
  # Reconstruct η → γγ via Kalman 1C fit
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
  end
  # Competing hypothesis: 4 pions + 2 photons
  .kinematic_fit([:pip, :pim, :pip, :pim, :gamma, :gamma]) do
    constrain_four_momentum
  end
  # Competing hypothesis: 4 kaons + 2 photons
  .kinematic_fit([:kp, :km, :kp, :km, :gamma, :gamma]) do
    constrain_four_momentum
  end
  # Competing hypothesis: 2 pions + 2 kaons + 2 photons
  .kinematic_fit([:pip, :pim, :kp, :km, :gamma, :gamma]) do
    constrain_four_momentum
  end
  # Competing hypothesis: 2 pions + p pbar + 2 photons
  .kinematic_fit([:pip, :pim, :prp, :prm, :gamma, :gamma]) do
    constrain_four_momentum
  end
  # Nominal 4C fit: e+e-π+π- + η
  # η mass constraint NOT applied in fit — only in ROOT
  .kinematic_fit([:ep, :em, :pip, :pim, :eta]) do
    constrain_four_momentum
    chi2_cut 10
    nominal
  end

alg_sig.with_decay_card(decay_card_sig).apply(sel_sig)
alg_sig.note(:signal_description, "Signal channel: J/ψ → ωη, ω → π+π+e-e- (LNV), η → γγ.")
alg_sig.note(:electron_pid, "CL_e > 0.001 and CL_e/(CL_e+CL_K+CL_π) > 0.8. Applied via PID probability method; likelihood ratio check in ROOT.")
alg_sig.note(:kinematic_fit, "4C fit (four-momentum conservation). χ² < 10. η mass window applied in ROOT.")
alg_sig.note(:competing_hypothesis, "χ²(ππππγγ) < χ²_nominal → reject; χ²(KKKKγγ) < χ²_nominal → reject; χ²(ππKKγγ) < χ²_nominal → reject; χ²(ππppγγ) < χ²_nominal → reject. Evaluated in ROOT by comparing stored χ² values.")
alg_sig.note(:mass_windows, "η mass window |M_γγ - m_η| < 40 MeV; M_π+π-e+e- > 1.0 GeV/c². Applied in ROOT via stored invariant masses.")
alg_sig.note(:lnv_signal, "Lepton number violation signal extracted from M_π+π-e+e- spectrum in ROOT above 1.0 GeV/c².")
alg_sig.execute_on([data_jpsi, incMC_jpsi, sig_sig])