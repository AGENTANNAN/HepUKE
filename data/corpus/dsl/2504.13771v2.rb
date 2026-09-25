# DSL: J/ψ → K_S0 K_S0  and  ψ(3686) → K_S0 K_S0  (CP violation search)
# Paper: 2504.13771v2
# Two separate algorithms — datasets differ (J/ψ vs ψ(2S)) and χ² cuts differ

# ============================================================
# Datasets
# ============================================================
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_jpsi  = DatasetManager.real_data.find("708_3097")
data_psip  = DatasetManager.real_data.find("709_3686")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")
incMC_psip = DatasetManager.inclusive_mc.find("709_3686")

# ============================================================
# Decay cards — K_S0 K_S0  (KKMC + J/ψ / psi(2S) top mothers)
# ============================================================
decay_card_jpsi_ksks = <<~DECAYCARD
  Decay J/psi
  1 K_S0 K_S0 PHSP;
  Enddecay
  Decay K_S0
  1 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

decay_card_psip_ksks = <<~DECAYCARD
  Decay psi(2S)
  1 K_S0 K_S0 PHSP;
  Enddecay
  Decay K_S0
  1 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# ============================================================
# Signal MC
# ============================================================
sig_jpsi_ksks = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_ksks"
  config.related_dataset = data_jpsi
  config.events          = 100_000
  config.decay_card      = decay_card_jpsi_ksks
  config.cross_section   = :default
end

sig_psip_ksks = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psip_ksks"
  config.related_dataset = data_psip
  config.events          = 100_000
  config.decay_card      = decay_card_psip_ksks
  config.cross_section   = :default
end

# ============================================================
# Algorithm 1: J/ψ → K_S0 K_S0
# Final state: 4 charged pions → 2 K_S0
# ============================================================
alg_jpsi = Algorithm.new("JpsiToKsKs")
alg_jpsi.set_header(["JpsiToKsKsAlg/JpsiToKsKs.h"])
         .set_constant({ "ECMS" => [:double, 3.097] })

sel_jpsi = Selection.new
  # 4 charged tracks, net zero charge, pion PID
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrg "==4"
    NetCharge 0
  end
  .select_photon do
    nGam ">=1"
    min_energy 0.025
    min_angle 10.0
  end
  .pid do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  # Reconstruct two K_S0 via secondary vertex fits
  # Mass window |M_ππ - m_K_S0| < 18 MeV/c², decay length > 2σ
  .secondary_vertex_fit(:pip, :pim) do
    build_virtual_particle :K_S0, by_minimizing_mass_difference
    chi2_cut 200
  end
  .secondary_vertex_fit(:pip, :pim) do
    build_virtual_particle :K_S0, by_minimizing_mass_difference
    chi2_cut 200
  end
  # Competing-hypothesis veto: γ K_S0 K_S0
  .kinematic_fit([:gamma, :K_S0, :K_S0]) do
    constrain_four_momentum
  end
  # Nominal 4C fit: K_S0 K_S0  (χ² < 15)
  .kinematic_fit([:K_S0, :K_S0]) do
    constrain_four_momentum
    chi2_cut 15
    nominal
  end

alg_jpsi.with_decay_card(decay_card_jpsi_ksks).apply(sel_jpsi)
alg_jpsi.note(:mass_window, "K_S0 mass window: |M_ππ - m_K_S0| < 18 MeV/c². Applied in ROOT via stored invariant masses.")
alg_jpsi.note(:decay_length, "K_S0 decay length > 2× vertex resolution. Applied in ROOT via stored vertex info.")
alg_jpsi.note(:competing_hypothesis, "χ²_4C(γK_S0 K_S0) < χ²_4C(K_S0 K_S0) → rejected. Evaluated in ROOT by comparing stored χ² values.")
alg_jpsi.note(:signal_region, "Signal extracted in ROOT from K_S0 momentum space distribution.")
alg_jpsi.execute_on([data_jpsi, incMC_jpsi, sig_jpsi_ksks])

# ============================================================
# Algorithm 2: ψ(3686) → K_S0 K_S0
# Final state: 4 charged pions → 2 K_S0
# Wider mass window (30 MeV/c²) and looser χ² (< 30)
# ============================================================
alg_psip = Algorithm.new("PsipToKsKs")
alg_psip.set_header(["PsipToKsKsAlg/PsipToKsKs.h"])
         .set_constant({ "ECMS" => [:double, 3.686] })

sel_psip = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrg "==4"
    NetCharge 0
  end
  .select_photon do
    nGam ">=1"
    min_energy 0.025
    min_angle 10.0
  end
  .pid do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  # Reconstruct two K_S0 via secondary vertex fits
  # Mass window |M_ππ - m_K_S0| < 30 MeV/c², decay length > 2σ
  .secondary_vertex_fit(:pip, :pim) do
    build_virtual_particle :K_S0, by_minimizing_mass_difference
    chi2_cut 200
  end
  .secondary_vertex_fit(:pip, :pim) do
    build_virtual_particle :K_S0, by_minimizing_mass_difference
    chi2_cut 200
  end
  # Competing-hypothesis veto: γ K_S0 K_S0
  .kinematic_fit([:gamma, :K_S0, :K_S0]) do
    constrain_four_momentum
  end
  # Nominal 4C fit: K_S0 K_S0  (χ² < 30)
  .kinematic_fit([:K_S0, :K_S0]) do
    constrain_four_momentum
    chi2_cut 30
    nominal
  end

alg_psip.with_decay_card(decay_card_psip_ksks).apply(sel_psip)
alg_psip.note(:mass_window, "K_S0 mass window: |M_ππ - m_K_S0| < 30 MeV/c². Applied in ROOT via stored invariant masses.")
alg_psip.note(:decay_length, "K_S0 decay length > 2× vertex resolution. Applied in ROOT via stored vertex info.")
alg_psip.note(:competing_hypothesis, "χ²_4C(γK_S0 K_S0) < χ²_4C(K_S0 K_S0) → rejected. Evaluated in ROOT by comparing stored χ² values.")
alg_psip.note(:signal_region, "Signal extracted in ROOT from K_S0 momentum space distribution.")
alg_psip.execute_on([data_psip, incMC_psip, sig_psip_ksks])