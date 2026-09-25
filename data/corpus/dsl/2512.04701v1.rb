require 'hep_script'

# ============================================================
# Paper: Search for Ξ⁰n → ΛΛX using J/ψ → Ξ⁰Ξ̄⁰ events at BESIII
# arXiv: 2512.04701v1
# Dataset: J/ψ (10087±44)×10⁶ events at √s = 3.097 GeV
# ============================================================

ds_real = DatasetManager.load_real_data(
  "assets/BES3_dataset.md"
)
ds_inc_mc = DatasetManager.load_inclusive_mc(
  "assets/BES3_incMC.md"
)

jpsi_data = ds_real.find("708_3097")
jpsi_inc_mc = ds_inc_mc.find("708_3097")

# Decay card: J/ψ → Ξ⁰ Ξ̄⁰, Ξ̄⁰ → Λ̄ π⁰ (tag side),
# Ξ⁰ → Λ Λ (signal: Ξ⁰ interacts with beam-pipe neutron),
# Λ → p π⁻, Λ̄ → p̄ π⁺, π⁰ → γγ
# Both Λ from Ξ⁰ decay identically → no alias needed.
# Λ and Λ̄ are distinct particles → no alias needed.
decay_card = <<~DECAYCARD
Decay J/psi
1.000 Xi0 anti-Xi0 PHSP;
Enddecay
Decay anti-Xi0
1.000 anti-Lambda pi0 PHSP;
Enddecay
Decay Xi0
1.000 Lambda Lambda PHSP;
Enddecay
Decay anti-Lambda
1.000 anti-p- pi+ PHSP;
Enddecay
Decay Lambda
1.000 p+ pi- PHSP;
Enddecay
Decay pi0
1.000 gamma gamma PHSP;
Enddecay
End
DECAYCARD

jpsi_exc_mc = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "Xi0nToLambdaLambdaX_Jpsi"
  c.related_dataset = jpsi_data
  c.events = 2_250_000
  c.decay_card = decay_card
end

alg = Algorithm.new("Xi0nToLambdaLambdaX")
alg.set_header(["Xi0nToLambdaLambdaXAlg/Xi0nToLambdaLambdaX.h"])
alg.set_constant(ECMS: 3.097)

sel = Selection.new
  # --- Charged tracks ---
  # Final state: 2p (from 2Λ), 1p̄ (from Λ̄), 2π⁻ (from 2Λ), 1π⁺ (from Λ̄)
  # → 3 positive (2p + 1π⁺), 3 negative (2π⁻ + 1p̄)
  .select_track do
    cos_theta 0.93
    nChrp "==3"
    nChrn "==3"
  end
  # --- Photon selection ---
  # E>25 MeV barrel (|cosθ|<0.8), E>50 MeV endcap (0.86<|cosθ|<0.92)
  # EMC time [0,700] ns → tdc_emc [0,14] in units of 50 ns
  # Angle to nearest charged track > 10°
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  # --- PID ---
  # Proton hypothesis has greatest likelihood: L(p) > L(π) and L(p) > L(K)
  # Pion hypothesis has greatest likelihood: L(π) > L(p) and L(π) > L(K)
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion, against: [:kaon, :proton]
    nprp "==2"
    nprm "==1"
    npip "==1"
    npim "==2"
  end
  # --- Tag side: Λ̄ → p̄ π⁺ (secondary vertex fit) ---
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # --- Reconstruct π⁰ → γγ (Kalman 1C fit) ---
  # Mass window [0.11, 0.15] GeV/c² is applied post-fit in ROOT
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  # --- Signal side: two Λ → p π⁻ candidates from remaining tracks ---
  # Build Λ candidates from all remaining p π⁻ pairs (combinatorial).
  # DO NOT remove tracks yet — the best ΛΛ pair is selected below.
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
  end
  # --- Full 4C kinematic fit ---
  # Reconstruct J/ψ → Ξ̄⁰ (→ Λ̄ π⁰) + Ξ⁰ (→ Λ Λ)
  .kinematic_fit([:Lambda_bar, :pi0, :Lambda, :Lambda]) do
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

# --- Inexpressible BOSS-side procedures ---
alg
  .note(:xi0bar_mass_window,
    "Ξ̄⁰ candidate selected as the Λ̄π⁰ combination minimising |M(Λ̄π⁰) − m_Ξ⁰|. " \
    "Ξ̄⁰ signal region: M(Λ̄π⁰) − m_Ξ⁰ ∈ [−15, 10] MeV/c².")
  .note(:xi0_recoil_mass,
    "Recoil mass against Ξ̄⁰ required in [1.295, 1.325] GeV/c² to select Ξ⁰. " \
    "M_recoil(Ξ̄⁰) ≡ √(E_beam² − |p⃗_Ξ̄⁰|²c²)/c².")
  .note(:background_veto,
    "Suppress J/ψ → Ξ⁰Ξ̄⁰, Ξ⁰ → Λπ⁰, Ξ̄⁰ → Λ̄π⁰ background: " \
    "require M_recoil(Ξ̄⁰ Λ_H(L)) < 0 GeV/c². " \
    "When E_miss < |p⃗_miss|, the recoil-mass squared is forced positive; " \
    "recoil mass itself is set negative.")
  .note(:lambda_pair_selection,
    "Two Λ candidates formed from all disjoint pπ⁻ pairs. " \
    "Combination minimising |M(pπ⁻)_H − m_Λ| + |M(pπ⁻)_L − m_Λ| retained. " \
    "Λ_H = higher-momentum Λ, Λ_L = lower-momentum Λ. " \
    "Λ signal region: |M(pπ⁻) − m_Λ| < 3 MeV/c².")
  .note(:rxy_beam_pipe,
    "ΛΛ vertex fit to select beam-pipe origin. " \
    "Signal region R_xy ∈ [2.9, 3.6] cm. " \
    "Inner MDC wall region [6.0, 6.8] cm also visible but not statistically significant.")
  .note(:cross_section_measurement,
    "Cross section σ(Ξ⁰ + ⁹Be → Λ + Λ + X) extracted via unbinned ML fit to R_xy. " \
    "Effective luminosity computed from Ξ⁰ flux, angular distribution, " \
    "attenuation, and target material distribution. ROOT-level only.")

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on([jpsi_data])