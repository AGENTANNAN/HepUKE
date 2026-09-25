# Paper 2503.18620v1: First observation of ψ(3686)→Σ⁰Σ̄⁰ω
# ψ(2S) data: 27.12×10⁸ events at √s=3.686 GeV
# Σ⁰→γΛ, Λ→pπ⁻, Σ̄⁰→γΛ̄, Λ̄→p̄π⁺, ω→π⁺π⁻π⁰, π⁰→γγ
# Significance: 8.9σ

# ============================================================
# Decay card
# ============================================================
decay_card = <<~DECAYCARD
  Decay psi(2S)
  1.0000 Sigma0 anti-Sigma0 omega PHSP;
  Enddecay
  Decay Sigma0
  1.0000 gamma Lambda0 PHSP;
  Enddecay
  Decay anti-Sigma0
  1.0000 gamma anti-Lambda0 PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi- PHSP;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+ PHSP;
  Enddecay
  Decay omega
  1.0000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# ============================================================
# Datasets
# ============================================================
data  = DatasetManager.real_data.find("709_3686")
incMC = DatasetManager.inclusive_mc.find("709_3686")

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psi2S_Sigma0Sigma0omega"
  config.related_dataset = data
  config.events          = 3_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# ============================================================
# Algorithm
# ============================================================
algorithm = Algorithm.new("PsipToSigma0Sigma0omega", "00-00-01")
  .set_header(["EventModel/Event.h", "EvtRecEvent/EvtRecTrack.h"])
  .set_constant(ECMS: 3.686109)

# ============================================================
# Event selection
# ============================================================
event_selection = Selection.new
  # Charged tracks: ≥3 pos, ≥3 neg, |cosθ| < 0.93
  .select_track do
    cos_theta 0.93
    Vz         10.0
    Vr         1.0
    nChrp      ">=3"
    nChrn      ">=3"
  end
  # Photons: ≥4, standard EMC cuts
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=4"
  end
  # Isolated photons: anti-proton opening angle > 20°
  .select_isolated_photon do
    angle_to_prm_track 20.0
    angle_to_prp_track 20.0
    nGam               ">=4"
  end
  # PID for protons: P(p) > P(K) and P(p) > P(π)
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  end
  # Remove identified protons from generic track lists
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  # PID for remaining tracks as pions
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  # Reconstruct Λ → pπ⁻ via secondary vertex fit
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # Reconstruct Λ̄ → p̄π⁺ via secondary vertex fit
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # Reconstruct π⁰ → γγ via Kalman kinematic fit
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  # 5C kinematic fit: Λ Λ̄ π⁺π⁻ γγ (π⁰ mass constraint + 4-momentum conservation)
  .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim, :gamma, :gamma]) do
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 60
    nominal
  end

algorithm
  .note(:sigma0_photon_selection,
    "Σ⁰ photon selected by minimizing Δ = (M_γⁱΛ - M_γʲΛ̄)² among photons not used in π⁰ reconstruction")
  .note(:lambda_mass_window,
    "Λ/Λ̄ candidates required: |M_pπ⁻(p̄π⁺) - m_Λ(Λ̄)| < 0.01 GeV/c² (~5σ); decay length > 0")
  .note(:jpsi_veto,
    "J/ψ veto: recoil mass of π⁺π⁻(π⁰γⁱγʲ) outside [3.089,3.105] GeV/c² for π⁺π⁻J/ψ and [3.090,3.145] GeV/c² for π⁰π⁰J/ψ")
  .note(:fourC_veto,
    "4C kinematic fit veto: χ²_4C(ΛΛ̄π⁺π⁻4γ) < χ²_4C(ΛΛ̄π⁺π⁻3γ) AND < χ²_4C(ΛΛ̄π⁺π⁻5γ) to reject wrong photon multiplicities")
  .note(:omega_pion_vertex,
    "Vertex fit performed on each π⁺π⁻ pair from ω to ensure common decay vertex")
  .note(:sigma0_mass_window,
    "2-D Σ⁰Σ̄⁰ signal region per Table I: M_γΛ ∈ [1.176,1.204], M_γΛ̄ ∈ [1.175,1.204] GeV/c²")
  .note(:fiveC_chi2_optimization,
    "χ²_5C < 60 optimized via figure-of-merit S/√(S+B) using signal and inclusive MC")
  .with_decay_card(decay_card)
  .apply(event_selection)
  .execute_on([data, incMC, exMC])