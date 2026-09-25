# Paper 2504.01823v1: First evidence of ηc→ωφ
# J/ψ data: 10087×10⁶ events at √s=3.097 GeV
# J/ψ→γηc, ηc→ωφ, ω→π⁺π⁻π⁰, φ→K⁺K⁻, π⁰→γγ
# Final state: 3γ π⁺π⁻ K⁺K⁻
# Significance: 4.0σ
# Q-weight method for non-ω/φ background subtraction

# ============================================================
# Decay card
# ============================================================
decay_card = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta_c JPE;
  Enddecay
  Decay eta_c
  1.0000 omega phi HELAMP;
  Enddecay
  Decay omega
  1.0000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay
  Decay phi
  1.0000 K+ K- VSS;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# ============================================================
# Datasets
# ============================================================
data  = DatasetManager.real_data.find("708_3097")
incMC = DatasetManager.inclusive_mc.find("708_3097")

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Jpsi_gamma_etac_omegaphi"
  config.related_dataset = data
  config.events          = 1_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# ============================================================
# Algorithm
# ============================================================
algorithm = Algorithm.new("EtaCToOmegaPhi", "00-00-01")
  .set_header(["EventModel/Event.h", "EvtRecEvent/EvtRecTrack.h"])
  .set_constant(ECMS: 3.0969)

# ============================================================
# Event selection
# ============================================================
event_selection = Selection.new
  # Charged tracks: exactly 4, net charge 0, |cosθ| < 0.93, Vz < 10 cm, Vxy < 1 cm
  .select_track do
    cos_theta 0.93
    Vz         10.0
    Vr         1.0
    nChrp      "==2"
    nChrn      "==2"
    nNet       "==0"
  end
  # Photons: ≥3, standard EMC cuts, |cos(θ_decay)| < 0.95 cut applied later in analysis
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=3"
  end
  # PID for kaons and pions
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkp ">=1"
    nkm ">=1"
  end
  # Remove identified kaons, then identify remaining as pions
  .remove([:kp <= :chrgp, :km <= :chrgn])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon]
    npip ">=1"
    npim ">=1"
  end
  # Reconstruct π⁰ → γγ via Kalman kinematic fit
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  # 5C kinematic fit: K⁺K⁻π⁺π⁻γγγ → ωφ with π⁰ mass constraint
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma]) do
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    nominal
  end

algorithm
  .note(:pi0_decay_angle,
    "|cos(θ_decay)| < 0.95 where cos(θ_decay) = (E_γ1 - E_γ2)/p_γγ. " \
    "Reduces wrong π⁰ combination rate to < 1%.")
  .note(:primary_vertex_fit,
    "Primary vertex fit performed on K⁺K⁻π⁺π⁻ before kinematic fit. " \
    "Successful vertex fit required.")
  .note(:omega_phi_signal_regions,
    "ω signal region: |M_π⁺π⁻π⁰ - M_ω^PDG| < 40.0 MeV/c² (~3σ resolution). " \
    "φ signal region: |M_K⁺K⁻ - M_φ^PDG| < 15.0 MeV/c² (~3σ resolution). " \
    "Applied as mass window cuts on the respective invariant masses.")
  .note(:etap_veto,
    "η′ veto: events with M_π⁺π⁻3γ ∈ [943.0, 969.0] MeV/c² rejected. " \
    "Suppresses J/ψ→φη′, η′→γω background.")
  .note(:qweight,
    "Q-weight method for non-ω/φ background subtraction. " \
    "Phase-space variables: cos(θ_γ), cos(θ_ω), φ_ω, cos(θ_φ), φ_φ, " \
    "cos(θ_K⁺), φ_K⁺, m²(ωφ), m²(γω), m²(γφ), λ_ω/λ_max. " \
    "2-D unbinned fit to M_K⁺K⁻ vs M_π⁺π⁻π⁰ in control sample (200 nearest events). " \
    "ηc signal extracted via fit to M_K⁺K⁻π⁺π⁻π⁰ using M1 transition line shape.")
  .note(:fiveC_optimization,
    "χ²_5C < 20 optimized via FOM S/√(S+B) using signal MC and data.")
  .note(:non_etac_background,
    "Non-ηc background (J/ψ→γωφ continuum) modeled by ARGUS function in ηc fit. " \
    "No interference between ηc signal and continuum 0⁻⁺ considered.")
  .with_decay_card(decay_card)
  .apply(event_selection)
  .execute_on([data, incMC, exMC])