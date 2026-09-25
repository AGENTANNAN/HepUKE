# =============================================================================
# ψ(3686) → Σ⁰ Σ̄⁰ ω   (Σ⁰→γΛ, Λ→pπ⁻; Σ̄⁰→γΛ̄, Λ̄→p̄π⁺; ω→π⁺π⁻π⁰, π⁰→γγ)
# BOSS part: dataset preparation + event selection up to the final 5C kinematic fit
# =============================================================================

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data, √s = 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

# Decay card for the signal process (EvtGen syntax)
decay_card_signal = <<~DECAYCARD
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
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  Decay omega
  1.0000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 3M-event exclusive signal MC with the full decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_Sigma0Sigmabar0omega"
  config.related_dataset = psip_data
  config.events          = 3_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "Sigma0Sigmabar0Omega"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
  # ---- Charged track selection ----
  .select_track {
    cos_theta 0.93     # |cosθ| < 0.93
    Vz        10.0     # |Vz| < 10 cm
    Vr        1.0      # Vr < 1 cm
    nChrp     ">=3"    # at least 3 positively charged tracks
    nChrn     ">=3"    # at least 3 negatively charged tracks
  }
  # ---- Photon selection ----
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025   # barrel  E > 25 MeV
    energyThreshold_e 0.050   # endcap  E > 50 MeV
    nGam              ">=4"   # at least 4 photons
  }
  # ---- PID: identify (anti-)protons against kaons and pions ----
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p+ and p̄ (charge-conjugation shorthand)
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])     # remove (anti-)protons from the generic charged lists
  # ---- PID: remaining tracks identified as pions against kaons and protons ----
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]   # π+ and π⁻
    npip ">=1"
    npim ">=1"
  }
  # ---- Isolated photons: suppress showers from (anti-)proton ----
  .select_isolated_photon {
    angle_to_prp_track 20.0
    angle_to_prm_track 20.0
    nGam               ">=4"
  }
  # ---- Λ → pπ⁻ secondary vertex fit ----
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # ---- Λ̄ → p̄π⁺ secondary vertex fit ----
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # ---- π⁰ → γγ reconstructed with a Kalman (1C) mass constraint ----
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # ---- Final fit: 5C (4-momentum conservation + π⁰ mass constraint), nominal ----
  .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim, :gamma, :gamma]) {
    nominal
    vertex_fit([2, 3])   # π⁺ (index 2) and π⁻ (index 3) from ω share a common vertex
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 60
  }
  # ---- Competing hypothesis: pure 4C fit (no χ² cut, not nominal) ----
  # Stores chi2_4c so that events better explained by a wrong photon multiplicity
  # can be rejected against chi2_5c at the ROOT level.
  .kinematic_fit([:Lambda, :Lambda_bar, :pip, :pim, :gamma, :gamma]) {
    constrain_four_momentum
  }

my_algorithm
  .note(:lambda_mass_window, "Λ / Λ̄ candidates are required to satisfy |M(pπ) − m_Λ| < 0.01 GeV/c² and a positive decay length; the secondary_vertex_fit selects the combination by minimum mass difference, the mass window and decay-length sign cut are applied downstream")
  .note(:background_veto, "J/ψ veto applied on the recoil mass of π⁺π⁻(π⁰γᵢγⱼ): events with the recoil mass in [3.089,3.105] GeV/c² (π⁺π⁻J/ψ) or [3.090,3.145] GeV/c² (π⁰π⁰J/ψ) are rejected")
  .note(:four_c_fit_veto, "the non-nominal 4C fit stores chi2_4c and is used to reject events with wrong photon multiplicity by comparing chi2_4c with chi2_5c in ROOT")
  .note(:sigma0_photon_selection, "the Σ⁰ and Σ̄⁰ photons are chosen among the photons not used in the π⁰ reconstruction by minimising Δ = (M(γΛ) − M(γΛ̄))²")
  .note(:sigma0_signal_region, "2-D Σ⁰Σ̄⁰ signal region M(γΛ) ∈ [1.176,1.204] GeV/c² and M(γΛ̄) ∈ [1.175,1.204] GeV/c² selected in ROOT")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC and the exclusive signal MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])