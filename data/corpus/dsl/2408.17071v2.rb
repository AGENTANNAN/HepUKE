# Search for hc → π⁺π⁻J/ψ via ψ(3686) → π⁰hc
# ψ(3686), 2712M events, arXiv:2408.17071v2

psi3686_data = DatasetManager.real_data.find("709_3686")
psi3686_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
  Decay psi(2S)
  1.000 pi0 hc PHSP;
  Enddecay

  Decay hc
  1.000 pi+ pi- J/psi PHSP;
  Enddecay

  Decay J/psi
  1.000 e+ e- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psi2s_pi0_hc_pipi_jpsi"
  config.related_dataset = psi3686_data
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

# Algorithm: ψ(3686)→π⁰hc, hc→π⁺π⁻J/ψ, J/ψ→ℓ⁺ℓ⁻
# Lepton PID via momentum threshold (p>1.0 GeV/c → lepton);
# e/μ separation via EMC energy (E>1.0 GeV → e, else → μ) applied in ROOT
alg = Algorithm.new("Psi2SPi0HcPiPiJpsi")
alg.set_header(["Psi2SPi0HcPiPiJpsiAlg/Psi2SPi0HcPiPiJpsi.h"])
   .set_constant({ "ECMS" => [:double, 3.686] })

# ── Selection ──
event_selection = Selection.new

event_selection
  .select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrg "==4"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  }
  # Hybrid PID: high-momentum leptons (p > 1.0 GeV/c) for J/ψ→ℓ⁺ℓ⁻,
  # pions (p<1.0 GeV/c) for hc→π⁺π⁻. e/μ separation by EMC energy in ROOT.
  .pid(method: :probability) {
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 1.0
    identify :pion, against: [:kaon]
  }
  # Reconstruct π⁰ → γγ via 1C kalman fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # 5C kinematic fit: 4-momentum conservation + π⁰ mass constraint
  # Uses lp/lm from hybrid PID (both e and μ hypotheses contribute)
  .kinematic_fit([:pip, :pim, :lp, :lm, :pi0]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 15
  }
  # Competing hypothesis: η → π⁺π⁻π⁰ (store χ² for ROOT veto)
  .kinematic_fit([:pip, :pim, :pi0]) {
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  }

alg.with_decay_card(decay_card).apply(event_selection)

# Additional cuts applied in ROOT:
# - J/ψ mass window: [3.085, 3.108] GeV/c² for ℓ⁺ℓ⁻ pair
# - e/μ separation: E>1.0 GeV → electron, E<0.4 GeV → muon
# - |M(π⁺π⁻π⁰) - M(η)| > 24.5 MeV/c² (η veto)
# - M(π⁺π⁻) > 0.3 GeV/c² (γ→e⁺e⁻ conversion suppression)
alg.note(:root_cuts, "ROOT-level: J/ψ mass window [3.085,3.108] GeV/c², e/μ separation (E>1.0 GeV→e, E<0.4 GeV→μ), η veto |M(π⁺π⁻π⁰)-m_η|>24.5 MeV/c², M(π⁺π⁻)>0.3 GeV/c²")
alg.note(:signal_extraction, "hc signal yield from 2D fit to M(π⁰hc) vs M(π⁺π⁻J/ψ); upper limits set via Bayesian method")

alg.execute_on([psi3686_data, psi3686_incMC, exMC_signal])