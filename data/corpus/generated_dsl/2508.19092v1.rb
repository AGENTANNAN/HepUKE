# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset preparation ###
psip_data     = DatasetManager.real_data.find("709_3686")     # ψ(2S) real data (2.712×10^9 ψ(2S))
psip_incMC    = DatasetManager.inclusive_mc.find("709_3686")  # Corresponding ψ(2S) inclusive MC
psip3770_data = DatasetManager.real_data.find("712_3773")     # ψ(3770) data used as continuum background

# Decay card for the signal process ψ(2S) → ω η η,
# with ω → π+π−π0 (Dalitz), π0 → γγ, η → γγ (both η), giving π+π−6γ
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 omega eta eta    PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0       OMEGA_DALITZ;
    Enddecay

    Decay eta
    1.000 gamma gamma       PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma       PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC for ψ(2S) → ωηη (200k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_omegaetaeta"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsipToOmegaEtaEta"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track {          # Charged track selection
    cos_theta 0.93                      # |cosθ| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     "==1"                     # Exactly one π+
    nChrn     "==1"                     # Exactly one π−
    nNet      "==0"                     # Net charge zero
  }
  .select_photon {                      # Photon selection
    tdc_emc_start     0                 # EMC timing window start
    tdc_emc_end       14                # EMC timing window end
    energyThreshold_b 0.025             # Barrel energy threshold (25 MeV)
    energyThreshold_e 0.050             # Endcap energy threshold (50 MeV)
    angle_to_track    10.0              # Min angle to any charged track (degrees)
    nGam              ">=6"             # At least six photons (π+π−6γ final state)
  }
  .pid(method: :probability) {          # PID: identify pions against kaons/protons
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit to π0 from γγ
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit to η from γγ
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta ">=2"
  }
  # Nominal 4C kinematic fit to π+π−π0ηη
  .kinematic_fit([:pip, :pim, :pi0, :eta, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 40
  }
  # Competing hypothesis: π+π−7γ (no chi2 cut, no nominal -> stores χ² for ROOT-level veto)
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }
  # Competing hypothesis: π+π−8γ
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }

my_algorithm
  .note(:mass_window_selection, "π0/η candidates built from γγ: require ≥1 π0 with |M_γγ − m_π0| < 20 MeV and ≥2 η with |M_γγ − m_η| < 25 MeV; the combination that minimizes the combined π0–ηη mass-deviation χ² is selected before the 4C fit")
  .note(:background_veto_two_pi0, "two-π0 veto: reject events in which any γγ pair not assigned to the reconstructed π0 has |M_γγ − m_π0| < 0.03 GeV")
  .note(:background_veto_psip_pipi_jpsi, "ψ(2S) → ππJ/ψ veto: reject events with |M(π+π−)_recoil − m_J/ψ| < 0.02 GeV")
  .note(:background_veto_psip_x_jpsi_omega_eta, "ψ(2S) → XJ/ψ, J/ψ → ωη veto: reject events with M(ωη) < 3.0 GeV")
  .note(:background_veto_alt_hypotheses, "alternative-hypothesis veto: reject events for which a 4C fit under the π0π0η, π0π0π0 or ηηη assignment yields a smaller χ² than the nominal π+π−π0ηη fit (χ² of these reassignments to be stored for the ROOT-level comparison)")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC, continuum (ψ(3770)) and exclusive signal MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, psip3770_data, exMC_signal])