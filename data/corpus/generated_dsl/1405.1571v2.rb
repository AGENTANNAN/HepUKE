# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
# ψ(3686) data and inclusive MC at 3.686 GeV
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")
# 3.65 GeV continuum data (44.49 pb^-1) and its inclusive MC, used for background subtraction
cont_data  = DatasetManager.real_data.find("709_3650")
cont_incMC = DatasetManager.inclusive_mc.find("709_3650")

# Decay card for the signal process ψ(3686) → ω K+ K-, ω → π+ π- π0, π0 → γγ
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 omega K+ K- PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the peaking background ψ(3686) → γ η_c(2S), η_c(2S) → ω K+ K-
decay_card_etac2S = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma eta_c(2S) PHSP;
    Enddecay

    Decay eta_c(2S)
    1.000 omega K+ K- PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 100k ψ(3686) → ω K+ K- events at 3.686 GeV
exMC_signal_3686 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_omegaKK_3686"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Exclusive MC: 50k ψ(3686) → ω K+ K- events at 3.65 GeV
exMC_signal_3650 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_omegaKK_3650"
  config.related_dataset = cont_data
  config.events          = 50000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Exclusive MC: 50k peaking-background events ψ(3686) → γ η_c(2S)
exMC_etac2S = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_etac2S_omegaKK"
  config.related_dataset = psip_data
  config.events          = 50000
  config.decay_card      = decay_card_etac2S
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Common event selection shared by both energy points
event_selection_common = Selection.new
  .select_track {                                   # charged track selection
    cos_theta 0.93                                  # |cosθ| < 0.93
    Vz        10.0                                  # |Vz| < 10 cm
    Vr        1.0                                   # Vr < 1 cm
    nChrp     "==2"                                 # two positively charged tracks
    nChrn     "==2"                                 # two negatively charged tracks
    nNet      "==0"                                 # zero net charge
  }
  .select_photon {                                  # photon selection
    tdc_emc_start     0                             # EMC hit time 0 ns
    tdc_emc_end       14                            # ... to 700 ns
    angle_to_track    20.0                          # angle to any charged track > 20 degrees
    energyThreshold_b 0.025                         # > 25 MeV in the EMC barrel
    energyThreshold_e 0.050                         # > 50 MeV in the EMC endcap
    nGam              ">=2"                         # at least two photons
  }
  .pid(method: :probability) {                      # PID by the probability method
    prob_cut 0.001                                  # P > 0.001
    identify :kaon, against: [:pion, :proton]       # K+ and K- separation
    identify :pion, against: [:kaon, :proton]       # π+ and π- separation
    nkp   ">=1"                                     # at least one K+
    nkm   ">=1"                                     # at least one K-
    npip  ">=1"                                     # at least one π+
    npim  ">=1"                                     # at least one π-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {         # π0 → γγ mass-constrained (1-C) fit
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"                                      # at least one π0 candidate
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :pi0]) {    # 4C kinematic fit to K+ K- π+ π- π0
    vertex_fit([0, 1, 2, 3])                        # common vertex for the four charged tracks
    invariant_mass_of(:pip, :pim, :pi0).within(0.752, 0.812)  # ω mass window
    constrain_four_momentum
    chi2_cut 200
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {  # 5C fit = 4C + π0 mass constraint
    nominal                                          # final analysis uses this fit
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).within(0.11, 0.15)                  # π0 window
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # π0 mass constraint
    chi2_cut 200
  }

# Algorithm at the ψ(3686) energy point (3.686 GeV)
alg_name_3686 = "OmegaKpKm"
alg_3686 = Algorithm.new(alg_name_3686)
alg_3686.set_header(["#{alg_name_3686}Alg/#{alg_name_3686}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .note(:efficiency_curve, "signal efficiencies are determined with a data-driven BODY3 generator at both energy points (3.686 and 3.65 GeV)")
alg_3686.with_decay_card(decay_card_signal).apply(event_selection_common.dup)
root_files_psip = alg_3686.execute_on([psip_data, psip_incMC, exMC_signal_3686, exMC_etac2S])

# Algorithm at the 3.65 GeV continuum energy point (for background subtraction)
alg_name_3650 = "OmegaKpKm3650"
alg_3650 = Algorithm.new(alg_name_3650)
alg_3650.set_header(["#{alg_name_3650}Alg/#{alg_name_3650}.h"])
        .set_constant({"ECMS" => [:double, 3.650]})
        .note(:efficiency_curve, "signal efficiencies are determined with a data-driven BODY3 generator at both energy points (3.686 and 3.65 GeV)")
alg_3650.with_decay_card(decay_card_signal).apply(event_selection_common.dup)
root_files_cont = alg_3650.execute_on([cont_data, cont_incMC, exMC_signal_3650])