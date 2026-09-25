# frozen_string_literal: true
### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data @ 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # ψ(3686) inclusive MC
cont_data  = DatasetManager.real_data.find("712_3773")      # 3.773 GeV data (2.93 fb^-1) for continuum subtraction
cont_incMC = DatasetManager.inclusive_mc.find("712_3773")   # 3.773 GeV inclusive MC

# Scan points around 3.670-3.710 GeV
scan_data = [
  DatasetManager.real_data.find("704_psip_scan_2"),   # 3670.2 MeV
  DatasetManager.real_data.find("704_psip_scan_3"),   # 3680.1 MeV
  DatasetManager.real_data.find("704_psip_scan_4"),   # 3682.8 MeV
  DatasetManager.real_data.find("704_psip_scan_5"),   # 3684.2 MeV
  DatasetManager.real_data.find("704_psip_scan_6"),   # 3685.3 MeV
  DatasetManager.real_data.find("704_psip_scan_7"),   # 3686.5 MeV
  DatasetManager.real_data.find("704_psip_scan_8"),   # 3691.4 MeV
  DatasetManager.real_data.find("704_psip_scan_9"),   # 3709.8 MeV
]

# Decay card: ψ(3686) -> p p̄ π0, π0 -> γγ
decay_card_ppi0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 p+ anti-p- pi0     PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma        PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: ψ(3686) -> p p̄ η, η -> γγ
decay_card_ppeta = <<~DECAYCARD
    Decay psi(2S)
    1.0000 p+ anti-p- eta     PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma        PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for ψ(3686) -> p p̄ π0 (π0 -> γγ)
exMC_ppi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_ppbar_pi0"
  config.related_dataset = psip_data
  config.events         = 500000
  config.decay_card     = decay_card_ppi0
  config.cross_section  = :default
end

# Exclusive MC for ψ(3686) -> p p̄ η (η -> γγ)
exMC_ppeta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_psip_ppbar_eta"
  config.related_dataset = psip_data
  config.events         = 200000
  config.decay_card     = decay_card_ppeta
  config.cross_section  = :default
end

### Event selection (BOSS) ###
alg_name     = "psipToPpbarPi0Eta"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # c.m. energy 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# Common selection: 1 p, 1 anti-p and >=2 photons
event_selection = Selection.new
event_selection
  .select_track {
    cos_theta 0.93    # |cos(theta)| < 0.93
    Vz        10.0    # |Vz| < 10 cm
    Vr        1.0     # Vr < 1 cm
    nChrp     "==1"   # exactly one positive track
    nChrn     "==1"   # exactly one negative track
    nNet      "==0"   # net charge zero
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14       # TDC 0-700 ns
    angle_to_track    20.0     # > 20 deg from nearest charged track
    energyThreshold_b 0.025    # E > 25 MeV in barrel (|cos(theta)| < 0.80)
    energyThreshold_e 0.050    # E > 50 MeV in endcap (0.86 < |cos(theta)| < 0.92)
    nGam              ">=2"    # at least two photons
  }
  # No PID in the main selection: assign the positive/negative track directly as p / anti-p
  .assign({:chrgp => :prp, :chrgn => :prm})
  # Reconstruct π0 from the photon pair (1C Kalman fit), χ² < 100
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 100
    npi0     ">=1"
  }
  # 5C kinematic fit of p p̄ π0 (4C energy-momentum + π0 mass constraint), χ² < 30
  .kinematic_fit([:prp, :prm, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 30
  }

# BOSS-side procedures that cannot be expressed in the DSL
my_algorithm
  .note(:pid_correction_method, "the main selection applies no PID — the single positive/negative track is assigned directly as p / anti-p; a proton-over-pi/K confidence-level requirement is applied only to the 3.670-3.710 GeV scan data")
  .note(:signal_separation, "the eta channel (psi(3686) -> p pbar eta, eta -> gamma gamma) shares this selection; the pi0/eta separation is performed at the ROOT level using the gamma-gamma invariant-mass window")

# Generate the algorithm for the p p̄ π0 signal process in the decay card
my_algorithm.with_decay_card(decay_card_ppi0).apply(event_selection)

# Execute on real data, inclusive MC, continuum data/MC, signal MC and scan points
root_files = my_algorithm.execute_on([psip_data, psip_incMC, cont_data, cont_incMC,
                                      exMC_ppi0, exMC_ppeta] + scan_data)