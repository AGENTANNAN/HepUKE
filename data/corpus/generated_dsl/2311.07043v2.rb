# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data = DatasetManager.real_data.find("708_3097")        # J/psi real data at 3.097 GeV (BOSS 7.0.8)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive MC (~10B J/psi events)

# Decay card for the signal process J/psi -> phi pi0 eta (EvtGen format)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 phi pi0 eta        PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-              VSS;
    Enddecay

    Decay pi0
    1.0000 gamma gamma        PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma        PHSP;
    Enddecay

    End
DECAYCARD

# Create the exclusive signal MC sample of the same decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_phi_pi0_eta"
  config.related_dataset = jpsi_data        # Associated real dataset
  config.events          = 100000           # Number of events to generate
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiPhiPi0Eta"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy = J/psi mass (GeV)
         .set_alias({"std::vector<double>" => "Vdouble"})
         # The pi0pi0 / etaeta competing-hypothesis vetoes are stored by the non-nominal
         # kinematic fits below; the actual windows are applied in the ROOT analysis.
         .note(:background_veto, "J/psi -> K+K-pi0pi0 and J/psi -> K+K-etaeta backgrounds " \
              "suppressed by chi2(pi0pi0) > 90 and chi2(etaeta) > 8; the competing-hypothesis " \
              "chi2 values are stored by the non-nominal kinematic fits and the veto is applied " \
              "in the ROOT analysis.")

# Build the event selection chain (tracks -> photons -> PID -> pi0/eta -> kinematic fits)
event_selection = Selection.new
  .select_track {                        # Charged track selection
    cos_theta 0.93                       # |cos(theta)| < 0.93
    Vz        10.0                       # |Vz| < 10 cm
    Vr        1.0                        # Vxy < 1 cm
    nChrp     "==1"                      # exactly one positive track
    nChrn     "==1"                      # exactly one negative track
    nNet      "==0"                      # net charge zero
  }
  .select_photon {                       # Photon selection
    tdc_emc_start     0                  # TDC start time
    tdc_emc_end       14                 # TDC end time
    angle_to_track    10.0               # min angle to nearest charged track (deg)
    energyThreshold_b 0.025              # E > 25 MeV (barrel)
    energyThreshold_e 0.050              # E > 50 MeV (endcap)
    nGam              ">=4"              # at least four photons
  }
  .pid(method: :probability) {           # Kaon identification
    prob_cut 0.001                       # PID probability > 0.001
    identify :kaon, against: [:pion]     # identify K+ and K- against pions
    nkp "==1"                            # exactly one K+
    nkm "==1"                            # exactly one K-
  }
  # 1C mass-constrained fit: reconstruct pi0 from a photon pair (gamma gamma -> pi0)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"                           # at least one pi0 candidate
  }
  # 1C mass-constrained fit: reconstruct eta from a photon pair (gamma gamma -> eta)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"                           # at least one eta candidate
  }
  # Nominal kinematic fit: 4C energy-momentum constraint together with the pi0 and eta mass
  # constraints (6C in total). The best (minimum chi2) gamma-gamma pairing into pi0/eta is
  # chosen automatically when the event contains more than four photons.
  .kinematic_fit([:kp, :km, :pi0, :eta]) {
    nominal
    constrain_four_momentum             # 4C energy-momentum conservation
    invariant_mass_of(:kp, :km).within(1.01, 1.03)   # phi mass window M(K+K-)
    chi2_cut 45                         # chi2 < 45
  }
  # Competing-hypothesis fit: K+K-pi0pi0 (no chi2_cut, no nominal -> stores chi2 for the veto)
  .kinematic_fit([:kp, :km, :pi0, :pi0]) {
    constrain_four_momentum
  }
  # Competing-hypothesis fit: K+K-etaeta (no chi2_cut, no nominal -> stores chi2 for the veto)
  .kinematic_fit([:kp, :km, :eta, :eta]) {
    constrain_four_momentum
  }

# Generate the algorithm for the signal process and apply the selection
algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on the real data, inclusive MC, and the exclusive signal MC
root_files = algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])