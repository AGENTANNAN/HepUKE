# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive MC sample

# Decay card for the signal process  J/psi -> gamma gamma phi, phi -> K+ K-
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma gamma phi PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal mode (500k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_ggphi"       # Signal exclusive MC sample name
  config.related_dataset = jpsi_data               # Associated real dataset (J/psi)
  config.events          = 500000                  # 500k generated events
  config.decay_card      = decay_card_signal       # Decay card defining the signal
  config.cross_section   = :default                # Use default cross-section
end

### Event selection (BOSS) ###
alg_name = "JpsiToGGPhi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])   # Header file for the algorithm
            .set_constant({"ECMS" => [:double, 3.097]})     # ECMS = 3.097 GeV (J/psi)
            .set_alias({"std::vector<double>" => "Vdouble"})

# Full event selection chain, all functions chained for readability
event_selection = Selection.new
  .select_track {                       # Charged track selection
      cos_theta 0.93                    # |cos(theta)| < 0.93
      Vz        10.0                    # |Vz| < 10 cm along beam direction
      Vr        1.0                     # Vr < 1 cm in transverse plane
      nChrp     "==1"                   # Exactly one positively charged track
      nChrn     "==1"                   # Exactly one negatively charged track
      nNet      "==0"                   # Net charge zero
  }
  .select_photon {                      # Photon selection
      tdc_emc_start     0               # EMC TDC start time
      tdc_emc_end       14              # EMC TDC end time
      angle_to_track    10.0            # Min angle to nearest charged track (deg)
      energyThreshold_b 0.025           # E > 25 MeV in the barrel region
      energyThreshold_e 0.050           # E > 50 MeV in the endcap region
      nGam              ">=2"           # At least two photons
  }
  .pid(method: :probability) {          # Kaon identification (probability method)
      prob_cut 0.001                    # PID probability > 0.001
      identify :kaon, against: [:pion, :proton]  # K+ and K- vs pi/p (charge-conjugation shorthand)
      nkp ">=1"                         # At least one K+
      nkm ">=1"                         # At least one K-
  }
  # Nominal 4C kinematic fit to gamma gamma K+ K-
  .kinematic_fit([:gamma, :gamma, :kp, :km]) {
      nominal                           # Mark as the nominal fit
      constrain_four_momentum           # 4C energy-momentum constraint
      # pi0 / eta / eta' rejection on the two-photon invariant mass
      invariant_mass_of(:gamma, :gamma).out_of(0.105, 0.165)   # reject pi0 (|M(g g) - m(pi0)| > 0.03)
      invariant_mass_of(:gamma, :gamma).out_of(0.50, 0.58)     # reject eta (M(g g) < 0.50 or > 0.58)
      invariant_mass_of(:gamma, :gamma).out_of(0.928, 0.988)   # reject eta' (|M(g g) - m(eta')| > 0.03)
      chi2_cut 40                       # chi^2 < 40
  }
  # Competing 4C fit to 3 gamma K+ K- (no chi2_cut, no nominal -> stores chi2_3gamma for ROOT-level veto)
  .kinematic_fit([:gamma, :gamma, :gamma, :kp, :km]) {
      constrain_four_momentum
  }
  # Competing 4C fit to 4 gamma K+ K- (no chi2_cut, no nominal -> stores chi2_4gamma for ROOT-level veto)
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :kp, :km]) {
      constrain_four_momentum
  }

# Generate the complete algorithm for the signal process defined in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC, and exclusive signal MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])