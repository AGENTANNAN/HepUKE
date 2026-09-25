# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data = DatasetManager.real_data.find("708_3097")       # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

# Decay card for the signal process J/psi -> e+ e- pi0, pi0 -> gamma gamma (phase space)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 e+ e- pi0        PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma      PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 500k events, EvtGen phase space
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_ee_pi0"
  config.related_dataset = jpsi_data                     # associated real dataset
  config.events         = 500000                         # 500k events
  config.decay_card     = decay_card_signal
  config.cross_section  = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiToEEPi0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])   # header for generated kinematic variables
            .set_constant({"ECMS" => [:double, 3.097]})     # J/psi center-of-mass energy

# Build the event selection chain up to the final (nominal) 4C kinematic fit
event_selection = Selection.new
event_selection.select_track {                    # charged track selection
                  cos_theta 0.93                  # |cos(theta)| < 0.93
                  Vz        10.0                  # |Vz| < 10 cm
                  Vr        1.0                   # Vr < 1 cm
                  nChrp     "==1"                 # exactly one positive track
                  nChrn     "==1"                 # exactly one negative track
                  nNet      "==0"                 # net charge zero
                }
               .select_photon {                   # photon selection
                  tdc_emc_start     0             # EMC TDC start (0)
                  tdc_emc_end       14            # EMC TDC end (14)
                  energyThreshold_b 0.025         # barrel energy threshold 25 MeV
                  energyThreshold_e 0.050         # endcap energy threshold 50 MeV
                  nGam   ">=2"                    # at least two photons
                }
               .pid(method: :probability) {       # particle identification
                  # high-momentum tracks (p > 1.0 GeV/c) treated as leptons;
                  # lepton with EMC eraw > 0.6 GeV identified as electron
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  nlp   "==1"                     # exactly one e+
                  nlm   "==1"                     # exactly one e-
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {   # 1-C fit to reconstruct pi0 -> gamma gamma
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25                     # chi2 < 25 for the mass-constrained fit
                  npi0   ">=1"                    # at least one pi0 candidate
                }
               .kinematic_fit([:lp, :lm, :pi0]) { # nominal 4C fit over e+ e- pi0
                  nominal                         # this is the nominal fit (corrected 4-momenta saved)
                  vertex_fit([0, 1])              # vertex fit constraining the e+ e- pair to a common vertex
                  constrain_four_momentum         # 4C energy-momentum constraint
                  chi2_cut 200                    # loose chi2 cut; tight (chi2<100) applied later in ROOT
                }

# Additional background-suppression criteria that have no dedicated DSL construct
my_algorithm
  .note(:pid_correction_method, "electron/pion separation: require E/p > 0.8 for electron candidates with track momentum p > 0.25 GeV/c; applied at BOSS selection level")
  .note(:background_veto, "pi0 signal region: require the diphoton invariant mass m(gamma gamma) in [0.09, 0.18] GeV/c^2")
  .note(:background_veto, "photon-conversion veto: reject photons with conversion transverse distance delta_xy < 2 cm")
  .note(:background_veto, "two-photon-process veto: require cos(theta)(e+) < 0.8 and cos(theta)(e-) > -0.8")
  .note(:background_veto, "radiative Bhabha veto: require electron/positron momentum p(e+/-) < 1.45 GeV/c")
  .note(:efficiency_curve, "require the lower-energy photon in the pi0 candidate to satisfy E_gamma > 0.14 GeV")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute the algorithm on real data, inclusive MC and signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])