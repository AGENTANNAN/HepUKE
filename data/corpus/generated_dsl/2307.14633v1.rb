# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi(3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC at 3.097 GeV

# Decay card for the signal process J/psi -> e+ e- eta(1405), eta(1405) -> pi0 f0(980),
# f0(980) -> pi+ pi-, pi0 -> gamma gamma (EvtGen format, EvtGen particle names)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000  e+  e-  eta(1405)          PHSP;
    Enddecay

    Decay eta(1405)
    1.0000  pi0  f_0(980)              PHSP;
    Enddecay

    Decay f_0(980)
    1.0000  pi+  pi-                   PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma               PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 500k events for the full decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_jpsi_ee_eta1405"
  config.related_dataset = jpsi_data      # Associated real dataset
  config.events          = 500000         # 500k signal MC events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiEta1405"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # ECMS = 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {              # Charged track selection
                  cos_theta 0.93            # |cos(theta)| < 0.93
                  Vz        10.0            # |Vz| < 10 cm
                  Vr        1.0             # Vr < 1 cm
                  nChrp     "==2"           # Exactly 2 positive tracks
                  nChrn     "==2"           # Exactly 2 negative tracks
                  nNet      "==0"           # Net charge zero
                }
               .select_photon {             # Photon selection
                  tdc_emc_start     0       # EMC timing window begins at 0
                  tdc_emc_end       14      # EMC timing window ends at 14 (i.e. 700 ns)
                  angle_to_track    10.0    # At least 10 degrees from any charged track
                  energyThreshold_b 0.025   # 25 MeV in the barrel
                  energyThreshold_e 0.050   # 50 MeV in the endcap
                  nGam              ">=2"   # At least two photons
                }
               .pid(method: :probability) { # Particle identification (probability method)
                  prob_cut 0.001            # PID probability > 0.001
                  # e/mu are never identified via `identify`; declare high-momentum leptons:
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  nlp "==1"                 # One positive lepton (electron)
                  nlm "==1"                 # One negative lepton (electron)
                }
               .remove([:lp <= :chrgp, :lm <= :chrgn])   # Remove identified leptons from the charged lists
               .assign({:chrgp => :pip, :chrgn => :pim}) # Remaining positives/negatives taken as pi+/pi-
               .kalman_kinematic_fit([:gamma, :gamma]) { # Kalman fit of the two photons to the pi0 mass
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 200              # chi^2 < 200
                  npi0 ">=1"                # At least one pi0 candidate
                }
               .kinematic_fit([:lp, :lm, :pip, :pim, :pi0]) { # 4C kinematic fit to e+ e- pi+ pi- pi0
                  nominal                   # Nominal fit: its corrected four-momenta are used
                  constrain_four_momentum   # 4C energy-momentum constraint
                  chi2_cut 20               # chi^2 < 20
                }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC and the exclusive signal MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])