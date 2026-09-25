# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi (3.097 GeV) real data sample
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC sample

# Decay card for the signal process J/psi -> Lambda anti-Sigma0 eta (+ c.c.)
#   anti-Sigma0 -> gamma anti-Lambda0 ; Lambda0 -> p+ pi- ; anti-Lambda0 -> anti-p- pi+ ; eta -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    0.5000 Lambda0 anti-Sigma0 eta       PHSP;
    0.5000 anti-Lambda0 Sigma0 eta       PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 gamma anti-Lambda0            PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0                 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                        HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                   HypWK;
    Enddecay

    Decay eta
    1.0000 gamma gamma                   PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal decay chain (500k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_LambdaSigma0Eta"
  config.related_dataset = jpsi_data                # Associated real dataset
  config.events          = 500000                   # Number of events to generate
  config.decay_card      = decay_card_signal        # Decay card defining the process
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiLambdaSigma0Eta"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})     # CMS energy of the J/psi
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  .select_track {                 # Charged track selection
    cos_theta   0.93              # |cos(theta)| < 0.93
    Vz          100.0             # |Vz| < 100 cm
    Vr          10.0              # |Vr| < 10 mm
    nChrp       ">=2"             # at least two positive tracks
    nChrn       ">=2"             # at least two negative tracks
    nNet        "==0"             # net charge zero
  }
  .select_photon {                # Photon selection
    tdc_emc_start     0           # TDC start time (700 ns unit)
    tdc_emc_end       14          # TDC end time
    angle_to_track    10.0        # opening angle > 10 deg to nearest charged track
    energyThreshold_b 0.050       # E > 50 MeV in the barrel region
    energyThreshold_e 0.050       # E > 50 MeV in the endcap region
    nGam              ">=3"       # at least three photons
  }
  # No explicit PID is applied: the charged tracks are used directly under both the
  # pion and the proton hypotheses (needed for the Lambda / anti-Lambda vertex fits).
  .assign({:chrgp => :pip, :chrgn => :pim})
  .assign({:chrgp => :prp, :chrgn => :prm})
  .secondary_vertex_fit([:prp, :pim]) {        # Lambda0 -> p+ pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {        # anti-Lambda0 -> anti-p- pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Nominal 5C fit: 4C energy-momentum conservation plus the eta mass constraint on gamma gamma
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma, :gamma]) {
    nominal                                                    # nominal fit (its corrected four-momenta are kept)
    constrain_four_momentum                                    # 4C constraint
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)  # 1C eta mass constraint
    chi2_cut 200                                               # loose chi2 cut; tightened in ROOT
  }
  # Competing hypothesis 1: 4C fit to Lambda anti-Lambda gamma gamma (no chi2_cut / no nominal -> chi2 stored)
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma]) {
    constrain_four_momentum
  }
  # Competing hypothesis 2: 4C fit to Lambda anti-Lambda gamma gamma gamma gamma
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }
  # Competing hypothesis 3: 5C fit with a pi0 mass constraint on the same photon/track assignment
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :gamma, :gamma]) {
    use_track_index_from_nominal_kmfit                        # same candidates as the nominal fit
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])