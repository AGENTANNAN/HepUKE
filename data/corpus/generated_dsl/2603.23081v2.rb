# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# J/psi real data (10087M events, 3.097 GeV) and its corresponding inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for J/psi -> gamma eta pi0, eta -> gamma gamma, pi0 -> gamma gamma (phase space)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta pi0 PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive phase-space MC sample: 1,000,000 events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_eta_pi0_exmc"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name     = "JpsiGammaEtaPi0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # fixed CMS energy 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# Build the event selection chain
event_selection = Selection.new
  .select_track {            # Require zero charged tracks
    nChrp "==0"              # no positive tracks
    nChrn "==0"              # no negative tracks
    nNet  "==0"              # net charge zero
  }
  .select_photon {           # At least five good photons
    tdc_emc_start     0      # TDC start time
    tdc_emc_end       14     # TDC end time
    angle_to_track    10.0   # min angle to nearest charged track (degrees)
    energyThreshold_b 0.025  # EMC barrel energy threshold (GeV)
    energyThreshold_e 0.050  # EMC endcap energy threshold (GeV)
    nGam ">=5"               # at least 5 photons
  }
  # 6C kinematic fit over the five photons: four-momentum conservation (4C)
  # plus nominal pi0 (1C) and eta (1C) mass constraints on the diphoton pairs;
  # the combination with the smallest chi2 is chosen automatically.
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200             # loose BOSS cut; paper's tighter chi2 < 25 applied later in ROOT
  }

# Inexpressible BOSS-side suppression criteria captured for downstream systematics handling
my_algorithm
  .note(:background_veto, "min(Delta_pi0^2) = min over all four-photon combinations of " \
    "[(m(gg)-m_pi0)^2 + (m(g'g'')-m_pi0)^2] required > 0.05 GeV^2/c^4; " \
    "omega -> gamma pi0 vetoes |m(pi0 gamma_eta)-m_omega| < 65 MeV/c^2 and " \
    "|m(gamma gamma_eta)-m_omega| < 65 MeV/c^2; " \
    "J/psi -> gamma eta eta veto |m(gamma gamma_pi0)-m_eta| < 39 MeV/c^2")
  .note(:best_candidate_selection, "final five-photon assignment selected so that the two " \
    "diphoton masses best match the nominal pi0 (resolution 5 MeV/c^2) and eta " \
    "(resolution 9 MeV/c^2)")

# Generate the complete algorithm for the process defined in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC, and exclusive signal MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])