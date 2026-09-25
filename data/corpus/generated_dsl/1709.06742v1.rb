# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# psi(3686) real data sample (448.1x10^6 events) and its corresponding inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# --- Decay cards (EvtGen format) for the three chi_cJ hypotheses ---
# chi_c0: psi(3686) -> gamma chi_c0, chi_c0 -> gamma gamma
decay_card_chi_c0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# chi_c1: Landau-Yang-forbidden J=1 hypothesis, generated hypothetically for an upper limit
decay_card_chi_c1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# chi_c2: nominal signal generated in a pure helicity-two configuration
decay_card_chi_c2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 PHSP;
    Enddecay

    Decay chi_c2
    1.000 gamma gamma HELAMP 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC samples: 1.2 million events per chi_cJ decay chain ---
exMC_chi_c0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gammachic0_gammagamma"
  config.related_dataset = psip_data          # match the psi(3686) data conditions
  config.events          = 1_200_000
  config.decay_card      = decay_card_chi_c0
  config.cross_section   = :default
end

exMC_chi_c1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gammachic1_gammagamma"
  config.related_dataset = psip_data
  config.events          = 1_200_000
  config.decay_card      = decay_card_chi_c1
  config.cross_section   = :default
end

exMC_chi_c2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gammachic2_gammagamma"
  config.related_dataset = psip_data
  config.events          = 1_200_000
  config.decay_card      = decay_card_chi_c2
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# All three chi_cJ hypotheses share the identical three-photon final state and the
# same selection chain, so a single Algorithm instance serves all of them (the
# decay card only defines the kinematic variables in the generated header).
alg_name = "GammaChicJToGG"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # psi(3686) center-of-mass energy
            .note(:photon_cos_theta_cut,
                  "each selected photon is required to satisfy |cos(theta)| < 0.75 with
                   respect to the positron beam direction; this polar-angle acceptance is
                   an additional requirement on the photon candidates at BOSS level")
            .note(:chic2_helicity_admixture,
                  "chi_c2 signal MC is generated in a pure helicity-two configuration; a
                   variant with a 2% helicity-zero admixture is produced for systematic
                   studies of the helicity content")

# Build the common event selection chain
event_selection = Selection.new
  .select_track {                  # no charged track in the final state
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==0"              # no positive charged track
      nChrn     "==0"              # no negative charged track
      nNet      "==0"              # net charge zero
  }
  .select_photon {                 # exactly three photons
      tdc_emc_start      0
      tdc_emc_end        14
      angle_to_track     10.0
      energyThreshold_b  0.070     # Egamma > 70 MeV in the barrel region
      energyThreshold_e  0.070     # Egamma > 70 MeV in the endcap region
      nGam               "==3"     # exactly three photon candidates
  }
  # No charged-particle PID is needed (no charged tracks in the event).
  # 4C kinematic fit over the three photons: four-momentum conservation.
  # Loose chi2 < 200 at BOSS level; the published chi2_4C <= 80 is applied downstream in ROOT.
  .kinematic_fit([:gamma, :gamma, :gamma]) {
      nominal                     # nominal fit; its corrected four-momenta are saved
      constrain_four_momentum     # 4C energy-momentum constraint
      chi2_cut 200                # loose BOSS-level chi2 cut
  }

# Generate the algorithm (one decay card drives the header for the shared final state)
my_algorithm.with_decay_card(decay_card_chi_c2).apply(event_selection)

# Execute on the real psi(3686) data, inclusive MC and the three signal exclusive MCs
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_chi_c0, exMC_chi_c1, exMC_chi_c2])