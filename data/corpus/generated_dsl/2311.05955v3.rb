# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# J/psi (3.097 GeV) real data sample and the corresponding inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal chain  J/psi -> pbar Sigma+ K_S0
#   Sigma+ -> p pi0, pi0 -> gamma gamma, K_S0 -> pi+ pi-
# (the charge-conjugate mode J/psi -> p Sigma- K_S0 has the identical final state
#  and shares the very same selection chain; the two yields are separated later
#  at the ROOT stage from the proton / antiproton charge)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 anti-p- Sigma+ K_S0    PHSP;
    Enddecay

    Decay Sigma+
    1.0000 p+ pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-    PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal (no explicit event count was requested)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_pbarsigmaplus_ks0"
  config.related_dataset = jpsi_data
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiPbarSigmaKShort"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # 3.097 GeV

# Build the event-selection chain (shared by both charge-conjugate modes)
event_selection = Selection.new
event_selection
  .select_track {
      cos_theta 0.93   # |cos(theta)| < 0.93
      Vz        20.0   # |Vz| < 20 cm  (applied to the K_S0 daughters)
      nChrp     "==2"  # exactly two positive tracks
      nChrn     "==2"  # exactly two negative tracks
      nNet      "==0"  # net charge zero
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=2"    # at least two photons
  }
  .pid(method: :probability) {
      prob_cut 0                              # loose probability cut (value 0)
      identify :proton, against: [:kaon, :pion]  # identify proton AND antiproton at once
      nprp ">=1"                              # one proton
      nprm ">=1"                              # one antiproton
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])   # take the (anti)protons out of the charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks (K_S0 daughters) assumed pions, no pion PID
  .secondary_vertex_fit([:pip, :pim]) {       # K_S0 -> pi+ pi-
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma (1-C mass constraint)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"
  }
  .kinematic_fit([:prm, :prp, :pi0, :K_S0]) { # 4C fit to pbar Sigma+ K_S0
      nominal
      constrain_four_momentum                 # four-momentum conservation
      invariant_mass_of(:prp, :pi0).constrain_to_nominal_mass_of(:"Sigma+")  # Sigma+ mass constraint
      chi2_cut 200
  }

# BOSS-side procedure that has no direct DSL construct
my_algorithm
  .note(:secondary_vertex_selection,
        "K_S0 secondary-vertex fit requires chi^2 < 100 and decay-length " \
        "significance > 2.0; the Vz < 20 cm cut on the K_S0 daughters is applied " \
        "in the charged-track selection.")

# Generate the algorithm for the decay card and apply the shared selection
my_algorithm
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])