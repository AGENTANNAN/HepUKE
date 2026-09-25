# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# 4.600 GeV data (BOSS 703) and the corresponding inclusive MC
data_4600  = DatasetManager.real_data.find("703_4600")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

# Decay card for e+e- -> Lambda_c+ anti-Lambda_c- (EvtGen format)
#   Lambda_c+      -> p+ K_S0 eta
#   K_S0           -> pi+ pi-
#   eta            -> gamma gamma
#   anti-Lambda_c- -> anti-p- pi+   (charge-conjugate / second Lambda_c mode)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0 Lambda_c+ anti-Lambda_c-   PHSP;
    Enddecay

    Decay Lambda_c+
    1.0 p+ K_S0 eta                PHSP;
    Enddecay

    Decay K_S0
    1.0 pi+ pi-                    PHSP;
    Enddecay

    Decay eta
    1.0 gamma gamma                PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0 anti-p- pi+                PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal decay chain: 1,000,000 events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4600_lambdac_pkseta"
  config.related_dataset = data_4600
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "LambdaCToPKSeta"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.600]})

event_selection = Selection.new
event_selection
  .select_track {                       # Charged track selection
      cos_theta 0.93                    # |cos(theta)| < 0.93
      Vz        10.0                    # |Vz| < 10 cm
      Vr        1.0                     # Vr < 1 cm in the transverse plane
      nChrp     ">=1"                   # at least one positive track
      nChrn     ">=3"                   # at least three negative tracks
      nNet      "==0"                   # net charge zero
  }
  .select_photon {                      # Photon selection
      tdc_emc_start     0               # EMC TDC start
      tdc_emc_end       14              # EMC TDC end (TDC 0-700 window)
      energyThreshold_b 0.025           # 25 MeV threshold in the barrel
      energyThreshold_e 0.050           # 50 MeV threshold in the endcap
      nGam              ">=2"           # at least two photons (eta -> gamma gamma)
  }
  .pid(method: :probability) {          # Particle identification
      prob_cut 0.001                    # PID probability > 0.001
      identify :proton, against: [:kaon, :pion]  # identify p+ / anti-p-
      nprp     ">=1"                    # at least one proton
  }
  .remove([:prp <= :chrgp])             # remove the identified proton from the positive tracks
  .assign({:chrgp => :pip, :chrgn => :pim})  # remaining positives -> pi+, negatives -> pi-
  .secondary_vertex_fit([:pip, :pim]) { # Reconstruct K_S0 from a pi+ pi- pair
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      invariant_mass_of(:pip, :pim).within(0.487, 0.511)  # K_S0 mass window 0.487-0.511 GeV
      remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {  # Reconstruct eta from two photons
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 25                             # chi2 < 25 for the 1-C mass constraint
      neta     ">=1"                          # at least one eta candidate
  }
  .kinematic_fit([:prp, :K_S0, :eta]) {       # Final 4C fit on the p K_S0 eta system
      nominal                    # nominal fit (corrected four-momenta are used)
      constrain_four_momentum    # constrain total four-momentum to the CMS energy
      chi2_cut 200               # loose chi2 cut; tight cut applied later in ROOT
  }

# Attach the decay card and render the selection into the BOSS algorithm
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([data_4600, incMC_4600, exMC_signal])