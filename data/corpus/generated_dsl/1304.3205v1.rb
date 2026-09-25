### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/ψ real data (225.3 M events)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive J/ψ MC

# Decay card for the LFV signal J/ψ → e±μ∓ (EvtGen format, EvtGen particle names)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 e+ mu- PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC for the e+ μ- final state (100k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_jpsi_emu"
  config.related_dataset = jpsi_data                 # associated real dataset J/ψ(3097)
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiToEMu"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # J/ψ center-of-mass energy (GeV)
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {          # charged track selection
                  cos_theta 0.8         # |cosθ| < 0.8
                  Vz        5.0         # |Vz| < 5 cm
                  Vr        1.0         # Vr < 1 cm
                  nChrp     "==1"       # exactly one positive track
                  nChrn     "==1"       # exactly one negative track
                  nNet      "==0"       # net charge zero
                }
               .select_photon {         # photon veto: require zero good photon candidates
                  tdc_emc_start      0
                  tdc_emc_end        14
                  energyThreshold_b  0.015   # 15 MeV, EMC barrel
                  energyThreshold_e  0.015   # 15 MeV, EMC endcap
                  angle_to_track     20.0    # > 20° from any extrapolated track
                  nGam               "==0"   # zero good photons
                }
               .pid(method: :probability) {      # lepton identification
                  # p > 1.0 GeV → lepton; EMC energy > 0.6 GeV → electron, else muon
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  nlp "==1"          # exactly one l+
                  nlm "==1"          # exactly one l-
                }
               .kinematic_fit([:lp, :lm]) {      # 4C fit of the two leptons to the J/ψ four-momentum
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }

# BOSS-side selection criteria with no corresponding DSL construct:
my_algorithm
  .note(:cosmic_rejection, "reject cosmic rays: require the TOF difference between the two charged tracks < 1.0 ns")
  .note(:event_shape,      "require acollinearity < 0.9 deg and acoplanarity < 1.4 deg between the two charged tracks")
  .note(:signal_region,    "signal region: 0.93 <= E_vis/sqrt(s) <= 1.10 and |sum p|/sqrt(s) <= 0.10")

# Generate the complete algorithm for the process in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC, and exclusive signal MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])