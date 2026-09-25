# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# sqrt(s) = 4.178 GeV  ->  BOSS 7.0.3 real data sample "703_4180" and its inclusive MC
data_4178  = DatasetManager.real_data.find("703_4180")
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")

# Decay card: e+e- -> Ds*+ Ds- ; Ds*+ -> Ds+ gamma ; Ds+ -> pi+ pi- pi+ (signal) ; Ds- -> K+ K- pi- (tag)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.0000 D_s+ gamma PHSP;
    Enddecay

    Decay D_s+
    1.0000 pi+ pi- pi+ PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the full decay chain (500k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "dsstar_ds_dspipipip"
  config.related_dataset = data_4178
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "DsStarDsPipipi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.178]})       # sqrt(s) = 4.178 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# BOSS-side procedures that have no dedicated DSL construct
my_algorithm.note(:background_veto,
                  "reject pi+/pi- originating from K_S0 -> pi+ pi- : the pi+pi- pair invariant mass must
                   fall within +/-12 MeV of the nominal K_S0 mass AND the pi+pi- secondary-vertex decay
                   length must exceed 2x the vertex resolution")
            .note(:background_veto_pi0,
                  "reject photons originating from pi0 -> gamma gamma : the gamma-gamma invariant mass
                   must lie outside [0.125, 0.145] GeV")
            .note(:tag_dsstar,
                  "tag Ds*+ -> Ds+ gamma against the reconstructed Ds- (K+K-pi-) via either the direct
                   recoil mass |M_rec - m(Ds*+)| < 0.02 GeV or the indirect mass difference
                   Delta M = M(Ds+ gamma) - M(Ds+) in [0.135, 0.150] GeV")
            .note(:tmva_suppression,
                  "TMVA neural-network classifier applied for background suppression; the optimal
                   classifier-output cut is determined outside the formal BOSS selection")

# Build the event selection chain up to (and including) the nominal kinematic fit
event_selection = Selection.new
event_selection.select_track {              # charged track quality cuts
                  cos_theta   0.93          # |cos(theta)| < 0.93
                  Vz          10.0          # |Vz| < 10 cm
                  Vr          1.0           # Vr < 1 cm
                  nTot        "==3"         # exactly three charged tracks
                  nChrp       ">=1"         # at least one positive track
                  nChrn       ">=1"         # at least one negative track
                }
               .select_photon {             # photon selection
                  tdc_emc_start      0      # EMC TDC start time
                  tdc_emc_end        14     # EMC TDC end time
                  angle_to_track     10.0   # > 10 degrees to nearest charged track
                  energyThreshold_b  0.025  # 25 MeV in the EMC barrel
                  energyThreshold_e  0.050  # 50 MeV in the EMC endcap
                  nGam               ">=1"  # at least one photon (for Ds*+ -> Ds+ gamma)
                }
               .pid(method: :probability) { # probability-based particle identification
                  prob_cut 0.001            # PID probability > 0.001
                  identify :pion, against: [:kaon, :proton]   # pi+ and pi- vs K / p
                  npip ">=1"                # at least one pi+
                  npim ">=1"                # at least one pi-
                }
               .kinematic_fit([:pip, :pim, :pip]) {          # nominal 4C fit to the signal-side pi+pi-pi+
                  nominal                                    # mark as the nominal fit
                  constrain_four_momentum                    # 4C energy-momentum constraint
                  invariant_mass_of(:pip, :pim, :pip)
                    .constrain_to_nominal_mass_of(:"D_s+")   # Ds+ mass constraint
                  chi2_cut 200                               # loose chi2; optimal tight cut applied in ROOT
                }

# Generate the complete algorithm for the process in the decay card and run it
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([data_4178, incMC_4178, exMC_signal])