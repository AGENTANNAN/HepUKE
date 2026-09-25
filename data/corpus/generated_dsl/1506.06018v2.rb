# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# Ten XYZ energy points from 4.190 to 4.420 GeV: real data + inclusive MC
xyz_point_names = %w[703_4190 703_4200 703_4210 703_4220 703_4230
                     703_4246 703_4260 703_4270 703_4280 703_4420]
xyz_data  = xyz_point_names.map { |name| DatasetManager.real_data.find(name) }     # Real data at the ten points
xyz_incMC = xyz_point_names.map { |name| DatasetManager.inclusive_mc.find(name) } # Inclusive MC at the ten points

# Decay card for the signal process (EvtGen format):
#   psi(4260) -> pi0 pi0 J/psi (phase space), J/psi -> e+e- / mu+mu- (50/50), pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi0 pi0 J/psi PHSP;
    Enddecay

    Decay J/psi
    0.5000 e+ e-    PHOTOS VLL;
    0.5000 mu+ mu-  PHOTOS VLL;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive signal MC at each of the ten energy points (same card, same cross section)
exMCs_signal = DatasetManager.create_exclusive_mc_for(xyz_data) do |config|
  config.sample_name   = "exmc_xyz_pi0pi0Jpsi"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "Pi0Pi0JpsiXYZ"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            # ECMS is the representative scan energy; the per-point beam energy is resolved at execution
            .set_constant({"ECMS" => [:double, 4.260]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
    .select_track {                 # Charged track selection
        cos_theta 0.93              # |cos(theta)| < 0.93
        Vz        10.0              # |Vz| < 10 cm
        Vr        1.0               # Vr < 1 cm in the transverse plane
        nChrp     "==1"             # Exactly one positively charged track
        nChrn     "==1"             # Exactly one negatively charged track
        nNet      "==0"             # Net charge zero
    }
    .select_photon {                # Photon selection
        tdc_emc_start     0         # TDC start time
        tdc_emc_end       14        # TDC end time
        angle_to_track    5.0       # Angle to the nearest charged track > 5 degrees
        energyThreshold_b 0.025     # 25 MeV threshold in the EMC barrel
        energyThreshold_e 0.050     # 50 MeV threshold in the EMC endcap
        nGam              ">=4"     # At least four photons
    }
    .pid(method: :probability) {    # PID by the probability method
        prob_cut 0.001              # PID probability > 0.001
        # High-momentum tracks (p > 1.0 GeV) treated as leptons; electron if EMC energy > 0.6 GeV, else muon.
        # Creates the combined lepton lists index_lp / index_lm.
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                        treat_as_electron_if_energy_above: 0.6
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C Kalman fit: reconstruct pi0 from photon pairs
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # Constrain m(gamma gamma) to the nominal pi0 mass
        invariant_mass_of(:gamma, :gamma).within(0.100, 0.160)                # Loose window 100 < M(gamma gamma) < 160 MeV/c^2
        chi2_cut 200                                                          # chi2_1C < 200
        npi0 ">=2"                                                            # At least two pi0 candidates
    }
    .kinematic_fit([:pi0, :pi0, :lp, :lm]) {    # Nominal 4C kinematic fit to the pi0 pi0 l+ l- hypothesis
        nominal                                 # Nominal fit: corrected four-momenta are the ones kept
        constrain_four_momentum                 # 4C energy-momentum constraint
        invariant_mass_of(:lp, :lm).within(2.95, 3.20)  # Pre-fit J/psi mass window 2.95 < M(l+l-) < 3.2 GeV/c^2
        chi2_cut 200                            # chi2_4C < 200
    }
    .kinematic_fit([:pi0, :pi0, :lp, :lm]) {    # 7C fit: 4C + the two pi0 masses + the J/psi mass
        constrain_four_momentum                 # 4C energy-momentum constraint
        invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)  # J/psi mass constraint
        chi2_cut 200                            # chi2_7C < 200
    }

# BOSS-side procedures that cannot be expressed with the current DSL constructs
my_algorithm
  .note(:pid_correction_method, "Refinement of the high-momentum lepton identification: " \
        "electrons are required to have E/p > 0.7 and muons E/p < 0.3; at least one muon " \
        "must have hits in more than six MUC layers.")
  .note(:pi0_pair_selection, "Among all gamma-gamma pairings passing the 1C Kalman fit and the " \
        "loose 100 < M(gamma gamma) < 160 MeV/c^2 window, the pi0 pi0 pair minimising " \
        "chi2_1C + chi2_1C + chi2_4C and sharing no photon is selected. The M(gamma gamma) window " \
        "is tightened to 120-150 MeV/c^2 to reject extra combinations.")
  .note(:seven_constraint_fit, "The 7C fit additionally constrains the mass of each of the two pi0 " \
        "mesons; since both pi0 candidates are already mass-constrained by the preceding 1C Kalman " \
        "fit, only the J/psi mass constraint is declared in the DSL on top of the 4C constraint.")
  .note(:background_veto, "Bhabha / two-photon veto: the event is rejected if either track is " \
        "identified as e+ or e- with |cos(theta)| > 0.5 and the opening angle between the two " \
        "tracks exceeds 175 degrees.")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute the algorithm on real data, inclusive MC and the ten signal MC samples
root_files = my_algorithm.execute_on(xyz_data + xyz_incMC + exMCs_signal)