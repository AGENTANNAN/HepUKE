# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive MC

# Decay card: J/psi -> Sigma+ anti-Sigma-,
#             Sigma+ -> p + (massless invisible, written as "nu"),
#             anti-Sigma- -> anti-p pi0, pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Sigma+ anti-Sigma-    PHSP;
    Enddecay

    Decay Sigma+
    1.0000 p+ nu                 PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000 anti-p- pi0           PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma           PHSP;
    Enddecay

    End
DECAYCARD

# 1M-event exclusive MC for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_sigma_p_invisible"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "SigmaInv"
my_alg = Algorithm.new(alg_name)
my_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})

# BOSS-side procedures that cannot be expressed with formal DSL constructs
my_alg
  .note(:pi0_mass_window, "pi0 candidates restricted to the 115-150 MeV/c^2 gamma-gamma invariant-mass window")
  .note(:background_veto, "double-tag competing hypotheses stored as chi2 for ROOT-level veto: 2C fit with the invisible mass set to the pi0 mass, 5C fit with an extra photon (Sigma+ -> p gamma), 6C fit with an extra pi0 (Sigma+ -> p pi0)")

selection = Selection.new
    .select_track {                       # at least two charged tracks
        cos_theta 0.93                    # |cos(theta)| < 0.93
        Vz        10.0                    # |Vz| < 10 cm
        Vr        2.0                     # Vr < 2 cm
        nChrp     ">=1"
        nChrn     ">=1"
    }
    .select_photon {                      # at least two photons
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0            # > 10 deg from any charged track
        energyThreshold_b 0.025           # 25 MeV (barrel)
        energyThreshold_e 0.050           # 50 MeV (endcap)
        nGam              ">=2"
    }
    .pid(method: :probability) {          # proton identification (probability method)
        prob_cut 0.001                    # probability cut 0.001
        identify :proton, against: [:pion, :kaon]   # against pi+ and K+
    }
    .select_isolated_photon {             # isolated photons
        angle_to_prp_track 20.0           # > 20 deg from the proton candidate
        nGam               ">=2"          # (the 10 deg from other tracks is in select_photon)
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit: gamma gamma mass -> m(pi0)
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0     ">=1"
    }
    # Single-tag: combine the pi0 with the anti-proton to form the anti-Sigma-,
    # taking the combination closest to the nominal anti-Sigma- mass
    .secondary_vertex_fit([:prm, :pi0]) {
        build_virtual_particle(:Sigma_bar).by_minimizing_mass_difference
    }
    # Double-tag nominal 2C fit: J/psi -> p anti-p pi0 + invisible
    # (pi0 mass and zero invisible mass constrained); anti-Sigma- mass window
    # |M(anti-p pi0) - M_Sigma| < 15 MeV/c^2 applied here.
    .kinematic_fit([:prp, :prm, :pi0]) {
        nominal
        miss_track_of(:nu)                            # massless invisible recoiling against p anti-p pi0
        constrain_four_momentum
        chi2_cut 200                                  # loose cut; tight cut applied in ROOT
        invariant_mass_of(:prm, :pi0).within(1.18245, 1.21245)
    }
    # Competing 2C hypothesis (invisible mass set to the pi0 mass); chi2 stored for ROOT veto
    .kinematic_fit([:prp, :prm, :pi0]) {
        miss_track_of(:nu)
        constrain_four_momentum
    }
    # 5C hypothesis with an extra photon (veto Sigma+ -> p gamma); chi2 stored for ROOT veto
    .kinematic_fit([:prp, :prm, :pi0, :gamma]) {
        constrain_four_momentum
    }
    # 6C hypothesis with an extra pi0 (veto Sigma+ -> p pi0); chi2 stored for ROOT veto
    .kinematic_fit([:prp, :prm, :pi0, :pi0]) {
        constrain_four_momentum
    }

my_alg.with_decay_card(decay_card_signal).apply(selection)
root_files = my_alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])