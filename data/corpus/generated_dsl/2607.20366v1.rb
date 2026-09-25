# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# J/psi data and inclusive MC at sqrt(s) = 3.097 GeV (BOSS 7.0.8)
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for J/psi -> gamma X(2370), X(2370) -> K*(892)0 anti-K0 (+ c.c.),
# K*(892)0 -> K_S0 pi0, anti-K0 -> K_S0, K_S0 -> pi+ pi-, pi0 -> gamma gamma.
# X(2370) is modelled as a ~2.36 GeV pseudoscalar decaying through phase space.
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma X(2370) PHSP;
    Enddecay

    Decay X(2370)
    1.0000 K*0 anti-K0 PHSP;
    Enddecay

    Decay K*0
    1.0000 K_S0 pi0 PHSP;
    Enddecay

    Decay anti-K0
    1.0000 K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal decay chain (500k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_X2370_KSKSpi0"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiGX2370KSKSpi0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
    .select_track {                       # Charged track selection
        cos_theta 0.93                    # |cos(theta)| < 0.93
        Vz        100.0                   # |Vz| < 100 cm
        Vr        10.0                    # Vr < 10 cm in the transverse plane
        nChrp     ">=2"                   # At least two positive tracks
        nChrn     ">=2"                   # At least two negative tracks
        nNet      "==0"                   # Net charge zero
    }
    .select_photon {                      # Photon selection
        tdc_emc_start     0               # TDC start time
        tdc_emc_end       14              # TDC end time
        angle_to_track    10.0            # At least 10 deg from any charged track
        energyThreshold_b 0.025           # Barrel energy threshold 25 MeV
        energyThreshold_e 0.050           # Endcap energy threshold 50 MeV
        nGam              ">=3"           # At least three photons
    }
    .assign({:chrgp => :pip, :chrgn => :pim})   # All positives as pi+, all negatives as pi-
    .secondary_vertex_fit([:pip, :pim]) {       # First K_S0 -> pi+ pi- (mass-difference optimized)
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:pip, :pim]) {       # Second K_S0 -> pi+ pi- (mass-difference optimized)
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # Reconstruct pi0 from gamma gamma
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25                             # chi^2 < 25 for the mass-constrained fit
        npi0     ">=1"                          # At least one pi0 candidate
    }
    # final 4C kinematic fit on gamma K_S0 K_S0 pi0
    .kinematic_fit([:gamma, :K_S0, :K_S0, :pi0]) {
        nominal                                 # nominal fit: corrected four-momenta are saved
        constrain_four_momentum                 # 4C energy-momentum constraint
        chi2_cut 200                            # loose chi^2 cut; tight cut applied in ROOT
    }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])