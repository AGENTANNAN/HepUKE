# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psip_data      = DatasetManager.real_data.find("709_3686")     # psi(2S) (3.686 GeV) real data
psip_incMC     = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC
continuum_data = DatasetManager.real_data.find("709_3650")     # 3.650 GeV off-resonance data (401 pb^-1) for continuum estimate

# Decay card for the signal process: psi(2S) -> gamma K_S0 anti-K_S0, both K_S0 -> pi+ pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma K_S0 anti-K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay anti-K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Signal exclusive MC: 500k events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_KS0KS0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "GammaKS0KS0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {           # Charged track selection
    cos_theta 0.93                        # |cos(theta)| < 0.93
    Vz        20.0                        # |Vz| < 20 cm
    Vr        1.0                         # Vr < 1 cm in the transverse plane
    nChrp     "==2"                       # exactly two positive tracks
    nChrn     "==2"                       # exactly two negative tracks
    nNet      "==0"                       # net charge zero
}
.select_photon {                          # Photon selection
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025               # 25 MeV in the barrel
    energyThreshold_e 0.050               # 50 MeV in the endcap
    angle_to_track    10.0                # at least 10 degrees from any charged track
    nGam              ">=1"               # at least one photon
}
.secondary_vertex_fit([:pip, :pim]) {     # First K_S0 from a pi+pi- pair (minimize mass difference)
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
}
.secondary_vertex_fit([:pip, :pim]) {     # Second K_S0 from the remaining pi+pi- pair
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
}
.kinematic_fit([:gamma, :K_S0, :K_S0]) {  # 4C fit to gamma pi+pi- pi+pi- (== gamma K_S0 anti-K_S0)
    nominal                               # nominal fit: corrected four-momenta are saved
    constrain_four_momentum               # constrain total four-momentum to the CMS energy
    chi2_cut 200                          # loose cut in BOSS; optimal tight cut applied in ROOT
}

# BOSS-side procedures that have no dedicated DSL construct, preserved as notes
my_algorithm
  .note(:ks0_selection, "Require exactly two K_S0 per event; each K_S0, built from a pi+pi- pair via the secondary-vertex fit, must satisfy |M(pi+pi-) - M(K_S0)| < 12 MeV/c^2 (2.5 sigma) and have a decay length exceeding 2 sigma of the vertex-fit resolution.")
  .note(:background_veto, "Events with M(K_S0 K_S0) < 2.8 GeV/c^2 are vetoed to suppress chi_c0 / chi_c2 backgrounds.")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC, off-resonance continuum data, and signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, continuum_data, exMC_signal])