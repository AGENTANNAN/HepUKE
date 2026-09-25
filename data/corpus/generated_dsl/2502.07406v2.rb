# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Real data and inclusive MC at 13 energy points from 4.600 to 4.951 GeV.
# Sample-name convention: [BOSS version]_[CMS energy in MeV]
energy_points = %w[
  703_4600 706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840 707_4914 707_4946
]
data_points   = energy_points.map { |name| DatasetManager.real_data.find(name) }     # real data
inc_mc_points = energy_points.map { |name| DatasetManager.inclusive_mc.find(name) } # inclusive MC

# Decay card for the signal process (EvtGen format, phase space)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 K_S0 K_S0 h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive signal MC per energy point (same card, one sample per energy point;
# sample names are auto-suffixed with the corresponding dataset name)
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_K0K0hc_etac"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "K0K0hcEtaC"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.600]})   # CMS energy (GeV); scan runs off this constant
   .set_alias({"std::vector<double>" => "Vdouble"})
   # BOSS-side criteria that have no dedicated DSL method — captured for serialisation
   .note(:ks0_mass_window,
         "K_S0 candidates from the secondary-vertex fit are kept only within "
         "+/-6 MeV/c^2 of the PDG K_S0 mass (0.497614 GeV/c^2); the pair-mass window "
         "has no reserved method inside the secondary_vertex_fit block")
   .note(:ks0_decay_length,
         "K_S0 decay length required to be greater than twice the vertex resolution "
         "(flight-length significance of the pi+pi- vertex), applied right after the "
         "secondary vertex fit; not expressible in the current DSL")

event_selection = Selection.new
event_selection.select_track {              # charged track selection
                 cos_theta 0.93             # |cos(theta)| < 0.93
                 Vz        10.0             # |Vz| < 10 cm
                 Vr        1.0              # Vr < 1 cm in the transverse plane
                 nChrp     "==2"            # exactly two positive tracks
                 nChrn     "==2"            # exactly two negative tracks
               }
               .select_photon {             # photon selection
                 tdc_emc_start     0        # EMC time window [0, 700] ns (50 ns per unit)
                 tdc_emc_end       14
                 energyThreshold_b 0.025    # 25 MeV threshold in the barrel
                 energyThreshold_e 0.050    # 50 MeV threshold in the endcap
                 nGam              ">=1"    # at least one good photon
               }
               .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks treated as pions
               .secondary_vertex_fit([:pip, :pim]) {       # K_S0 from pi+ pi- pairs
                 build_virtual_particle(:K_S0).by_minimizing_mass_difference
                 remove_used_particle_from_candidate_list
               }
               # Partial reconstruction: K_S0 K_S0 gamma are reconstructed while the
               # eta_c is NOT reconstructed — it is identified as the recoil.
               # Decay-card rec IDs: 0 psi(4260) | 1 K_S0 | 2 pi+ | 3 pi- | 4 K_S0 |
               #                     5 pi+ | 6 pi- | 7 h_c | 8 gamma(h_c) | 9 eta_c |
               #                    10 gamma(eta_c) | 11 gamma(eta_c)
               .partial_miss([9]) {
                 # E1 photon: the gamma giving M(gamma K_S0 K_S0) closest to the eta_c mass
                 best_combination_by_mass :eta_c, 2.9836
                 # Recoil mass against the reconstructed K_S0 K_S0 gamma system (GeV/c^2)
                 require_recoil_mass 2.94, 3.06
               }

# Generate the complete algorithm for the process defined in the decay card
alg.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC and the per-energy-point signal MC samples
root_files = alg.execute_on(data_points + inc_mc_points + exMCs_signal)