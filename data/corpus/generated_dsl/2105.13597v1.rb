### Dataset preparation ###
# 15 R-scan center-of-mass energies (2.000-3.080 GeV) for the
# e+e- -> K_S^0 K_L^0 Born cross-section measurement.
rscan_energies = [2000, 2050, 2100, 2150, 2175, 2200, 2396, 2500,
                  2700, 2800, 2900, 2981, 3000, 3020, 3080]

rscan_data  = rscan_energies.map { |e| DatasetManager.real_data.find("713_#{e}") }     # real data per point
rscan_incMC = rscan_energies.map { |e| DatasetManager.inclusive_mc.find("713_#{e}") } # inclusive MC per point

# ConExc decay card: continuum e+e- -> K_S^0 K_L^0 (ConExc mode 46).
# No KKMC top mother and no `Particle vpho` line: the DSL detects the literal
# `ConExc` token, switches to the no-KKMC template and injects
# `Particle vpho <ECMS> 0.0` at every energy point of the scan.
decay_card_conexc = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 46;
    Enddecay
    End
DECAYCARD

# 100k-event exclusive MC per scan point (same card, different energies)
exMC_signal = DatasetManager.create_exclusive_mc_for(rscan_data) do |config|
  config.sample_name   = "exmc_kskl_scan"
  config.events        = 100_000
  config.decay_card    = decay_card_conexc
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "KsKlScan"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 2.000]})

event_selection = Selection.new
event_selection
  .select_track {                              # charged-track selection
      cos_theta 0.93                           # |cos(theta)| < 0.93
      Vz        20.0                           # |Vz| < 20 cm (no Vr cut)
      nChrp     "==1"                          # exactly one positive track
      nChrn     "==1"                          # exactly one negative track
      nNet      "==0"                          # net charge zero
  }
  .assign({:chrgp => :pip, :chrgn => :pim})    # positive -> pi+, negative -> pi- (no PID)
  # Reconstruct K_S^0 -> pi+ pi- with a secondary-vertex fit
  .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference  # closest to K_S^0 mass
      remove_used_particle_from_candidate_list
  }
  # 4C kinematic fit with K_L^0 treated as the missing particle
  .kinematic_fit([:K_S0]) {
      nominal                     # this is the nominal fit -> corrected 4-momenta are saved
      miss_track_of :K_L0         # K_L^0 undetected
      constrain_four_momentum     # constrain the visible system to the CMS four-momentum
      chi2_cut 200                # loose chi^2 cut (tight cut optimised in the ROOT step)
  }

# BOSS-side procedures with no DSL primitive
my_algorithm
  .note(:electron_veto,
        "Electron veto: require E/p < 0.8 for every charged track; the DSL has "
        "no E/p primitive (ep_ratio_of is unimplemented), so the cut is applied "
        "in BOSS directly on the EMC energy over the MDC momentum of each track.")
  .note(:ks_mass_window,
        "Require |m(pi+ pi-) - m(K_S0)| < 35 MeV/c^2 after the secondary-vertex fit.")
  .note(:ks_decay_length_significance,
        "Require K_S^0 decay-length significance L/sigma_L > 2; no DSL primitive "
        "exists for this vertex-quality variable.")
  .note(:ks_momentum_match,
        "Require |p(pi+ pi-) - sqrt(s)/4| < 15 MeV/c, i.e. the reconstructed K_S "
        "momentum must match the expected value sqrt(s)/4 at each energy point.")
  .note(:conexc_mode_range,
        "ConExc mode 46 (K_S K_L) table is defined only up to 2.14 GeV; the "
        "higher-energy scan points (2.15-3.08 GeV) require validation of the "
        "ConExc model/table extrapolation before use.")
  .with_decay_card(decay_card_conexc)
  .apply(event_selection)

# Execute the algorithm on real data, inclusive MC and exclusive MC for every scan point
root_files = my_algorithm.execute_on(rscan_data + rscan_incMC + exMC_signal)