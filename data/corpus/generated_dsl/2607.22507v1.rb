# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# Reference/placeholder energy point: psi(3770) (712_3773). The identical selection
# chain is executed on all 56 fine scan points (3.510 - 4.951 GeV).
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the signal channel e+e- -> K_S0 anti-Xi+ Sigma- (charge conjugate implied)
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 K_S0 anti-Xi+ Sigma- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay anti-Xi+
    1.000 anti-Lambda0 pi+ PHSP;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ PHSP;
    Enddecay

    Decay Sigma-
    1.000 n0 pi- PHSP;
    Enddecay

    End
DECAYCARD

# 400k PHSP exclusive-MC events per energy point (reference point shown)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "KsXibarSigma_excl"
  config.related_dataset = data_3773
  config.events          = 400000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "KsXibarSigma"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.773]})   # reference point sqrt(s) = 3.773 GeV

# Build the (common) event selection chain: applied identically to real data,
# inclusive MC and signal MC on every scan point.
event_selection = Selection.new
  .select_track {                       # charged track selection
      cos_theta 0.93                    # |cos(theta)| < 0.93
      Vz        100.0                   # |Vz| < 100 cm
      Vr        10.0                    # Vr < 10 mm
      nChrp     ">=3"                   # at least 3 positively charged tracks
      nChrn     ">=2"                   # at least 2 negatively charged tracks
  }
  .pid(method: :probability) {          # PID by the probability method
      prob_cut 0.001                    # probability > 0.001
      identify :proton, against: [:kaon, :pion]   # p and anti-p vs K and pi
      identify :pion,   against: [:kaon, :proton] # pi+ and pi- vs K and p
      nprm ">=1"                        # at least 1 anti-proton
      npip ">=3"                        # at least 3 pi+
      npim ">=1"                        # at least 1 pi-
  }
  # Sequential secondary-vertex fits; each built by minimizing the mass difference
  # and removing the used tracks from their candidate lists.
  .secondary_vertex_fit([:pip, :pim]) {          # K_S0 -> pi+ pi-
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {          # anti-Lambda -> anti-p pi+
      build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:Lambda_bar, :pip]) {   # anti-Xi+ -> anti-Lambda pi+
      build_virtual_particle(:Xi_bar_plus).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  # Sigma- is NOT reconstructed: inferred by partial reconstruction from the recoil
  # against the reconstructed K_S0 anti-Xi+ system (recID 3 = Sigma- in the card).
  .partial_miss([3]) {
      require_recoil_mass 1.1, 1.3      # 1.1 < M_recoil < 1.3 GeV/c^2
  }

# Capture BOSS-side criteria the DSL cannot express as formal constructs.
my_algorithm
  .note(:KS0_selection,
        "K_S0 candidate required to satisfy |M(pi+pi-) - m_K_S0| < 12 MeV/c^2 and " \
        "decay-length significance L/dL > 2")
  .note(:antiLambda_selection,
        "anti-Lambda candidate required to satisfy |M(anti-p pi+) - m_Lambda| < 5 MeV/c^2 " \
        "and positive decay length")
  .note(:antiXi_selection,
        "anti-Xi+ candidate required to satisfy |M(anti-Lambda pi+) - m_Xi| < 8 MeV/c^2 " \
        "and positive decay length")
  .note(:best_combination,
        "best anti-Lambda / anti-Xi+ combination chosen by minimizing " \
        "delta_min = sqrt(|M(anti-p pi+) - m_Lambda|^2 + |M(anti-Lambda pi+) - m_Xi|^2)")
  .note(:energy_scan,
        "identical selection chain executed on each of the 56 c.m. energy scan points " \
        "(3.510 - 4.951 GeV); 712_3773 (psi(3770)) used as the reference/placeholder point")

# Generate the algorithm for the signal decay card and attach the selection.
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and signal MC (all scan points handled identically).
root_files = my_algorithm.execute_on([data_3773, incMC_3773, exMC_signal])