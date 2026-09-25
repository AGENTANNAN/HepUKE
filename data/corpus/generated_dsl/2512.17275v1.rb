# =============================================================================
# e+e- -> Xi(1530)^0 anti-Xi^0 across a cross-section scan (3.51 - 4.95 GeV)
# Single-baryon tagging: only the Xi(1530)^0 side is reconstructed, the
# anti-Xi^0 is extracted from the recoil-mass spectrum.
#   Xi(1530)^0 -> Xi- pi+
#   Xi-        -> Lambda pi-
#   Lambda     -> p pi-
#   (anti-Xi^0 -> anti-Lambda pi0, anti-Lambda -> anti-p pi+, pi0 -> gamma gamma)
# =============================================================================

### Dataset description ###
# Real-data scan points between 3.51 and 4.95 GeV (44.2 fb^-1 total)
scan_data = DatasetManager.real_data.where(cms_energy: { value: 3510.0..4950.0 })

# Decay card (EvtGen format) for the full signal chain
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 Xi(1530)0 anti-Xi0 PHSP;
    Enddecay

    Decay Xi(1530)0
    1.000 Xi- pi+ PHSP;
    Enddecay

    Decay Xi-
    1.000 Lambda0 pi- PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Xi0
    1.000 anti-Lambda0 pi0 PHSP;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 100k events at every scan point (no inclusive MC is used)
exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_Xi1530_Xibar0"   # auto-suffixed per energy point
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "Xi1530Xibar0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({ "ECMS" => [:double, 3.510] })   # beam energy (set per scan point)
            .set_alias({ "std::vector<double>" => "Vdouble" })

# Capture BOSS-side criteria that cannot be expressed with the DSL primitives
my_algorithm
  .note(:lambda_selection_cuts, "Lambda (p pi-): |M(p pi-) - m_Lambda| < 7 MeV; "
        "secondary-vertex fit chi2 < 500; decay length > 0")
  .note(:xi_minus_selection_cuts, "Xi- (Lambda pi-): |M(Lambda pi-) - m_Xi-| < 6.5 MeV; "
        "decay length > 0; best candidate taken closest to the nominal Xi- mass")

# Build the event-selection chain
event_selection = Selection.new
    .select_track {                    # Charged-track selection
        cos_theta 0.93                 # |cos(theta)| < 0.93
        nChrp     ">=2"                # at least 2 positive tracks
        nChrn     ">=2"                # at least 2 negative tracks
    }
    .pid(method: :probability) {       # Particle identification (probability method)
        prob_cut 0.001                 # PID probability > 0.001
        identify :proton, against: [:kaon, :pion]     # p+ / anti-p-
        identify :pion,   against: [:kaon, :proton]   # pi+ / pi-
        nprp ">=1"                     # at least 1 p
        nprm ">=1"                     # at least 1 anti-p
        npip ">=1"                     # at least 1 pi+
        npim ">=2"                     # at least 2 pi-
    }
    # Lambda -> p pi- (secondary-vertex fit)
    .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    # Xi- -> Lambda pi- (secondary-vertex fit)
    .secondary_vertex_fit([:Lambda, :pim]) {
        build_virtual_particle(:Xi_minus).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    # 4C kinematic fit of the tagged Xi(1530)^0 -> Xi- pi+;
    # the anti-Xi^0 is the recoiling missing particle
    .kinematic_fit([:Xi_minus, :pip]) {
        nominal
        constrain_four_momentum
        miss_track_of(:anti_Xi0)
        chi2_cut 200
    }

# Generate the algorithm for the process in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on the real-data scan points and the exclusive MC
root_files = my_algorithm.execute_on(scan_data + exMC_signal)