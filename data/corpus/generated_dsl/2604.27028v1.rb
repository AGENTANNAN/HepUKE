### Dataset description ###
# 3.097 GeV J/psi real data and the corresponding inclusive MC sample
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal process J/psi -> K- Sigma0 anti-Xi+ (charge conjugate included)
#   Sigma0 -> Lambda0 gamma,  anti-Xi+ -> anti-Lambda0 pi+,
#   anti-Lambda0 -> anti-p- pi+,  Lambda0 -> p+ pi-
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 K- Sigma0 anti-Xi+ PHSP;
    Enddecay

    Decay Sigma0
    1.0000 Lambda0 gamma PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 anti-Lambda0 pi+ PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive MC sample for the single decay mode (200k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3097_K_Sigma0_Xibar"
    config.related_dataset = jpsi_data
    config.events          = 200000
    config.decay_card      = decay_card_signal
    config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "KSigma0Xi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            # Procedures that cannot be expressed with a dedicated DSL primitive:
            .note(:recoil_mass_window,
                  "BOSS-level pre-select: recoil (missing) mass of the K- anti-Xi+ system " \
                  "required in [1.174, 1.204] GeV (m_Sigma0 +/- 15 MeV) before the 1C kinematic " \
                  "fit; no recoil_mass_of primitive available in the DSL")
            .note(:charge_conjugate,
                  "charge-conjugate channel J/psi -> K+ anti-Sigma0 Xi- (anti-Sigma0 -> " \
                  "anti-Lambda0 gamma, Xi- -> Lambda0 pi-, Lambda0 -> p+ pi-, anti-Lambda0 -> " \
                  "anti-p- pi+) analysed with the mirrored selection")

# Build the event-selection chain
event_selection = Selection.new
    .select_track {                     # loose charged-track selection (long-lived Lambda/Xi)
        cos_theta 0.93                  # |cos(theta)| < 0.93
        Vz        20.0                  # |Vz| < 20 cm
        Vr        10.0                  # Vxy < 10 cm
        nChrp     ">=2"                 # at least two positive tracks
        nChrn     ">=1"                 # at least one negative track
    }
    .pid(method: :probability) {        # probability-method PID, prob > 0.001
        prob_cut 0.001
        identify :proton, against: [:kaon,   :pion]   # p / anti-p vs K and pi
        identify :kaon,   against: [:proton, :pion]   # K+ / K- vs p and pi
        identify :pion,   against: [:proton, :kaon]   # pi+ / pi- vs p and K
        nprm ">=1"                      # at least one anti-proton
        nkm  ">=1"                      # at least one K-
        npip ">=2"                      # at least two pi+
    }
    .remove(:km) { condition "Vxy > 1.0 || Vz > 5.0" }   # K- must originate from the IP
    .secondary_vertex_fit([:prm, :pip]) {                # anti-Lambda -> anti-p pi+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:Lambda_bar, :pip]) {         # anti-Xi+ -> anti-Lambda pi+
        build_virtual_particle(:Xi_bar_plus).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:km, :Xi_bar_plus, :Sigma0]) {       # 1C fit: missing Sigma0 constrained to nominal mass
        nominal                                          # nominal fit (four-momenta taken from this fit)
        miss_track_of(:Sigma0)                           # Sigma0 treated as missing
        constrain_four_momentum
        # no chi2_cut: the fit is used only to improve the missing-mass resolution
    }

# Attach the decay card, render the selection, execute on data + MC
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])