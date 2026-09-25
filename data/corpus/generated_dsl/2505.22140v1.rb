# ============================================================================
# Search for the dark baryon chi in J/psi -> Xi- Xibar+,
#   Xi-  -> pi- chi (chi invisible),
#   Xibar+ -> pi+ Lambdabar,  Lambdabar -> pbar pi+
# J/psi (3.097 GeV), 1.0087e10 real events.
# ============================================================================

### ------------------------- Dataset preparation ------------------------- ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi real data (3.097 GeV, 1.0087e10 evts)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

# ---- Decay cards: full J/psi -> Xi- Xibar+ chain with the invisible chi. ----
# One card per m_chi hypothesis.  The dark baryon chi is declared as a stable
# (undecayed) EvtGen particle so that it leaves the detector invisible.
decay_card_mchi_107 = <<~DECAYCARD
    Particle chi 1.07 0.0 0.0 0.0 0 0 0 0 0 0 0

    Decay J/psi
    1.0000 Xi- anti-Xi+       PHSP;
    Enddecay

    Decay Xi-
    1.0000 pi- chi            PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 pi+ anti-Lambda0   PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+        HypWK;
    Enddecay

    End
DECAYCARD

decay_card_mchi_110 = <<~DECAYCARD
    Particle chi 1.10 0.0 0.0 0.0 0 0 0 0 0 0 0

    Decay J/psi
    1.0000 Xi- anti-Xi+       PHSP;
    Enddecay

    Decay Xi-
    1.0000 pi- chi            PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 pi+ anti-Lambda0   PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+        HypWK;
    Enddecay

    End
DECAYCARD

decay_card_mchi_1157 = <<~DECAYCARD
    Particle chi 1.1157 0.0 0.0 0.0 0 0 0 0 0 0 0

    Decay J/psi
    1.0000 Xi- anti-Xi+       PHSP;
    Enddecay

    Decay Xi-
    1.0000 pi- chi            PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 pi+ anti-Lambda0   PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+        HypWK;
    Enddecay

    End
DECAYCARD

decay_card_mchi_113 = <<~DECAYCARD
    Particle chi 1.13 0.0 0.0 0.0 0 0 0 0 0 0 0

    Decay J/psi
    1.0000 Xi- anti-Xi+       PHSP;
    Enddecay

    Decay Xi-
    1.0000 pi- chi            PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 pi+ anti-Lambda0   PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+        HypWK;
    Enddecay

    End
DECAYCARD

decay_card_mchi_116 = <<~DECAYCARD
    Particle chi 1.16 0.0 0.0 0.0 0 0 0 0 0 0 0

    Decay J/psi
    1.0000 Xi- anti-Xi+       PHSP;
    Enddecay

    Decay Xi-
    1.0000 pi- chi            PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 pi+ anti-Lambda0   PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+        HypWK;
    Enddecay

    End
DECAYCARD

# ---- Exclusive MC: 500k events of the full chain, one sample per m_chi ----
exMC_mchi_107 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_jpsi_pichi_mchi_1p07"
    config.related_dataset = jpsi_data
    config.events          = 500_000
    config.decay_card      = decay_card_mchi_107
    config.cross_section   = :default
end

exMC_mchi_110 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_jpsi_pichi_mchi_1p10"
    config.related_dataset = jpsi_data
    config.events          = 500_000
    config.decay_card      = decay_card_mchi_110
    config.cross_section   = :default
end

exMC_mchi_1157 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_jpsi_pichi_mchi_mLambda"
    config.related_dataset = jpsi_data
    config.events          = 500_000
    config.decay_card      = decay_card_mchi_1157
    config.cross_section   = :default
end

exMC_mchi_113 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_jpsi_pichi_mchi_1p13"
    config.related_dataset = jpsi_data
    config.events          = 500_000
    config.decay_card      = decay_card_mchi_113
    config.cross_section   = :default
end

exMC_mchi_116 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_jpsi_pichi_mchi_1p16"
    config.related_dataset = jpsi_data
    config.events          = 500_000
    config.decay_card      = decay_card_mchi_116
    config.cross_section   = :default
end

### ------------------------- Event selection (BOSS) ------------------------- ###
alg_name  = "XiDarkBaryon"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})   # J/psi centre-of-mass energy (GeV)

event_selection = Selection.new

event_selection
    .select_track {                       # charged track quality cuts
        cos_theta  0.93                   # |cos(theta)| < 0.93
        Vz         10.0                   # |Vz| < 10 cm
        Vr         1.0                    # Vr   < 1 cm
        nChrp      ">=2"                  # at least two positive tracks
        nChrn      ">=3"                  # at least three negative tracks
        nNet       ">=0"                  # net charge >= 0
    }
    .select_photon {                      # good photon (shower) definition
        tdc_emc_start     0               # EMC timing window 0 - 700 ns
        tdc_emc_end       14              # (unit = 50 ns, i.e. 14 -> 700 ns)
        angle_to_track    10.0            # > 10 deg from nearest charged track
        energyThreshold_b 0.025           # 25 MeV in the barrel (|cos(theta)| < 0.80)
        energyThreshold_e 0.050           # 50 MeV in the endcap (0.86 < |cos(theta)| < 0.92)
    }
    .pid(method: :probability) {          # probability-method PID
        prob_cut 0.001
        identify :prm, against: [:pim, :km]   # pbar separated from pi- and K-
        identify :pim, against: [:prm, :km]   # pi- separated from pbar and K-
        npim "==1"                            # double tag: exactly one additional pi-
    }
    .select_isolated_photon {             # reject showers from (anti)proton interactions
        angle_to_prm_track 20.0           # > 20 deg from the antiproton track
    }
    .assign({:chrgp => :pip})             # remaining positive tracks are pi+
    # --- Single-tag Xibar+ -> Lambdabar pi+, Lambdabar -> pbar pi+ ---
    .secondary_vertex_fit([:prm, :pip]) {                 # Lambdabar -> pbar pi+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:Lambda_bar, :pip]) {          # Xibar+ -> Lambdabar pi+
        build_virtual_particle(:Xi_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    # --- Double tag: nominal kinematic fit (4C + mass constraints) ---
    .kinematic_fit([:Xi_bar, :pim]) {
        nominal                                   # nominal fit -> corrected four-momenta are saved
        constrain_four_momentum                   # four-momentum conservation
        chi2_cut 20
    }
    # --- Competing Xi- -> pi- Lambda hypothesis (non-nominal, no chi2 cut:
    #     only the chi2 value is stored for a ROOT-level veto) ---
    .kinematic_fit([:Xi_bar, :pim]) {
        miss_track_of(:Lambda)                    # Lambda left undetected (competing hypothesis)
        constrain_four_momentum
    }

# --- BOSS-side procedures not expressible with the current DSL ---
algorithm
    .note(:vertex_fit_chi2, "Both secondary vertex fits (Lambdabar -> pbar pi+ and Xibar+ -> Lambdabar pi+) are required to have a vertex-fit chi2 < 200; the secondary_vertex_fit DSL block has no chi2_cut option, so this cut is imposed in the generated BOSS code and validated on signal MC.")
    .note(:tag_best_candidate, "Single-tag best-candidate selection: choose the candidate minimising |M(pbar pi+) - M_Lambda| + |M(Lambdabar pi+) - M_Xibar|, and require the Lambdabar candidate to satisfy |M(pbar pi+) - M_Lambda| < 0.004 GeV/c^2. The combined metric / mass window is applied at the ROOT level.")
    .note(:xi_recoil_mass, "Xibar+ recoil mass (equivalent to M(Xi-)) required to lie in (1.290, 1.345) GeV/c^2; not expressible as a dedicated BOSS-side cut, applied at the ROOT level.")
    .note(:missing_chi_constraint, "1C kinematic constraint of the missing mass (the invisible chi) to m_chi, with chi2_pi-chi < 20 - applied per mass hypothesis; the present kinematic_fit DSL cannot constrain an invisible particle's mass.")
    .note(:background_veto, "Competing Xi- -> pi- Lambda hypothesis rejected by chi2_pi-Lambda > chi2_pi-chi (chi2 of the second, non-nominal kinematic fit), except for the m_chi = m_Lambda hypothesis where the two interpretations coincide; the comparison is done at the ROOT level.")
    .note(:pi_momentum_window, "The pi- momentum window (0.070-0.192 GeV/c, different for each m_chi hypothesis) is optimised with the Punzi figure of merit and applied in the ROOT analysis.")

# Attach a representative decay card (the generated header defines the same
# visible particles for all five m_chi hypotheses).
algorithm
    .with_decay_card(decay_card_mchi_107)
    .apply(event_selection)

# --- Execute on real data, inclusive MC and all five signal MC hypotheses ---
root_files = algorithm.execute_on(
    [jpsi_data, jpsi_incMC,
     exMC_mchi_107, exMC_mchi_110, exMC_mchi_1157, exMC_mchi_113, exMC_mchi_116]
)