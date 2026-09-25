### Dataset description ###
# R = sigma(e+e- -> hadrons) / sigma(e+e- -> mu+mu-) measured at the 14 R-scan
# c.m. energies between 2.2324 and 3.6710 GeV.
# Real data AND inclusive MC are used at every scan point; NO exclusive MC is generated.
scan_sample_names = %w[713_2232 713_2309 713_2386 713_2396 713_2500
                       713_2644 713_2700 713_2800 713_2900 713_2950
                       713_2981 713_3000 713_3020 713_3080]

rscan_data  = scan_sample_names.map { |s| DatasetManager.real_data.find(s) }     # real data at each point
rscan_incMC = scan_sample_names.map { |s| DatasetManager.inclusive_mc.find(s) }  # inclusive MC at each point

# ConExc decay card (continuum / R-value measurement): models ISR (up to 2nd order)
# and the measured sigma_0(m). The DSL auto-detects the literal `ConExc` token, switches
# to the no-KKMC simulation template and injects `Particle vpho <ECMS>` per energy point,
# so `Particle vpho` is deliberately omitted for this multi-energy scan.
decay_card_rscan = <<~DECAYCARD
    Decay vpho
    1 ConExc 6;
    Enddecay
    End
DECAYCARD

### Event selection (BOSS) ###
alg_name = "RScanHadronic"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 2.7000]})   # varies per scan point; representative value
   .set_alias({"std::vector<double>" => "Vdouble"})

# BOSS-side steps with no dedicated DSL construct are recorded as notes.
alg.note(:background_veto,
         "reject e+e- -> e+e- and e+e- -> gamma gamma: require >= 2 EMC showers with " \
         "|Delta theta| < 10 deg and the second-most-energetic shower E > 0.65 * E_beam")
   .note(:gamma_conversion_veto,
         "remove photon conversions: opposite-charge track pair with E/(pc) > 0.8 for both " \
         "tracks, invariant mass < 0.1 GeV/c^2 and opening angle < 15 deg")
   .note(:deuteron_removal,
         "remove deuterons with chi_p < 10")
   .note(:track_quality_cuts,
         "charged tracks required p < 0.94 * p_beam; a track is removed if E/(pc) > 0.8 " \
         "and p > 0.65 * p_beam")
   .note(:two_prong_topology,
         "for 2-prong events the two tracks must not be back-to-back " \
         "(|Delta theta| < 10 deg and ||Delta phi| - 180 deg| < 15 deg); require N_iso_2prg > 1 " \
         "with isolated photon E > 100 MeV and angle to nearest track > 20 deg")
   .note(:three_prong_topology,
         "for 3-prong events the two highest-momentum tracks must not be back-to-back; " \
         "require fewer than 2 tracks with E/(pc) > 0.8 or r_PID > 0.25")
   .note(:multi_prong_topology,
         "events with more than 3 prongs are counted directly, without the 2/3-prong topology cuts")
   .note(:beam_background,
         "beam background estimated from the V_z sideband (5,10) cm versus the signal region (0,5) cm")

# Inclusive hadronic event selection.
event_selection = Selection.new
event_selection
    .select_track {                    # good charged tracks
        cos_theta   0.93               # |cos(theta)| < 0.93
        Vz          5.0                # |Vz| < 5 cm
        Vr          0.5                # Vr < 0.5 cm
        nTot        ">=2"              # at least two charged tracks (positive + negative) in total
    }
    .select_photon {                   # good photons
        tdc_emc_start     0            # EMC TDC window start (0-700)
        tdc_emc_end       700          # EMC TDC window end
        energyThreshold_b 0.1          # E > 100 MeV, barrel
        energyThreshold_e 0.1          # E > 100 MeV, endcap
        nGam             ">=1"         # at least one photon
    }
    .pid(method: :probability) {       # hadron PID (probability method)
        prob_cut 0.001                 # PID probability > 0.001
        identify :pion, against: [:kaon, :proton]   # pi+ and pi- against K and p
    }
    .assign({:chrgp => :pip, :chrgn => :pim})       # treat all remaining tracks as pions
    # Placeholder nominal 4C fit on pi+pi- (chi2 < 200) -- NOT used in the actual counting analysis.
    .kinematic_fit([:pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
    }

# Generate the algorithm for the ConExc continuum process and run over the scan.
alg.with_decay_card(decay_card_rscan).apply(event_selection)
root_files = alg.execute_on(rscan_data + rscan_incMC)