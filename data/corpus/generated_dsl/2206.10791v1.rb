# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
data = DatasetManager.real_data.find("708_3097")     # J/psi (3.097 GeV) real data
incMC = DatasetManager.inclusive_mc.find("708_3097") # Corresponding inclusive MC

# Decay card for J/psi -> Lambda anti-Lambda (tag: anti-Lambda -> anti-p- pi+, signal: Lambda -> n gamma)
decay_card_for_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay Lambda0
    1.0000 n0 gamma PHSP;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive MC for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3097_LLbar_ngamma"
  config.related_dataset = data
  config.events = 200000
  config.decay_card = decay_card_for_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiLLbarNgamma"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
        .set_constant({"ECMS" => [:double, 3.097]}) # double ECMS = 3.097 GeV
        .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
    .select_track {                     # Charged track selection
        cos_theta        0.93           # |cos(theta)| < 0.93
        nTot             ">=2"          # at least two charged tracks
    }
    .select_photon {                    # Photon selection
        tdc_emc_start    0              # EMC timing 0-700 ns
        tdc_emc_end      14
        energyThreshold_b 0.025         # 25 MeV in the barrel region
        energyThreshold_e 0.050         # 50 MeV in the endcap region
        angle_to_track   10.0           # > 10 deg from nearest charged track
        nGam             ">=1"          # at least one photon
    }
    .pid(method: :probability) {        # PID by the probability method
        prob_cut         0.001          # probability > 0.001
        identify :prm, against: [:kaon, :pion]  # anti-p vs K and pi
        nprm             ">=1"          # at least one anti-proton
    }
    .remove([:prm <= :chrgn])           # remove identified anti-protons from negative charged list
    .select_isolated_photon {           # Isolated photon selection
        angle_to_prm_track 20.0         # > 20 deg from the anti-proton track
        nGam             ">=1"          # at least one isolated photon
    }
    .assign({:chrgp => :pip, :chrgn => :pim}) # treat remaining tracks as pions
    .secondary_vertex_fit([:prm, :pip]) {     # Reconstruct tag anti-Lambda from anti-p pi+
        build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    # Nominal 1C kinematic fit: J/psi -> Lambda_bar + (Lambda -> n gamma) with the neutron missing
    .kinematic_fit([:Lambda_bar, :gamma]) {
        nominal                          # mark as the nominal fit
        miss_track_of(:neutron)          # neutron is missing
        constrain_four_momentum          # constrain total four-momentum to CMS energy
        chi2_cut 10                      # chi2 < 10
    }
    # Competing 3C hypothesis J/psi -> Lambda_bar n gamma gamma used to veto pi0
    .kinematic_fit([:Lambda_bar, :gamma, :gamma]) {
        miss_track_of(:neutron)
        constrain_four_momentum
        invariant_mass_of(:gamma, :gamma).out_of(0.115, 0.155) # reject M(gamma gamma) within 20 MeV of the pi0 mass
    }

# Capture BOSS-side procedures that cannot be expressed in the DSL
algorithm
    .note(:secondary_vertex_cuts, "anti-Lambda -> anti-p pi+ secondary-vertex fit requires chi2 < 20; M(anti-p pi+) within 8 MeV/c^2 of the Lambda mass; and decay-length significance > 2 sigma")
    .note(:neutral_recoil_cuts, "recoil mass against the anti-Lambda tag required in 1.03-1.18 GeV/c^2")
    .note(:photon_bdt, "signal photon required to have a photon BDT score > 0.3")
    .note(:photon_energy, "signal photon EMC energy required > 150 MeV and < 400 MeV")
    .note(:neutron_opening_angle, "opening angle between the signal photon and the missing-neutron candidate required > 20 deg")
    .with_decay_card(decay_card_for_signal)
    .apply(event_selection)

# Execute on real data, inclusive MC, and signal exclusive MC
root_files = algorithm.execute_on([data, incMC, exMC_signal])