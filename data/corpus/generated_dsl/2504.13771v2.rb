### Dataset description ###
# Real data and corresponding inclusive MC for the two resonances
jpsi_data  = DatasetManager.real_data.find("708_3097")        # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")     # J/psi inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")        # psi(3686) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")     # psi(3686) inclusive MC

# Decay card: J/psi -> K_S0 K_S0, K_S0 -> pi+ pi-
decay_card_jpsi = <<~DECAYCARD
    Decay J/psi
    1.000 K_S0 K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: psi(3686) -> K_S0 K_S0, K_S0 -> pi+ pi-
decay_card_psip = <<~DECAYCARD
    Decay psi(2S)
    1.000 K_S0 K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples, 100k events each
exMC_jpsi = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_jpsi_ksks"
    config.related_dataset = jpsi_data
    config.events          = 100000
    config.decay_card      = decay_card_jpsi
    config.cross_section   = :default
end

exMC_psip = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_ksks"
    config.related_dataset = psip_data
    config.events          = 100000
    config.decay_card      = decay_card_psip
    config.cross_section   = :default
end

### Event selection (BOSS) ###
# Two independent algorithms: same final state but different K_S0 mass window and nominal chi2 cut
alg_name_jpsi = "JpsiKSKS"
alg_jpsi = Algorithm.new(alg_name_jpsi)
alg_jpsi.set_header(["#{alg_name_jpsi}Alg/#{alg_name_jpsi}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})
        # BOSS-side procedures that cannot be expressed in the DSL structure, kept for later systematic studies
        .note(:ks0_vertex_quality, "each K_S0 candidate is required to satisfy the secondary vertex fit chi2 < 200; this quality requirement is applied on the K_S0 built by build_virtual_particle and is enforced outside the DSL secondary_vertex_fit block")
        .note(:ks0_mass_window, "each reconstructed K_S0 is required to satisfy |M(pi+pi-) - m_K_S0| < 18 MeV/c^2; the window is imposed on the secondary-vertex K_S0 mass")
        .note(:ks0_decay_length, "each reconstructed K_S0 is required to have a decay length significance > 2 sigma; applied on the K_S0 flight length from the secondary vertex fit")

alg_name_psip = "PsipKSKS"
alg_psip = Algorithm.new(alg_name_psip)
alg_psip.set_header(["#{alg_name_psip}Alg/#{alg_name_psip}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .note(:ks0_vertex_quality, "each K_S0 candidate is required to satisfy the secondary vertex fit chi2 < 200; this quality requirement is applied on the K_S0 built by build_virtual_particle and is enforced outside the DSL secondary_vertex_fit block")
        .note(:ks0_mass_window, "each reconstructed K_S0 is required to satisfy |M(pi+pi-) - m_K_S0| < 30 MeV/c^2; the window is imposed on the secondary-vertex K_S0 mass")
        .note(:ks0_decay_length, "each reconstructed K_S0 is required to have a decay length significance > 2 sigma; applied on the K_S0 flight length from the secondary vertex fit")

# Common selection chain shared by both modes
event_selection_common = Selection.new
event_selection_common.select_track {      # four good charged tracks, net charge zero
        cos_theta 0.93                     # |cos(theta)| < 0.93
        Vz        10.0                     # |Vz| < 10 cm
        Vr        1.0                      # Vr < 1 cm
        nChrp     "==2"                    # exactly 2 positive tracks
        nChrn     "==2"                    # exactly 2 negative tracks
        nNet      "==0"                    # net charge zero
    }
    .select_photon {                       # at least one good photon (veto source)
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025            # E > 25 MeV in the barrel
        energyThreshold_e 0.050            # E > 50 MeV in the endcap
        angle_to_track    10.0             # angle to nearest charged track > 10 degrees
        nGam              ">=1"            # at least one photon
    }
    .pid(method: :probability) {           # pion identification
        prob_cut 0.001                     # PID probability > 0.001
        identify :pion, against: [:kaon, :proton]  # pi+ and pi- against K and p
        npip     "==2"                     # exactly 2 pi+ candidates
        npim     "==2"                     # exactly 2 pi- candidates
    }

# --- J/psi -> K_S0 K_S0 chain ---
jpsi_selection = event_selection_common.dup
jpsi_selection
    .secondary_vertex_fit([:pip, :pim]) {                 # first K_S0 from a pi+ pi- pair
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:pip, :pim]) {                 # second K_S0 from the remaining pi+ pi- pair
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:K_S0, :K_S0]) {                      # nominal 4C fit to K_S0 K_S0
        nominal
        constrain_four_momentum
        chi2_cut 15                                       # J/psi mode: chi2 < 15
    }
    .kinematic_fit([:gamma, :K_S0, :K_S0]) {              # competing hypothesis: gamma K_S0 K_S0
        constrain_four_momentum                           # chi2 stored for later ROOT-level veto
    }

# --- psi(3686) -> K_S0 K_S0 chain ---
psip_selection = event_selection_common.dup
psip_selection
    .secondary_vertex_fit([:pip, :pim]) {                 # first K_S0 from a pi+ pi- pair
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:pip, :pim]) {                 # second K_S0 from the remaining pi+ pi- pair
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:K_S0, :K_S0]) {                      # nominal 4C fit to K_S0 K_S0
        nominal
        constrain_four_momentum
        chi2_cut 30                                       # psi(3686) mode: chi2 < 30
    }
    .kinematic_fit([:gamma, :K_S0, :K_S0]) {              # competing hypothesis: gamma K_S0 K_S0
        constrain_four_momentum
    }

# Generate the algorithms for the two decays
alg_jpsi.with_decay_card(decay_card_jpsi).apply(jpsi_selection)
alg_psip.with_decay_card(decay_card_psip).apply(psip_selection)

# Execute on real data, inclusive MC and signal exclusive MC
root_files_jpsi = alg_jpsi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi])
root_files_psip = alg_psip.execute_on([psip_data, psip_incMC, exMC_psip])