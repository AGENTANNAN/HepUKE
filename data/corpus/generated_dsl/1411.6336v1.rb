# =====================================================================
# BOSS DSL : e+e- -> gamma chi_cJ (J=0,1,2)
#            chi_cJ -> gamma J/psi ,  J/psi -> mu+ mu-
#            at sqrt(s) = 4.009, 4.230, 4.260 and 4.360 GeV
# =====================================================================

### Dataset preparation ###
# Real data at the four centre-of-mass energies (sample name = [BOSS version]_[Ecms MeV])
data_4009 = DatasetManager.real_data.find("703_4009")
data_4230 = DatasetManager.real_data.find("703_4230")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4360 = DatasetManager.real_data.find("703_4360")
data_points = [data_4009, data_4230, data_4260, data_4360]

# Corresponding inclusive MC samples at each energy
incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360")

# --- Decay cards (EvtGen format); one card per chi_cJ (J = 0, 1, 2) ---
decay_card_chi0 = <<~DECAYCARD
    Decay psi(4260)
    1.000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

decay_card_chi1 = <<~DECAYCARD
    Decay psi(4260)
    1.000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

decay_card_chi2 = <<~DECAYCARD
    Decay psi(4260)
    1.000 gamma chi_c2 PHSP;
    Enddecay

    Decay chi_c2
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# --- 200k-event exclusive MC for each chi_cJ at each CME point ---
exMC_chi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "exmc_chi_c0_to_gamma_jpsi_mumu"
    config.events        = 200_000
    config.decay_card    = decay_card_chi0
    config.cross_section = :default
end

exMC_chi1 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "exmc_chi_c1_to_gamma_jpsi_mumu"
    config.events        = 200_000
    config.decay_card    = decay_card_chi1
    config.cross_section = :default
end

exMC_chi2 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "exmc_chi_c2_to_gamma_jpsi_mumu"
    config.events        = 200_000
    config.decay_card    = decay_card_chi2
    config.cross_section = :default
end

### Event selection (BOSS) ###
# chi_c0, chi_c1 and chi_c2 share identical final states (gamma gamma mu+ mu-)
# and identical selection criteria -> a single Algorithm + Selection chain.
alg_name = "GammaChiCJToGammaJpsiMumu"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.26]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
    .select_track {                       # charged-track selection
        cos_theta 0.93                    # |cos(theta)| < 0.93
        Vz        10.0                    # |Vz| < 10 cm
        Vr        1.0                     # Vr < 1 cm
        nChrp     "==2"                   # exactly two positive tracks
        nChrn     "==2"                   # exactly two negative tracks
        nNet      "==0"                   # net charge zero
    }
    .select_photon {                      # photon selection
        tdc_emc_start     0               # EMC time window 0 ...
        tdc_emc_end       14              # ... to 700 ns
        angle_to_track    20.0            # photons at least 20 deg from any charged track
        energyThreshold_b 0.025           # E > 25 MeV in barrel   (|cos(theta)| < 0.80)
        energyThreshold_e 0.050           # E > 50 MeV in endcap   (0.86 < |cos(theta)| < 0.92)
        nGam              ">=2"           # at least two photon candidates
    }
    .pid(method: :probability) {           # probability-method PID
        prob_cut 0.001
        # high-momentum tracks (p > 1.0 GeV/c) treated as leptons;
        # electron if EMC energy > 0.6 GeV, otherwise muon
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.6
        identify :pion, against: [:kaon]   # pi/K separation
        nlp "==1"                          # exactly one positive lepton
        nlm "==1"                          # exactly one negative lepton
    }
    # 5C kinematic fit to gamma gamma mu+ mu- : 4C energy-momentum conservation
    # + 1C mass constraint M(mu+mu-) = M(J/psi).  The candidate with the smallest
    # chi2 is retained automatically when several survive.
    .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
        # additional background suppression: veto the gamma-gamma mass around pi0 / eta / eta'
        invariant_mass_of(:gamma, :gamma).out_of(0.110, 0.160)   # |M(gg) - M(pi0)| > 0.025
        invariant_mass_of(:gamma, :gamma).out_of(0.518, 0.578)   # |M(gg) - M(eta)| > 0.03
        invariant_mass_of(:gamma, :gamma).out_of(0.938, 0.978)   # |M(gg) - M(eta')| > 0.02
        chi2_cut 200                       # loose BOSS chi2 cut; tight published chi2 < 40 applied later in ROOT
    }

# Inexpressible BOSS-side procedures preserved for the systematic-uncertainty stage
my_algorithm
    .note(:muon_pid_criterion, "E/p < 0.35 is required for muon candidates in PID; this track-level E/p criterion is not expressible through identify_high_momentum_leptons and must be applied explicitly in the generated selection code")
    .note(:background_veto, "at sqrt(s) = 4.009 GeV the chi_c1,2 decay photon energy must be below 0.403 GeV to suppress background; this energy-dependent cut applies only to the 4.009 GeV sample and uses the measured (pre-ROOT) photon energy")
    .note(:sideband_study, "J/psi sidebands 2.917 < M(mu+mu-) < 3.057 and 3.137 < M(mu+mu-) < 3.277 GeV/c^2 are studied by performing additional 5C fits with M(mu+mu-) constrained to 3.047 or 3.147 GeV/c^2; the arbitrary sideband mass values cannot be expressed with constrain_to_nominal_mass_of")
    .note(:multi_energy_ecms, "four centre-of-mass energies are analysed; the ECMS constant is set to the nominal 4.260 GeV value and must be overridden per energy point for the 5C fit")

# Shared selection chain applied for chi_c0, chi_c1 and chi_c2
my_algorithm.with_decay_card(decay_card_chi0).apply(event_selection)

# Execute on real data + inclusive MC at each energy and all exclusive MC samples
root_files = my_algorithm.execute_on([data_4009, data_4230, data_4260, data_4360,
                                      incMC_4009, incMC_4230, incMC_4260, incMC_4360,
                                      *exMC_chi0, *exMC_chi1, *exMC_chi2])