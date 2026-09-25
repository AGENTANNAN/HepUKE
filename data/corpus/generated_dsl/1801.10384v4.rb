# Core DSL classes and dependencies are loaded automatically at execution.
#
# Search for Z_s in e+e- -> phi pi pi at sqrt(s) = 2.125 GeV (R-scan point 713 Rscan_2125)
#   Mode I : e+e- -> phi pi+ pi-   (phi -> K+ K-)
#   Mode II: e+e- -> phi pi0 pi0   (phi -> K+ K-, pi0 -> gamma gamma)

### Dataset description ###
data_2125  = DatasetManager.real_data.find("713_Rscan_2125")     # 2.125 GeV real data
incMC_2125 = DatasetManager.inclusive_mc.find("713_Rscan_2125")  # matching inclusive MC

# ---- Decay card: Mode I  (e+e- -> phi pi+ pi-, phi -> K+ K-) ----
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.000 phi pi+ pi- PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# ---- Decay card: Mode II (e+e- -> phi pi0 pi0, phi -> K+ K-, pi0 -> gamma gamma) ----
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.000 phi pi0 pi0 PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---- Exclusive MC: 100k events for each channel ----
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_2125_phi_pip_pim"
    config.related_dataset = data_2125
    config.events          = 100000
    config.decay_card      = decay_card_modeI
    config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_2125_phi_pi0_pi0"
    config.related_dataset = data_2125
    config.events          = 100000
    config.decay_card      = decay_card_modeII
    config.cross_section   = :default
end

### Event selection (BOSS) ###

# ============ Mode I: e+e- -> phi pi+ pi- (phi -> K+ K-) ============
alg_name_modeI = "ZsPhiPiPi"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 2.125]})       # sqrt(s) = 2.125 GeV
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
sel_modeI.select_track {
        cos_theta 0.93    # |cos(theta)| < 0.93
        Vz        10.0    # |Vz| < 10 cm
        Vr        1.0     # Vr < 1 cm
        nTot      ">=3"   # at least three charged tracks
    }
    .pid(method: :probability) {
        prob_cut 0.001                             # PID probability > 0.001
        identify :kaon, against: [:pion, :proton]  # pi/K separation: K+/K-
        identify :pion, against: [:kaon, :proton]  # pi/K separation: pi+/pi-
        nkp  ">=1"   # at least one K
        npip ">=1"   # at least one pi+
        npim ">=1"   # at least one pi-
    }
    # Nominal 1C fit: K+ pi+ pi- with the K- allowed to be missing, 4-momentum constrained
    .kinematic_fit([:kp, :pip, :pim]) {
        nominal
        miss_track_of :km
        constrain_four_momentum
        chi2_cut 10
    }
    # Complementary 1C fit (charge-conjugate hypothesis): K- pi+ pi- with the K+ missing
    # (no chi2_cut / no nominal -> stores the competing chi2 for a ROOT-level veto)
    .kinematic_fit([:km, :pip, :pim]) {
        miss_track_of :kp
        constrain_four_momentum
    }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)


# ============ Mode II: e+e- -> phi pi0 pi0 (phi -> K+ K-, pi0 -> gamma gamma) ============
alg_name_modeII = "ZsPhiPi0Pi0"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 2.125]})      # sqrt(s) = 2.125 GeV
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
sel_modeII.select_track {
        cos_theta 0.93   # |cos(theta)| < 0.93
        Vz        10.0   # |Vz| < 10 cm
        Vr        1.0    # Vr < 1 cm
        nTot      ">=1"  # at least one charged track
    }
    .select_photon {
        tdc_emc_start     0      # EMC timing start (0 -> 700 ns)
        tdc_emc_end       14     # EMC timing end
        angle_to_track    10.0   # angle to nearest charged track > 10 deg
        energyThreshold_b 0.025  # E > 25 MeV in the barrel
        energyThreshold_e 0.050  # E > 50 MeV in the endcap
        nGam              ">=4"  # at least four photons
    }
    .pid(method: :probability) {
        prob_cut 0.001                             # PID probability > 0.001
        identify :kaon, against: [:pion, :proton]  # K+/K-
        nkp ">=1"   # at least one K
    }
    # 2C Kalman fit on the four photons: two gamma-gamma pairs constrained to the pi0 mass
    .kalman_kinematic_fit([:gamma, :gamma, :gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=2"   # at least two pi0 candidates
    }
    # Nominal 1C fit: K+ pi0 pi0 with the K- allowed to be missing, 4-momentum constrained
    .kinematic_fit([:kp, :pi0, :pi0]) {
        nominal
        miss_track_of :km
        constrain_four_momentum
        chi2_cut 20
    }
    # Complementary 1C fit (charge-conjugate hypothesis): K- pi0 pi0 with the K+ missing
    # (no chi2_cut / no nominal -> stores the competing chi2 for a ROOT-level veto)
    .kinematic_fit([:km, :pi0, :pi0]) {
        miss_track_of :kp
        constrain_four_momentum
    }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)


### Execution: both modes on the same data, inclusive MC, and their signal MC ###
root_files_modeI  = alg_modeI.execute_on([data_2125, incMC_2125, exMC_modeI])
root_files_modeII = alg_modeII.execute_on([data_2125, incMC_2125, exMC_modeII])