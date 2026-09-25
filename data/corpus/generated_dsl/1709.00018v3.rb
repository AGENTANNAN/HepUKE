# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# Real data and inclusive MC at both energy points
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) (3.686 GeV) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # Corresponding inclusive MC

# Decay card: J/psi -> pi+ pi- eta', eta' -> eta pi+ pi-, eta -> gamma gamma
decay_card_jpsi = <<~DECAYCARD
    Decay J/psi
    1.0000 pi+ pi- eta'      PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi-       PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma       PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: psi(2S) -> pi+ pi- eta', eta' -> eta pi+ pi-, eta -> gamma gamma
decay_card_psip = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- eta'      PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi-       PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma       PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for J/psi -> pi+ pi- eta' (eta' -> eta pi+ pi-, eta -> gamma gamma)
exMC_jpsi = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_3097_pipi_etap_etapipigamgam"
    config.related_dataset = jpsi_data          # Associated real dataset
    config.events = 500_000                     # 500k events
    config.decay_card = decay_card_jpsi
    config.cross_section = :default
end

# Exclusive MC for psi(2S) -> pi+ pi- eta' (eta' -> eta pi+ pi-, eta -> gamma gamma)
exMC_psip = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_3686_pipi_etap_etapipigamgam"
    config.related_dataset = psip_data          # Associated real dataset
    config.events = 500_000                     # 500k events
    config.decay_card = decay_card_psip
    config.cross_section = :default
end

### Event selection (BOSS) ###
# J/psi (3.097 GeV) algorithm
alg_name_jpsi = "JpsiPiPiEtaP"
alg_jpsi = Algorithm.new(alg_name_jpsi)
alg_jpsi.set_header(["#{alg_name_jpsi}Alg/#{alg_name_jpsi}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})   # sqrt(s) = 3.097 GeV
        .set_alias({"std::vector<double>" => "Vdouble"})

# psi(2S) (3.686 GeV) algorithm
alg_name_psip = "PsipPiPiEtaP"
alg_psip = Algorithm.new(alg_name_psip)
alg_psip.set_header(["#{alg_name_psip}Alg/#{alg_name_psip}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})   # sqrt(s) = 3.686 GeV
        .set_alias({"std::vector<double>" => "Vdouble"})

# Common event selection chain shared by both energy points
common_selection = Selection.new
    .select_track {                     # Charged track selection
        cos_theta   0.93                # |cos(theta)| < 0.93
        Vz          10.0                # |Vz| < 10 cm
        Vr          1.0                 # Vr < 1 cm
        nChrp       "==2"               # Exactly 2 positive tracks
        nChrn       "==2"               # Exactly 2 negative tracks
        nNet        "==0"               # Net charge zero
    }
    .select_photon {                    # Photon selection
        tdc_emc_start     0             # EMC TDC start
        tdc_emc_end       14            # EMC TDC end
        angle_to_track    10.0          # At least 10 degrees from any charged track
        energyThreshold_b 0.025         # 25 MeV in the barrel region
        energyThreshold_e 0.050         # 50 MeV in the endcap region
        nGam              ">=2"         # At least two photons (from eta -> gamma gamma)
    }
    .pid(method: :probability) {        # Particle identification (probability method)
        prob_cut 0.001                  # PID probability > 0.001
        identify :pion, against: [:kaon]  # pi+ and pi- (charge-conjugation shorthand), separated from kaons
        npip ">=2"                      # At least 2 pi+
        npim ">=2"                      # At least 2 pi-
    }
    .kinematic_fit([:gamma, :gamma, :pip, :pim, :pip, :pim]) {  # 4C fit to gamma gamma pi+ pi- pi+ pi-
        nominal                                     # Nominal fit: corrected four-momenta are used downstream
        constrain_four_momentum                     # 4C energy-momentum constraint
        chi2_cut 40                                 # chi^2 < 40
    }

# Generate the complete algorithms for the two energy points (each requires its own decay card),
# sharing the identical selection chain.
alg_jpsi.with_decay_card(decay_card_jpsi).apply(common_selection)
alg_psip.with_decay_card(decay_card_psip).apply(common_selection)

# Execute the algorithms on the corresponding datasets, producing ROOT files
root_files_jpsi = alg_jpsi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi])
root_files_psip = alg_psip.execute_on([psip_data, psip_incMC, exMC_psip])