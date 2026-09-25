# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Four real-data samples; the 4.226 GeV point is the nominal one
data_4226 = DatasetManager.real_data.find("703_4230")   # sqrt(s) = 4.226 GeV (nominal point)
data_4360 = DatasetManager.real_data.find("703_4360")   # sqrt(s) = 4.360 GeV
data_4420 = DatasetManager.real_data.find("703_4420")   # sqrt(s) = 4.420 GeV
data_4600 = DatasetManager.real_data.find("703_4600")   # sqrt(s) = 4.600 GeV

# Signal decay card (EvtGen format): e+e- -> pi+ pi- pi0 eta_c, eta_c -> p pbar, pi0 -> gamma gamma
# (the representative p pbar channel)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- pi0 eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0000 p+ anti-p- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Peaking-background decay card (EvtGen format): e+e- -> pi+ pi- h_c
decay_card_bkg = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 p+ anti-p- pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event signal exclusive MC, generated at the nominal 4.226 GeV point
exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_4226_pipipi0_etac_ppbar"
    config.related_dataset = data_4226
    config.events = 500000
    config.decay_card = decay_card_signal
    config.cross_section = :default
end

# 600k-event peaking-background exclusive MC e+e- -> pi+ pi- h_c
exMC_bkg = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_4226_pipi_hc"
    config.related_dataset = data_4226
    config.events = 600000
    config.decay_card = decay_card_bkg
    config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PipPimPi0EtaC"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
    .set_constant({"ECMS" => [:double, 4.226]})          # nominal energy point
    .set_alias({"std::vector<double>" => "Vdouble"})

# Build the event selection chain
event_selection = Selection.new
event_selection
    .select_track {                       # charged track selection
        cos_theta 0.93                    # |cos(theta)| < 0.93
        Vz 10.0                           # |Vz| < 10 cm
        Vr 1.0                            # Vr < 1 cm
        nChrp ">=1"                       # at least one positive track
        nChrn ">=1"                       # at least one negative track
        nTot ">=4"                        # at least four charged tracks in total
        nNet "==0"                        # net charge zero
    }
    .select_photon {                      # photon selection
        tdc_emc_start 0                   # EMC timing window 0-14
        tdc_emc_end 14
        energyThreshold_b 0.025           # barrel  E > 25 MeV
        energyThreshold_e 0.050           # endcap  E > 50 MeV
        angle_to_track 10.0               # at least 10 degrees from any charged track
        nGam ">=2"                        # at least two photons
    }
    .pid(method: :probability) {          # probability-based PID
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]   # p and pbar (charge-conjugation shorthand)
        identify :pion, against: [:kaon, :proton]   # pi+ and pi-
        nprp ">=1"
        nprm ">=1"
        npip ">=1"
        npim ">=1"
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {         # 1C fit: reconstruct pi0 -> gamma gamma
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
    }
    .kinematic_fit([:pip, :pim, :prp, :prm, :pi0]) {  # nominal 4C fit to pi+ pi- p pbar pi0
        nominal
        constrain_four_momentum
        chi2_cut 200
    }

# BOSS-side selection procedures that have no formal DSL construct
my_algorithm
    .note(:background_veto, "D-meson, K*(892), omega and eta vetoes applied to suppress
        peaking backgrounds; veto windows are optimised on signal and inclusive MC in ROOT")

# Attach the representative (p pbar) decay card and render the selection
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on the four real-data points and both exclusive MC samples
root_files = my_algorithm.execute_on([data_4226, data_4360, data_4420, data_4600, exMC_signal, exMC_bkg])