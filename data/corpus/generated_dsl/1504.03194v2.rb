### Dataset description ###
# J/psi resonance at 3.097 GeV -> BOSS 7.0.8 sample "708_3097"
jpsi_data = DatasetManager.real_data.find("708_3097")        # 1.311e9-event J/psi real data sample
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Matching inclusive J/psi MC

### Decay cards (EvtGen format) ###
# Signal: J/psi -> phi pi0 -> K+ K- gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 phi pi0 PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Coherent (irreducible) background: J/psi -> K+ K- pi0
decay_card_bkg_KKpi0 = <<~DECAYCARD
    Decay J/psi
    1.0000 K+ K- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Background: J/psi -> phi pi0 pi0
decay_card_bkg_phipi0pi0 = <<~DECAYCARD
    Decay J/psi
    1.0000 phi pi0 pi0 PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Background: J/psi -> phi gamma gamma
decay_card_bkg_phigammagamma = <<~DECAYCARD
    Decay J/psi
    1.0000 phi gamma gamma PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Background: J/psi -> gamma eta_c(1S), eta_c(1S) -> K+ K- pi0
decay_card_bkg_eta_c = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta_c(1S) PHSP;
    Enddecay

    Decay eta_c(1S)
    1.0000 K+ K- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Background: J/psi -> gamma K+ K-
decay_card_bkg_gammaKK = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma K+ K- PHSP;
    Enddecay

    End
DECAYCARD

# Background: J/psi -> gamma pi0 K+ K-
decay_card_bkg_gammapi0KK = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma pi0 K+ K- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Background: J/psi -> pi0 pi0 K+ K-
decay_card_bkg_pi0pi0KK = <<~DECAYCARD
    Decay J/psi
    1.0000 pi0 pi0 K+ K- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Background: J/psi -> phi eta, eta -> gamma gamma
decay_card_bkg_phieta = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples ###
# Signal
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_phiPi0_KKgg"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Coherent background J/psi -> K+ K- pi0
exMC_KKpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_KKpi0"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_bkg_KKpi0
  config.cross_section   = :default
end

# J/psi -> phi pi0 pi0
exMC_phipi0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_phiPi0Pi0"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_phipi0pi0
  config.cross_section   = :default
end

# J/psi -> phi gamma gamma
exMC_phigammagamma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_phiGammaGamma"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_phigammagamma
  config.cross_section   = :default
end

# J/psi -> gamma eta_c(1S), eta_c(1S) -> K+ K- pi0
exMC_eta_c = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaEtac_KKpi0"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_eta_c
  config.cross_section   = :default
end

# J/psi -> gamma K+ K-
exMC_gammaKK = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaKK"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_gammaKK
  config.cross_section   = :default
end

# J/psi -> gamma pi0 K+ K-
exMC_gammapi0KK = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammapi0KK"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_gammapi0KK
  config.cross_section   = :default
end

# J/psi -> pi0 pi0 K+ K-
exMC_pi0pi0KK = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_pi0pi0KK"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_pi0pi0KK
  config.cross_section   = :default
end

# J/psi -> phi eta, eta -> gamma gamma
exMC_phieta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_phiEta_gg"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_phieta
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiToPhiPi0"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
    .select_track {                       # Charged track selection
        cos_theta 0.93                    # |cos(theta)| < 0.93
        Vz        10.0                    # |Vz| < 10 cm
        Vr        1.0                     # Vr < 1 cm
        nChrp     "==1"                   # exactly one positive track
        nChrn     "==1"                   # exactly one negative track
        nNet      "==0"                   # net charge zero
    }
    .select_photon {                      # Photon selection
        tdc_emc_start     0               # EMC timing start 0
        tdc_emc_end       14              # EMC timing end 14
        angle_to_track    20.0            # angle to nearest charged track > 20 degrees
        energyThreshold_b 0.025           # > 25 MeV in the barrel
        energyThreshold_e 0.050           # > 50 MeV in the endcap
        nGam              ">=2"           # at least two photons
    }
    .pid(method: :probability) {          # Particle identification (probability method)
        prob_cut 0.001                    # PID probability > 0.001
        identify :kaon, against: [:pion, :proton]  # K+ and K- (charge-conjugation shorthand)
        nkp "==1"                         # exactly one K+
        nkm "==1"                         # exactly one K-
    }
    # 4C kinematic fit to the K+K-gamma gamma final state; best (smallest chi2) combination chosen automatically
    .kinematic_fit([:kp, :km, :gamma, :gamma]) {
        nominal                           # nominal fit: corrected four-momenta are saved
        constrain_four_momentum           # 4C energy-momentum constraint
        chi2_cut 200                      # loose chi2 < 200 (tight chi2 < 30 applied in ROOT after the fit)
    }

# Generate the algorithm for the signal process
my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC, and all exclusive MC samples
root_files = my_Algorithm.execute_on([
  jpsi_data, jpsi_incMC,
  exMC_signal,
  exMC_KKpi0,
  exMC_phipi0pi0,
  exMC_phigammagamma,
  exMC_eta_c,
  exMC_gammaKK,
  exMC_gammapi0KK,
  exMC_pi0pi0KK,
  exMC_phieta
])