### Dataset description ###
data_3097  = DatasetManager.real_data.find("708_3097")      # J/psi real data at 3.097 GeV
incMC_3097 = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC at 3.097 GeV

# ---------------- Decay cards (EvtGen format) ----------------
# Signal: J/psi -> gamma eta', eta' -> gamma gamma pi0, pi0 -> gamma gamma (five photons)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma gamma pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Background: J/psi -> gamma eta', eta' -> gamma omega, omega -> gamma pi0
decay_card_bkg_etap_gomega = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma omega PHSP;
    Enddecay

    Decay omega
    1.0000 gamma pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Background: J/psi -> gamma eta', eta' -> pi0 pi0 eta, eta -> gamma gamma
decay_card_bkg_etap_2pi0eta = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 pi0 pi0 eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Background: J/psi -> gamma eta', eta' -> pi0 pi0 pi0
decay_card_bkg_etap_3pi0 = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 pi0 pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Background: J/psi -> gamma pi0 pi0
decay_card_bkg_g2pi0 = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Background: J/psi -> omega eta, omega -> gamma pi0, eta -> gamma gamma
decay_card_bkg_omega_eta = <<~DECAYCARD
    Decay J/psi
    1.0000 omega eta PHSP;
    Enddecay

    Decay omega
    1.0000 gamma pi0 PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---------------- Exclusive MC samples ----------------
# Signal: 1,000,000 events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_jpsi_gamma_etap_ggpi0"
  config.related_dataset = data_3097
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Background samples: 500,000 events each
exMC_bkg_etap_gomega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_etap_gomega"
  config.related_dataset = data_3097
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_etap_gomega
  config.cross_section   = :default
end

exMC_bkg_etap_2pi0eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_etap_2pi0eta"
  config.related_dataset = data_3097
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_etap_2pi0eta
  config.cross_section   = :default
end

exMC_bkg_etap_3pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_etap_3pi0"
  config.related_dataset = data_3097
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_etap_3pi0
  config.cross_section   = :default
end

exMC_bkg_g2pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_jpsi_gamma_2pi0"
  config.related_dataset = data_3097
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_g2pi0
  config.cross_section   = :default
end

exMC_bkg_omega_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_jpsi_omega_eta"
  config.related_dataset = data_3097
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_omega_eta
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiGammaEtaP"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy 3.097 GeV

# Common event selection applied to data, inclusive MC, signal MC and all background MC
event_selection = Selection.new
  .select_track {          # Fully neutral final state: require zero charged tracks
    nChrp "==0"            # no positive tracks (no |cos(theta)|, |Vz|, Vr track cuts applied)
    nChrn "==0"            # no negative tracks
    nNet  "==0"            # net charge zero
  }
  .select_photon {         # Photon selection
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0   # > 10 deg from any charged track
    energyThreshold_b 0.025  # > 25 MeV in barrel (|cos(theta)| < 0.80)
    energyThreshold_e 0.050  # > 50 MeV in end-cap (0.86 < |cos(theta)| < 0.92)
    nGam              "==5"  # exactly five photon candidates
  }
  # No charged-particle PID: the final state is fully neutral
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma]) {
    nominal                   # nominal fit; only its corrected four-momenta are saved
    constrain_four_momentum   # 4C energy-momentum conservation
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # +1C: pi0 mass constraint on one photon pair -> 5C
    chi2_cut 200              # loose chi2 cut in BOSS (paper chi2 < 30 is applied later outside BOSS)
  }

# Attach the signal decay card and render the common selection into the BOSS algorithm
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC, signal MC and all background MC samples
root_files = my_algorithm.execute_on([
  data_3097, incMC_3097,
  exMC_signal,
  exMC_bkg_etap_gomega, exMC_bkg_etap_2pi0eta, exMC_bkg_etap_3pi0,
  exMC_bkg_g2pi0, exMC_bkg_omega_eta
])