### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi (3.097 GeV) real data sample
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding inclusive MC sample

# Decay card for the signal process: J/psi -> p pbar a0(980), a0(980) -> pi0 eta, pi0 -> gg, eta -> gg
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 p+ anti-p- a0(980)   PHSP;
    Enddecay

    Decay a0(980)
    1.0000 pi0 eta              PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma          PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma          PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the phase-space sample J/psi -> p pbar pi0 eta (a0(980) resonance removed)
decay_card_phsp = <<~DECAYCARD
    Decay J/psi
    1.0000 p+ anti-p- pi0 eta   PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma          PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma          PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 100k signal events (through a0(980))
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_ppbar_a0980_pi0eta"
  config.related_dataset = jpsi_data
  config.events          = 100_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Exclusive MC: 500k phase-space events (used to derive the mass-dependent efficiency curve)
exMC_phsp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_ppbar_pi0eta_phsp"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_phsp
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "Jpsi2ppbarA0"
my_alg = Algorithm.new(alg_name)
my_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy of the J/psi
       .note(:pid_dedx_only, "p/pbar identification uses the probability method with dE/dx information only; a track is taken as a proton when prob(p) exceeds both prob(K) and prob(pi)")

# Build the event selection chain (stops at the 4C kinematic fit; photon pairing into pi0/eta,
# 3-sigma mass windows, tighter chi2<35 and the pi0pi0 veto are post-fit and applied in ROOT).
event_selection = Selection.new
  .select_track {          # Charged track selection
    cos_theta 0.93         # |cos(theta)| < 0.93
    Vr        1.0          # |Vr| < 1 cm in the transverse plane
    Vz        10.0         # |Vz| < 10 cm along the beam direction
    nChrp     ">=1"        # at least one positive track
    nChrn     ">=1"        # at least one negative track
  }
  .select_photon {         # Photon selection
    tdc_emc_start     0    # EMC timing start (50 ns units)
    tdc_emc_end       14   # EMC timing end (50 ns units)
    energyThreshold_b 0.025 # E > 25 MeV in the barrel
    energyThreshold_e 0.050 # E > 50 MeV in the end caps
    nGam              ">=4" # at least four photons
  }
  .pid(method: :probability) {               # probability-method PID
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion] # p+ and pbar, requires prob(p) > prob(K) and prob(pi)
    nprp ">=1"                                # at least one proton
    nprm ">=1"                                # at least one antiproton
  }
  # 4C kinematic fit to p pbar gamma gamma gamma gamma; the best four-photon combination
  # (smallest chi2) is chosen automatically. Loose chi2<200 here; tighter chi2<35 applied later in ROOT.
  .kinematic_fit([:prp, :prm, :gamma, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

# Generate the algorithm for the signal decay card (same final state p pbar 4gamma as the phase-space MC)
my_alg.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the two exclusive MC samples
root_files = my_alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_phsp])