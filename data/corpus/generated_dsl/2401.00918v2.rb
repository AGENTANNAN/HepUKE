### Dataset preparation ###
# J/psi real data and inclusive MC at sqrt(s) = 3.097 GeV
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal process J/psi -> gamma gamma phi, phi -> K+ K-
# PHSP generator models the gamma gamma phi system for the partial-wave analysis
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 gamma gamma phi        PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# 1M-event exclusive MC for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_ggphi_kk"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "GGPhiKK"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {
    cos_theta 0.93     # |cos(theta)| < 0.93
    Vz        10.0     # |Vz| < 10 cm
    Vr        1.0      # Vr < 1 cm
    nChrp     "==1"    # exactly one positive track
    nChrn     "==1"    # exactly one negative track
    nNet      "==0"    # net charge zero
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025   # EMC barrel energy > 25 MeV
    energyThreshold_e 0.050   # EMC endcap energy > 50 MeV
    nGam              ">=2"   # at least two photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]   # K+ and K- (charge-conjugation shorthand)
    nkp "==1"
    nkm "==1"
  }
  # Nominal 4C kinematic fit to gamma gamma K+ K-
  .kinematic_fit([:gamma, :gamma, :kp, :km]) {
    nominal
    constrain_four_momentum
    chi2_cut 40
  }
  # Competing 4C hypotheses (no nominal, no chi2_cut) — chi2 values are stored
  # and the "nominal chi2 < all competing chi2" veto is applied later in ROOT
  .kinematic_fit([:gamma, :kp, :km]) {
    constrain_four_momentum
  }
  .kinematic_fit([:gamma, :gamma, :gamma, :kp, :km]) {
    constrain_four_momentum
  }
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :kp, :km]) {
    constrain_four_momentum
  }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])