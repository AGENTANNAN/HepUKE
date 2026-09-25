# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# J/psi (3.097 GeV) real data and its inclusive MC (2x10^8 J/psi)
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi real data sample
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # corresponding inclusive MC

# Decay card for the signal process: J/psi -> omega eta pi+ pi- (omega -> pi+ pi- pi0 Dalitz chain)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 omega eta pi+ pi- PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the phase-space MC: same final state, no omega resonance
# (used to build the phase-space correction curve)
decay_card_phsp = <<~DECAYCARD
    Decay J/psi
    1.000 pi+ pi- pi+ pi- pi0 eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the peaking b_1(1235)a_0(980) background (b_1 -> omega pi, a_0 -> eta pi)
decay_card_bkg = <<~DECAYCARD
    Decay J/psi
    1.000 b_1(1235)+ a_0(980)- PHSP;
    Enddecay

    Decay b_1(1235)+
    1.000 omega pi+ PHSP;
    Enddecay

    Decay a_0(980)-
    1.000 eta pi- PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive signal MC (omega eta pi+ pi- Dalitz chain)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_omega_eta_pipi"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# 2M-event phase-space MC of the same final state without the omega resonance
exMC_phsp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_omega_eta_pipi_phsp"
  config.related_dataset = jpsi_data
  config.events          = 2000000
  config.decay_card      = decay_card_phsp
  config.cross_section   = :default
end

# 100k-event exclusive MC for the peaking b_1(1235)a_0(980) background
exMC_bkg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_b1a0_peaking_bkg"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_bkg
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "omegaetapipi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy = 3.097 GeV

event_selection = Selection.new
event_selection
  .select_track {                     # four charged tracks, 2 positive + 2 negative, net charge 0
    cos_theta   0.93                  # |cos(theta)| < 0.93
    Vz          20.0                  # |Vz| < 20 cm
    Vr          2.0                   # Vr < 2 cm
    nChrp       "==2"                 # exactly 2 positive tracks
    nChrn       "==2"                 # exactly 2 negative tracks
    nNet        "==0"                 # net charge zero
  }
  .select_photon {                    # at least four good photons
    tdc_emc_start     0               # TDC window 0-14
    tdc_emc_end       14
    angle_to_track    10.0            # at least 10 degrees from any charged track
    energyThreshold_b 0.025           # EMC energy > 25 MeV (barrel)
    energyThreshold_e 0.050           # EMC energy > 50 MeV (endcap)
    nGam              ">=4"           # at least 4 photons
  }
  .assign({:chrgp => :pip, :chrgn => :pim})   # no PID: all charged tracks treated as pions
  .kalman_kinematic_fit([:gamma, :gamma]) {   # Kalman fit: one gamma-gamma pair -> pi0 mass
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"                        # require at least one pi0
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # Kalman fit: another gamma-gamma pair -> eta mass
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"                        # require at least one eta
  }
  .kinematic_fit([:pip, :pim, :pip, :pim, :pi0, :eta]) {   # nominal 4C fit on pi+pi-pi+pi-pi0 eta
    nominal
    vertex_fit([0, 1, 2, 3])          # common-vertex constraint on the four charged tracks
    constrain_four_momentum           # 4-momentum conservation
    chi2_cut 50                       # chi2_4C < 50
  }

# The IP constraint on the four charged tracks and the vertex-fit chi2 (chi2_V < 100)
# have no dedicated DSL expression; capture them for the systematics skill.
my_algorithm
  .note(:vertex_fit_constraint, "the four charged tracks are constrained to a common
    vertex that is also constrained to the interaction point (IP); the vertex-fit
    chi2 is required to be < 100 (chi2_V < 100) - the IP constraint and the vertex
    chi2 value are not expressible in the DSL")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC and the three exclusive MC samples
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_phsp, exMC_bkg])