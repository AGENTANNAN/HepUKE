# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) real data at sqrt(s) = 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # matching inclusive MC sample

# Decay card for signal mode I: psi(2S) -> pi+ pi- J/psi, J/psi -> gamma gamma gamma
decay_card_3gamma = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 gamma gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for signal mode II: psi(2S) -> pi+ pi- J/psi, J/psi -> gamma eta_c, eta_c -> gamma gamma
decay_card_gamma_etac = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for each J/psi decay mode, with default cross sections
exMC_3gamma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pipi_jpsi_3gamma"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_3gamma
  config.cross_section   = :default
end

exMC_gamma_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pipi_jpsi_gamma_etac_gg"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_gamma_etac
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Both J/psi decay modes share the identical final state pi+ pi- gamma gamma gamma
# and the same selection criteria, so a single Algorithm instance covers both.
alg_name = "pipiJpsi3gamma"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # ECMS = 3.686 GeV

event_selection = Selection.new
event_selection
    .select_track {                   # Charged track selection
        cos_theta 0.93                # |cos(theta)| < 0.93
        Vz        10.0                # |Vz| < 10 cm
        Vr        1.0                 # Vr < 1 cm
        nChrp     "==1"               # exactly one positively charged track
        nChrn     "==1"               # exactly one negatively charged track
        nNet      "==0"               # net charge zero
    }
    .select_photon {                  # Photon selection
        tdc_emc_start     0           # shower time window start (0*, *0-14 ~ 700 ns units)
        tdc_emc_end       14          # shower time window end
        angle_to_track    5.0         # > 5 degrees from any charged track
        energyThreshold_b 0.025       # E > 25 MeV in the barrel
        energyThreshold_e 0.050       # E > 50 MeV in the endcap
        nGam              ">=3"       # at least three photon candidates
    }
    .assign({:chrgp => :pip, :chrgn => :pim})   # assign the two tracks as pi+ and pi- (no explicit PID)
    .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma]) {  # 4C fit to pi+ pi- gamma gamma gamma
        nominal                    # nominal fit: only its corrected four-momenta are saved
        vertex_fit([0, 1])         # vertex constraint on the pi+ (0) pi- (1) pair
        constrain_four_momentum    # four-momentum (4C) constraint
        chi2_cut 200               # loose chi2 cut; tighter chi2 < 50 applied later in ROOT
    }

# Inexpressible BOSS-side criterion: the photon-multiplicity upper cap.
my_algorithm
    .note(:photon_multiplicity_cap,
          "only events with 3 or 4 photon candidates are retained; nGam expresses the " \
          "minimum (>=3), but the DSL cannot express the <=4 upper bound directly")

my_algorithm.with_decay_card(decay_card_3gamma).apply(event_selection)

# Execute on real data, inclusive MC and both exclusive MC samples
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_3gamma, exMC_gamma_etac])