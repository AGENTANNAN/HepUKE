# Core DSL classes and dependencies are loaded automatically at execution

### Dataset description ###
# Eleven real data points spanning 4.178 - 4.280 GeV (energy scan)
scan_names = ["703_4180", "703_4190", "703_4200", "703_4210", "703_4220",
              "703_4230", "703_4237", "703_4246", "703_4260", "703_4270", "703_4280"]
data_points  = scan_names.map { |n| DatasetManager.real_data.find(n) }   # real data at each point
incMC_points = scan_names.map { |n| DatasetManager.inclusive_mc.find(n) } # matching inclusive MC

# Decay card: e+e- -> gamma chi_c1(3872), chi_c1(3872) -> gamma psi_2(3823),
#             psi_2(3823) -> gamma chi_c1, chi_c1 -> gamma J/psi, J/psi -> e+ e-
signal_decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000 gamma chi_c1(3872) PHSP;
    Enddecay

    Decay chi_c1(3872)
    1.000 gamma psi_2(3823) PHSP;
    Enddecay

    Decay psi_2(3823)
    1.000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# 200k-event signal exclusive MC, generated for every scan point
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_gamchic1psi2_cascade"
  config.events        = 200000
  config.decay_card    = signal_decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "GamChiC1Psi2Cascade"
my_algorithm = Algorithm.new(alg_name)
my_algorithm
  .set_header(["#{alg_name}Alg/#{alg_name}.h"])
  .set_constant({"ECMS" => [:double, 4.260]})
  .note(:ecms_scan, "11-point energy scan 4.178-4.280 GeV; ECMS is set to a representative value (4.260 GeV) here, but the per-run beam energy of each scan point must be used in the actual kinematic fit")
  .note(:pid_correction_method, "lepton identification uses EMC-based separation: electron if EMC deposit >= 0.8 GeV, muon if deposit < 0.4 GeV; the asymmetric window is not fully captured by identify_high_momentum_leptons (single electron energy threshold), which instead treats tracks above the momentum threshold as leptons and uses treat_as_electron_if_energy_above for the e/mu split")
  .note(:decay_model, "the initial radiative production e+e- -> gamma chi_c1(3872) is an E1 transition; modelled with the PHSP generator in the decay card")

event_selection = Selection.new
event_selection
  .select_track {                 # exactly two charged tracks, one positive and one negative
        cos_theta   0.93            # |cos(theta)| < 0.93
        Vz          10.0            # |Vz| < 10 cm
        Vr          1.0             # Vr < 1 cm
        nChrp       "==1"           # 1 positive track
        nChrn       "==1"           # 1 negative track
        nNet        "==0"           # net charge zero
  }
  .select_photon {                # at least four good photons
        tdc_emc_start     0         # EMC timing 0 ...
        tdc_emc_end       14        # ... to 700 ns
        angle_to_track    10.0      # > 10 degrees from any charged track
        energyThreshold_b 0.025     # 25 MeV threshold in the barrel
        energyThreshold_e 0.025     # 25 MeV threshold in the endcap
        nGam              ">=4"     # at least 4 photons
  }
  .pid(method: :probability) {    # lepton id (e/mu) via probability method, separated from pions
        prob_cut 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.8
        nlp "==1"                   # one lepton of positive charge
        nlm "==1"                   # one lepton of negative charge
  }
  # nominal 4C fit: l+l- + four photons to the beam four-momentum (best combination by minimum chi2)
  .kinematic_fit([:lp, :lm, :gamma, :gamma, :gamma, :gamma]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:lp, :lm).within(3.067, 3.127)        # m(l+l-) within 30 MeV of J/psi
        invariant_mass_of(:gamma, :gamma).out_of(0.120, 0.150)  # veto any photon pair at the pi0 mass
        chi2_cut 200
  }
  # 7C fit: adds the J/psi, chi_c1 and psi_2(3823) mass constraints to resolve the photon assignment
  .kinematic_fit([:lp, :lm, :gamma, :gamma, :gamma, :gamma]) {
        constrain_four_momentum
        invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
        invariant_mass_of(:lp, :lm, :gamma).constrain_to_nominal_mass_of(:chi_c1)
        invariant_mass_of(:lp, :lm, :gamma, :gamma).constrain_to_nominal_mass_of(:psi_2)
        chi2_cut 100
  }

my_algorithm.with_decay_card(signal_decay_card).apply(event_selection)
root_files = my_algorithm.execute_on(data_points + incMC_points + exMCs_signal)