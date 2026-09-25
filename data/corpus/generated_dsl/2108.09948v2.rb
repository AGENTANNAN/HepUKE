### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding inclusive MC

# Decay card for the signal process (EvtGen format):
# J/psi -> Xi- anti-Xi+ ; Xi- -> Xi0 e- anti-nu_e (semileptonic);
# tag side anti-Xi+ -> anti-Lambda0 pi+, anti-Lambda0 -> anti-p- pi+
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Xi- anti-Xi+ PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 anti-Lambda0 pi+ PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay Xi-
    1.0000 Xi0 e- anti-nu_e PHSP;
    Enddecay

    Decay Xi0
    1.0000 Lambda0 pi0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 1 million events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "jpsi_XimXip_semileptonic"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "XiSemileptonic"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
event_selection
  .select_track {                       # charged track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        100.0                     # |Vz| < 100 cm
    Vr        10.0                      # Vr < 10 mm
    nChrp     ">=3"                     # at least 3 positive tracks
    nChrn     ">=2"                     # at least 2 negative tracks
    nTot      "==5"                     # exactly 5 tracks in total
  }
  .select_photon {                      # photon selection
    tdc_emc_start      0
    tdc_emc_end        14
    angle_to_track     10.0            # > 10 degrees from any charged track
    energyThreshold_b  0.025           # barrel > 25 MeV
    energyThreshold_e  0.050           # endcap > 50 MeV
    nGam               ">=2"           # at least 2 photons
  }
  .pid(method: :probability) {          # PID by probability method, CL > 0.001
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # identify p+ and anti-p- against K/pi
    nprp ">=1"                                  # at least one proton
    nprm ">=1"                                  # at least one anti-proton
    identify :pion, against: [:kaon, :proton]   # identify pi+ and pi- against K/p
    # electron: identified separately (highest e/pi/K confidence from dE/dx and TOF);
    # represented here through the high-momentum lepton identification
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])     # remove (anti-)protons from charged lists
  .kalman_kinematic_fit([:gamma, :gamma]) {     # reconstruct pi0 from two photons (1-C mass fit)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25                                 # chi2 < 25
    npi0 ">=1"                                  # at least one pi0
  }
  .secondary_vertex_fit([:prp, :pim]) {         # Lambda -> p pi- secondary vertex
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {         # anti-Lambda -> anti-p pi+ secondary vertex
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Nominal kinematic fit: anti-Lambda pi+ Lambda pi0 e- with a missing neutrino, 4C constraint
  .kinematic_fit([:Lambda_bar, :pip, :Lambda, :pi0, :lm]) {
    nominal
    constrain_four_momentum
    miss_track_of :nu_e                         # missing anti-nu_e
    chi2_cut 200
    # tag / signal mass windows (mass-difference windows translated to absolute ranges)
    invariant_mass_of(:Lambda_bar).within(1.1107, 1.1207)        # |M(pbar pi+) - M(Lambda)| < 5 MeV
    invariant_mass_of(:Lambda_bar, :pip).within(1.3167, 1.3267)  # |M(Lambdabar pi+) - M(Xi+)| < 5 MeV
    invariant_mass_of(:Lambda, :pi0).within(1.3004, 1.3294)      # |M(Lambda pi0) - M(Xi0)| < 14.5 MeV
  }

my_algorithm
  # electron identification is not the standard momentum-based lepton ID
  .note(:pid_correction_method, "electron selected as the charged track with the highest combined e/pi/K confidence from dE/dx and TOF; represented in the DSL via identify_high_momentum_leptons")
  # single-tag recoil-mass window is not expressible in the BOSS DSL (recoil_mass_of not available)
  .note(:tag_selection_window, "single-tag Xi+ recoil mass required in [1.290, 1.342] GeV (applied as an event window)")
  # signal Xi0 momentum window
  .note(:signal_selection_window, "signal Xi0 momentum required in [0.79, 0.84] GeV")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute the algorithm on real data, inclusive MC and the exclusive signal MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])