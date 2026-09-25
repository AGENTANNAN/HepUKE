# J/psi -> K+ K- pi0 (pi0 -> gamma gamma) — partial-wave analysis of the K+ K- pi0 final state
# BOSS part: dataset preparation + event selection up to and including the 5C kinematic fit.

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # 2009 J/psi real data @ 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # corresponding inclusive MC sample

# Decay card for the signal process (EvtGen format); exclusive phase-space MC
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 K+ K- pi0   PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Create the exclusive signal MC sample (1,000,000 phase-space events, J/psi -> K+ K- pi0, pi0 -> gamma gamma)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_kkpi0"
  config.related_dataset = jpsi_data      # associated real dataset
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "JpsiKKpi0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93   # |cos(theta)| < 0.93
                  Vz        10.0   # |Vz| < 10 cm
                  Vr        1.0    # Vr < 1 cm (transverse plane)
                  nChrp     ">=1"  # at least one positively charged track
                  nChrn     ">=1"  # at least one negatively charged track
                  nNet      "==0"  # net charge zero
                }
               .select_photon {
                  tdc_emc_start     0      # EMC timing window 0-700 ns
                  tdc_emc_end       14
                  angle_to_track    10.0   # more than 10 deg from any charged track
                  energyThreshold_b 0.025  # E > 25 MeV in the barrel region
                  energyThreshold_e 0.050  # E > 50 MeV in the endcap region
                  nGam              ">=2"  # at least two photons
                }
               .pid(method: :probability) {
                  prob_cut 0.001                         # probability cut 0.001
                  identify :kaon, against: [:pion, :proton]  # assign each track among pi, K, p; K+ and K-
                  nkp "==1"                              # exactly one K+
                  nkm "==1"                              # exactly one K-
                }
               # 5C kinematic fit: 4C four-momentum conservation + 1C pi0 mass constraint on the gamma gamma pair
               .kinematic_fit([:kp, :km, :gamma, :gamma]) {
                  nominal
                  constrain_four_momentum
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 60
                }

my_algorithm
  .note(:track_pt_cut, "any charged track with transverse momentum pT < 120 MeV/c is rejected")
  .note(:pid_correction_method, "full PID uses the highest combined dE/dx + TOF confidence level")
  .note(:photon_pair_preselection, "diphoton pairs are preselected with M(gamma gamma) < 300 MeV/c^2")
  .note(:pi0_mass_window, "pi0 candidates must satisfy 110 < M(gamma gamma) < 150 MeV/c^2 after the 4C preselection")
  .note(:background_veto, "events are rejected if a background hypothesis (gamma gamma pi+ pi-, gamma K+ K-, or gamma gamma gamma K+ K-) gives a lower chi^2 than the signal hypothesis")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute the algorithm on the real data, inclusive MC and exclusive signal MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])