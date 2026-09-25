# Core DSL classes and dependencies are loaded automatically at execution

### Dataset preparation ###
data_4416 = DatasetManager.real_data.find("703_4420")        # 4.416 GeV XYZ-scan data point (Ecms = 4415.58 MeV)
incMC_4416 = DatasetManager.inclusive_mc.find("703_4420")    # Corresponding inclusive MC sample

# ---- Signal decay card: e+e- -> pi0 pi0 psi(3686), psi(3686) -> pi+ pi- J/psi, J/psi -> e+ e- ----
decay_card_signal_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi0 pi0 psi(3686)   PHSP;
    Enddecay

    Decay psi(3686)
    1.0 pi+ pi- J/psi       PHSP;
    Enddecay

    Decay J/psi
    1.0 e+ e-               PHOTOS VLL;
    Enddecay

    Decay pi0
    1.0 gamma gamma         PHSP;
    Enddecay

    End
DECAYCARD

# ---- Signal decay card: same topology, J/psi -> mu+ mu- ----
decay_card_signal_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi0 pi0 psi(3686)   PHSP;
    Enddecay

    Decay psi(3686)
    1.0 pi+ pi- J/psi       PHSP;
    Enddecay

    Decay J/psi
    1.0 mu+ mu-             PHOTOS VLL;
    Enddecay

    Decay pi0
    1.0 gamma gamma         PHSP;
    Enddecay

    End
DECAYCARD

# ---- Dominant background: e+e- -> pi+ pi- psi(3686), psi(3686) -> pi0 pi0 J/psi, J/psi -> e+ e- ----
decay_card_background = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- psi(3686)   PHSP;
    Enddecay

    Decay psi(3686)
    1.0 pi0 pi0 J/psi       PHSP;
    Enddecay

    Decay J/psi
    1.0 e+ e-               PHOTOS VLL;
    Enddecay

    Decay pi0
    1.0 gamma gamma         PHSP;
    Enddecay

    End
DECAYCARD

# ---- Exclusive MC samples (100k events each) ----
exMC_signal_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_4416_pi0pi0psip_ee"
  config.related_dataset = data_4416
  config.events          = 100000
  config.decay_card      = decay_card_signal_ee
  config.cross_section   = :default
end

exMC_signal_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_4416_pi0pi0psip_mumu"
  config.related_dataset = data_4416
  config.events          = 100000
  config.decay_card      = decay_card_signal_mumu
  config.cross_section   = :default
end

exMC_background = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_4416_pipipsip_pi0pi0jpsi_ee"
  config.related_dataset = data_4416
  config.events          = 100000
  config.decay_card      = decay_card_background
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "Pi0Pi0Psip"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.416]})   # CMS energy in GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# Common event selection (applies to both signal modes, which share the same final state)
event_selection = Selection.new
event_selection.select_track {              # charged-track quality + event topology
                  cos_theta 0.93            # |cos(theta)| < 0.93
                  Vz        10.0            # |Vz| < 10 cm
                  Vr        1.0             # Vr < 1 cm
                  nChrp     "==2"           # exactly two positive tracks
                  nChrn     "==2"           # exactly two negative tracks
                  nNet      "==0"           # net charge zero
                }
               .select_photon {              # photon (EMC shower) selection
                  tdc_emc_start     0       # EMC timing window: start
                  tdc_emc_end       14      # EMC timing window: end
                  angle_to_track    10.0    # at least 10 deg from any charged track
                  energyThreshold_b 0.025   # > 25 MeV in the barrel
                  energyThreshold_e 0.050   # > 50 MeV in the endcap
                  nGam              ">=4"   # at least four photons
                }
               .pid(method: :probability) {  # PID via per-track probability
                  prob_cut 0.001
                  # tracks with p > 1.0 GeV are treated as leptons;
                  # electron if EMC energy > 0.6 GeV, otherwise muon
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  identify :pion, against: [:kaon]   # pi/K separation
                  npip "==1"
                  npim "==1"
                  nlp  "==1"
                  nlm  "==1"
                }
               # First fit: 4C fit on gamma gamma gamma gamma pi+ pi- l+ l-
               .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :lp, :lm]) {
                  constrain_four_momentum
                  chi2_cut 120
                }
               # Nominal 7C fit: 4C + two pi0 mass constraints + J/psi mass constraint
               .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :lp, :lm]) {
                  nominal
                  constrain_four_momentum
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
                  chi2_cut 200
                }

# Additional BOSS-side requirements that have no dedicated DSL primitive -> captured as notes
my_algorithm
  .note(:muon_counter, "a muon-counter (MUC) hit-depth condition is applied to the l+l- pair")
  .note(:pi0_pairing, "the four photons are paired into the two pi0 by minimizing the summed mass "
                    + "difference, each pair required within 20 MeV/c^2 of the pi0 mass")
  .note(:jpsi_mass_window, "the dilepton invariant mass is required to lie within 3.05-3.15 GeV/c^2")
  .with_decay_card(decay_card_signal_ee).apply(event_selection)

# Execute on real data, inclusive MC and all exclusive MC samples
root_files = my_algorithm.execute_on([data_4416, incMC_4416,
                                      exMC_signal_ee, exMC_signal_mumu, exMC_background])