### Dataset selection ###
# Real data at 18 energy points from 4.288 to 4.951 GeV (BOSS 705 / 706 / 707)
data_points = [
  DatasetManager.real_data.find("705_4290"),
  DatasetManager.real_data.find("705_4315"),
  DatasetManager.real_data.find("705_4340"),
  DatasetManager.real_data.find("705_4380"),
  DatasetManager.real_data.find("705_4400"),
  DatasetManager.real_data.find("705_4440"),
  DatasetManager.real_data.find("706_4610"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4740"),
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
  DatasetManager.real_data.find("707_4840"),
  DatasetManager.real_data.find("707_4914"),
  DatasetManager.real_data.find("707_4946")
]

# Inclusive MC samples available in the same energy range (706 and 707 points)
incMC_points = [
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
  DatasetManager.inclusive_mc.find("707_4740"),
  DatasetManager.inclusive_mc.find("707_4750"),
  DatasetManager.inclusive_mc.find("707_4780"),
  DatasetManager.inclusive_mc.find("707_4840"),
  DatasetManager.inclusive_mc.find("707_4914"),
  DatasetManager.inclusive_mc.find("707_4946")
]

### Decay cards (EvtGen format) ###
# Signal A: e+e- -> eta psi(2S) (HELAMP production), psi(2S) -> pi+ pi- J/psi, J/psi -> e+ e-, eta -> gamma gamma
decay_card_eta_psi2s = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta psi(2S) HELAMP 1.0 0.0 1.0 0.0 1.0 0.0;
  Enddecay

  Decay psi(2S)
  1.0000 pi+ pi- J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Signal B: e+e- -> eta X(3872) (PHSP production), X(3872) -> pi+ pi- J/psi, J/psi -> e+ e-, eta -> gamma gamma
decay_card_eta_X = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.0000 pi+ pi- J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive signal MC: 200k-event sample per energy point for each signal ###
exMC_eta_psi2s = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_eta_psi2s_ee"
  config.events        = 200000
  config.decay_card    = decay_card_eta_psi2s
  config.cross_section = :default
end

exMC_eta_X = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_eta_X3872_ee"
  config.events        = 200000
  config.decay_card    = decay_card_eta_X
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "EtaPsi2sLepton"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.6]})   # representative CMS energy of the 4.288-4.951 GeV scan
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93      # |cos(theta)| < 0.93
                  Vz        100.0     # |Vz| < 100 mm along the beam direction
                  Vr        10.0      # Vr < 10 mm in the transverse plane
                  nChrp     "==2"     # exactly two positive charged tracks
                  nChrn     "==2"     # exactly two negative charged tracks
                  nNet      "==0"     # net charge zero
                }
               .select_photon {
                  tdc_emc_start     0      # EMC timing window 0-14 (units of 700 ns)
                  tdc_emc_end       14
                  angle_to_track    10.0   # > 10 degrees from the nearest charged track
                  energyThreshold_b 0.025  # 25 MeV in the barrel
                  energyThreshold_e 0.050  # 50 MeV in the endcap
                  nGam              ">=2"  # at least two photons (eta -> gamma gamma)
                }
               .pid(method: :probability) {
                  prob_cut 0.001  # PID probability > 0.001
                  # Tracks with p > 1.0 GeV/c treated as leptons; electron if EMC energy > 0.6 GeV, else muon
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  identify :pion, against: [:kaon]  # pi/K separation
                  npip "==1"   # exactly one pi+
                  npim "==1"   # exactly one pi-
                  nlp  "==1"   # exactly one l+
                  nlm  "==1"   # exactly one l-
                }
               # Nominal 4C fit to gamma gamma pi+ pi- e+ e- (l+l- under the e hypothesis)
               .kinematic_fit([:gamma, :gamma, :pip, :pim, :lp, :lm]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 40
                }
               # Competing 4C fit: same tracks under the mu+ mu- hypothesis; no chi2 cut,
               # the chi2 is stored for the ROOT-level choice of the lepton-pair hypothesis
               .assign({:lp => :mup, :lm => :mum})
               .kinematic_fit([:gamma, :gamma, :pip, :pim, :mup, :mum]) {
                  constrain_four_momentum
                }

# Both signals (eta psi(2S) and eta X(3872)) share the same final state and selection,
# so a single Algorithm instance is used.
my_algorithm.with_decay_card(decay_card_eta_psi2s).apply(event_selection)

root_files = my_algorithm.execute_on(data_points + incMC_points + exMC_eta_psi2s + exMC_eta_X)