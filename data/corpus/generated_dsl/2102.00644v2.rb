# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###

# --- Real data: energy points of the psi(4260) scan (4.178 - 4.600 GeV), BOSS 7.0.3 ---
data_points = [
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4245"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("703_4310"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("703_4390"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("703_4470"),
  DatasetManager.real_data.find("703_4530"),
  DatasetManager.real_data.find("703_4575"),
  DatasetManager.real_data.find("703_4600")
]

# --- Corresponding inclusive MC samples (same energy points) ---
incmc_points = [
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230"),
  DatasetManager.inclusive_mc.find("703_4237"),
  DatasetManager.inclusive_mc.find("703_4246"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4270"),
  DatasetManager.inclusive_mc.find("703_4280"),
  DatasetManager.inclusive_mc.find("703_4360"),
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("703_4600")
]

# Decay card (EvtGen) for the signal process:
#   e+e- -> pi0 X(3872) gamma,  X(3872) -> pi+ pi- J/psi,
#   J/psi -> e+ e-,  pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000 pi0 X(3872) gamma PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi+ pi- J/psi PHSP;
  Enddecay

  Decay J/psi
  1.000 e+ e- PHOTOS VLL;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# One 100k-event exclusive MC sample per energy point (same signal mode / decay card)
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_pi0_x3872_gamma"
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "Pi0X3872Gamma"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.26]})

# Helix-parameter correction applied to charged tracks before the kinematic fit
my_algorithm.note(:helix_correction,
  "helix-parameter correction applied to all charged tracks before the 4C kinematic fit; efficiency difference between with/without correction estimated by re-running the BOSS selection on signal MC")

# Build the event selection chain (all functions called in chain)
event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta 0.93        # |cos(theta)| < 0.93
                  Vz        10.0        # |Vz| < 10 cm
                  Vr        1.0         # Vr < 1 cm in the transverse plane
                  nChrp     ">=2"       # at least two positive tracks
                  nChrn     ">=2"       # at least two negative tracks
                  nNet      "==0"       # net charge zero
                }
               .select_photon {         # Photon selection
                  tdc_emc_start     0     # TDC start (0 - 700 ns window)
                  tdc_emc_end       14
                  angle_to_track    10.0  # min angle to nearest charged track (degrees)
                  energyThreshold_b 0.025 # barrel energy threshold (GeV)
                  energyThreshold_e 0.050 # endcap energy threshold (GeV)
                  nGam              ">=3" # at least three photons
                }
               .pid(method: :probability) {
                  prob_cut 0.001        # PID probability > 0.001
                  # high-momentum tracks (p > 1.0) treated as leptons;
                  # lepton with EMC energy > 0.6 -> electron, else muon
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  nlp "==1"             # exactly one l+
                  nlm "==1"             # exactly one l-
                }
               .remove([:lp <= :chrgp, :lm <= :chrgn])   # lepton lists populated by the lepton PID
               .assign({:chrgp => :pip, :chrgn => :pim}) # remaining tracks treated as pi+/pi-
               .kalman_kinematic_fit([:gamma, :gamma]) { # reconstruct pi0 from photon pairs
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25           # chi2 < 25
                  npi0     ">=1"        # at least one pi0 candidate
                }
               .kinematic_fit([:lp, :lm, :pip, :pim, :pi0, :gamma]) {
                  nominal               # nominal fit (four-momenta taken from here)
                  constrain_four_momentum # 4C energy-momentum constraint
                  chi2_cut 60           # chi2 < 60
                }

# Generate the algorithm for the signal process in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC, and the exclusive MC signal samples
root_files = my_algorithm.execute_on(data_points + incmc_points + exMCs_signal)