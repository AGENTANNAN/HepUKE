# ============================================================
# BESIII R-scan (BOSS 713): e+e- -> K+K-pi0pi0 for a partial-wave analysis
# Ten energy points: 2.000, 2.100, 2.125, 2.175, 2.200, 2.232,
#                    2.309, 2.386, 2.396, 2.644 GeV  (total ~300 pb^-1)
# ============================================================

### Dataset preparation ###
# Real R-scan data at the ten energy points (BOSS version 713, sample names Rscan_*)
data_points = [
  DatasetManager.real_data.find("713_Rscan_2000"),
  DatasetManager.real_data.find("713_Rscan_2100"),
  DatasetManager.real_data.find("713_Rscan_2125"),
  DatasetManager.real_data.find("713_Rscan_2175"),
  DatasetManager.real_data.find("713_Rscan_2200"),
  DatasetManager.real_data.find("713_Rscan_2232"),
  DatasetManager.real_data.find("713_Rscan_2309"),
  DatasetManager.real_data.find("713_Rscan_2386"),
  DatasetManager.real_data.find("713_Rscan_2396"),
  DatasetManager.real_data.find("713_Rscan_2644"),
]

# ConExc decay card: ISR production of K+K-pi0pi0 (continuum / R-scan Born cross section).
# The DSL detects the literal token "ConExc", switches to the no-KKMC template and injects
# `Particle vpho <ECMS> 0.0` per energy point — so `Particle vpho` must be omitted here.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 11;
    Enddecay
    End
DECAYCARD

# 200,000-event signal MC for EACH energy point (shared card, one ExclusiveMC per dataset)
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "KKpi0pi0_signal_mc"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "KKpi0pi0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 2.125]})  # nominal; per-point ECMS is injected for the scan

event_selection = Selection.new
event_selection.select_track {          # charged track selection
                  cos_theta 0.93         # |cos(theta)| < 0.93
                  Vz        10.0         # |Vz| < 10 cm
                  Vr        1.0          # Vr < 1 cm (x-y transverse plane)
                  nChrp     ">=1"        # at least one positive track
                  nChrn     ">=1"        # at least one negative track
                  nNet      "==0"        # net charge zero
                }
               .select_photon {          # photon selection
                  tdc_emc_start     0    # EMC TDC window start
                  tdc_emc_end       14   # EMC TDC window end
                  energyThreshold_b 0.025 # barrel energy threshold (25 MeV)
                  energyThreshold_e 0.050 # endcap energy threshold (50 MeV)
                  angle_to_track    10.0  # angle to nearest charged track > 10 deg
                  nGam              ">=4" # at least four photons
                }
               .pid(method: :probability) {   # kaon identification, probability method
                  prob_cut 0.001
                  identify :kaon, against: [:pion]  # K+ and K- (charge-conjugation shorthand), separated from pions
                  nkp "==1"               # exactly one K+
                  nkm "==1"               # exactly one K-
                }
               # Reconstruct two pi0 from photon pairs (1C Kalman fit, gamma-gamma -> pi0 mass)
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 ">=2"              # require at least two pi0
                }
               # Main 6C kinematic fit to K+K-pi0pi0 (built from four photons):
               # 4-momentum conservation + the two pi0 mass constraints = 6C
               .kinematic_fit([:kp, :km, :pi0, :pi0]) {
                  nominal                 # nominal fit — corrected four-momenta come from here
                  constrain_four_momentum
                  chi2_cut 80
                }

# Generate the BOSS algorithm for the signal decay card and render the selection chain
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Run over the ten data points and their matching signal MC samples;
# the surviving events (ntuple) are passed on to the GPUPWA partial-wave analysis
root_files = my_algorithm.execute_on(data_points + exMCs_signal)