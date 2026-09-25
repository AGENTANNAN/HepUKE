# Core DSL classes/dependencies are loaded automatically at execution.

### Dataset preparation ###
# e+e- -> phi pi+ pi- (phi -> K+ K-), measured at 22 R-scan points (BOSS 713), 2.00 - 3.08 GeV.
rscan_points = %w[
  713_Rscan_2000 713_Rscan_2050 713_Rscan_2100 713_Rscan_2125 713_Rscan_2150
  713_Rscan_2175 713_Rscan_2200 713_Rscan_2232 713_Rscan_2309 713_Rscan_2386
  713_Rscan_2396 713_Rscan_2500 713_Rscan_2644 713_Rscan_2646 713_Rscan_2700
  713_Rscan_2800 713_Rscan_2900 713_Rscan_2950 713_Rscan_2981 713_Rscan_3000
  713_Rscan_3020 713_Rscan_3080
]
rscan_data  = rscan_points.map { |name| DatasetManager.real_data.find(name) }     # 22 real-data points
rscan_incMC = rscan_points.map { |name| DatasetManager.inclusive_mc.find(name) } # matching inclusive-MC points

# Decay card for the continuum signal, generated with the ConExc model (ISR up to second
# order). Mode 35 corresponds to gamma* -> pi+ pi- K+ K- in phase space.
# For a multi-energy scan the DSL injects "Particle vpho <ECMS> 0.0" per point, so it is omitted.
decay_card_signal = <<~DECAYCARD
  Decay vpho
  1.0 ConExc 35;
  Enddecay

  End
DECAYCARD

# One signal exclusive-MC sample per energy point, 1M events each (same channel / decay card).
exMCs_signal = DatasetManager.create_exclusive_mc_for(rscan_data) do |config|
  config.sample_name   = "sig_phi_pipi_scan"   # auto-suffixed per energy point
  config.events        = 1_000_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PhiPiPiScan"
phi_pipi_alg = Algorithm.new(alg_name)
phi_pipi_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.080]})  # per-point energy driven by the scan datasets

event_selection = Selection.new
event_selection
  .select_track {                       # charged track selection (no photon selection in this analysis)
    cos_theta 0.93                       # |cos(theta)| < 0.93
    Vz        10.0                       # |Vz| < 10 cm
    Vr        1.0                        # Vr < 1 cm
    nTot      ">=3"                      # at least three charged tracks
  }
  .pid(method: :probability) {           # probability-method PID
    prob_cut 0.001                                             # PID probability > 0.001
    identify :pion, against: [:kaon, :proton]                  # separate pi from K/p (pi+ and pi-)
    identify :kaon, against: [:pion, :proton]                  # identify K+ and K- (K- may be missing)
    npip ">=1"                                                 # at least one pi+
    npim ">=1"                                                 # at least one pi-
    nkp  ">=1"                                                 # at least one K+
  }
  # Nominal missing-kaon hypothesis: pi+ pi- K+ measured, K- missing.
  # Vertex fit on pi+ pi- K+ (indices 0,1,2), then four-momentum (4C) fit; chi2 < 10.
  .kinematic_fit([:pip, :pim, :kp, :km]) {
    nominal
    miss_track_of :km                                          # K- defined as the missing track
    vertex_fit([0, 1, 2])                                      # common vertex for pi+, pi-, K+
    constrain_four_momentum
    chi2_cut 10
  }
  # Charge-conjugate pi+ pi- K- hypothesis (K+ missing); chi2 stored so the missing-kaon
  # hypothesis with the smaller chi2 can be kept downstream.
  .kinematic_fit([:pip, :pim, :km, :kp]) {
    miss_track_of :kp
    vertex_fit([0, 1, 2])
    constrain_four_momentum
  }

# Generate the BOSS algorithm from the signal decay card and the selection chain.
phi_pipi_alg.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on the 22 real-data points, their inclusive MC, and the per-point signal MC.
root_files = phi_pipi_alg.execute_on(rscan_data + rscan_incMC + exMCs_signal)