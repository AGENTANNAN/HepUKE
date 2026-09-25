# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description: 713 R-scan, 22 energy points from 2.00 to 3.08 GeV ###
rscan_names = %w[
  Rscan_2000 Rscan_2050 Rscan_2100 Rscan_2125 Rscan_2150 Rscan_2175 Rscan_2200
  Rscan_2232 Rscan_2309 Rscan_2386 Rscan_2396 Rscan_2500 Rscan_2644 Rscan_2646
  Rscan_2700 Rscan_2800 Rscan_2900 Rscan_2950 Rscan_2981 Rscan_3000 Rscan_3020
  Rscan_3080
]
# 713 real data and the matching inclusive MC, one dataset per scan point
rscan_data  = rscan_names.map { |n| DatasetManager.real_data.find("713_#{n}") }
rscan_incMC = rscan_names.map { |n| DatasetManager.inclusive_mc.find("713_#{n}") }

# ConExc decay card: e+e- -> phi eta (ConExc mode 23), phi -> K+K-, eta -> gamma gamma.
# `Particle vpho` is omitted: the DSL injects it per energy point for a multi-point scan,
# and switches to the no-KKMC simulation template (auto-detected from the `ConExc` token).
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 23;
    Enddecay
    End
DECAYCARD

# 2.5M-event ConExc exclusive signal MC generated at EVERY energy point
# (one ExclusiveMC per data set, sharing the decay card / cross section).
exMCs_signal = DatasetManager.create_exclusive_mc_for(rscan_data) do |config|
  config.sample_name   = "sig_phi_eta_conexc"
  config.events        = 2_500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PhiEtaScan"
phi_eta = Algorithm.new(alg_name)
phi_eta.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.080]})  # representative; per-point ECMS injected by the scan
       .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  .select_track {                  # exactly two charged tracks of opposite charge
     cos_theta 0.93                # |cos(theta)| < 0.93
     Vz        10.0                # |Vz| < 10 cm
     Vr        1.0                 # Vr < 1 cm
     nChrp     "==1"               # exactly one positive track
     nChrn     "==1"               # exactly one negative track
     nNet      "==0"               # net charge zero
   }
  .select_photon {                 # at least two good photons
     tdc_emc_start     0           # EMC time within [0, 700] ns of the event start
     tdc_emc_end       14
     energyThreshold_b 0.070       # E > 70 MeV in the barrel
     energyThreshold_e 0.070       # E > 70 MeV in the endcap
     nGam              ">=2"
   }
  .pid(method: :probability) {     # identify both tracks as kaons against pi+/pi-
     prob_cut 0.001
     identify :kaon, against: [:pion]
     nkp "==1"
     nkm "==1"
   }
   # Nominal 4C kinematic fit to K+K-gamma-gamma with chi2 cut 100. When more than two
   # good photons are found the fitter automatically keeps the gamma-gamma pairing
   # (and kaon assignment) giving the smallest chi2.
  .kinematic_fit([:kp, :km, :gamma, :gamma]) {
     nominal
     constrain_four_momentum
     chi2_cut 100
   }
   # Competing ISR-phi hypothesis (no chi2_cut, no nominal): stores chi2(K+K-gamma) so the
   # veto chi2(K+K-gamma) < chi2(K+K-gamma gamma) can be applied in the ROOT analysis.
  .kinematic_fit([:kp, :km, :gamma]) {
     constrain_four_momentum
   }

# Attach the decay card and run over real data, inclusive MC and the signal MC.
phi_eta.with_decay_card(decay_card_signal).apply(event_selection)
root_files = phi_eta.execute_on(rscan_data + rscan_incMC + exMCs_signal)