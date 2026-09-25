# frozen_string_literal: true
### Dataset preparation ###
# 2015 R-scan (BOSS 713) real-data points inside the 2.00-3.08 GeV range.
# The two ~1 pb^-1 scan points (2.70, 2.80 GeV) are not usable for a
# cross-section measurement and are dropped, leaving 19 c.m. energy points.
rscan_names = %w[
  713_Rscan_2000 713_Rscan_2050 713_Rscan_2100 713_Rscan_2150 713_Rscan_2175
  713_Rscan_2200 713_Rscan_2232 713_Rscan_2309 713_Rscan_2386 713_Rscan_2396
  713_Rscan_2500 713_Rscan_2644 713_Rscan_2646 713_Rscan_2900 713_Rscan_2950
  713_Rscan_2981 713_Rscan_3000 713_Rscan_3020 713_Rscan_3080
]
rscan_data = rscan_names.map { |name| DatasetManager.real_data.find(name) }

# ConExc decay card for the continuum process e+e- -> pi+ pi- pi0 (mode 7).
# ConExc models ISR up to second order and the measured sigma_0(m); the DSL
# auto-detects the ConExc token and injects `Particle vpho <ECMS> 0.0` per
# energy point, so no vpho particle line is written here.
decay_card_conexc = <<~DECAYCARD
    Decay vpho
    1.0000 ConExc 7;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# One 50k-event exclusive MC sample per R-scan energy point, all sharing the
# same ConExc decay card and cross-section treatment.
exMC_signal = DatasetManager.create_exclusive_mc_for(rscan_data) do |config|
  config.sample_name   = "exmc_rscan_pipipimpi0"  # auto-suffixed per energy point
  config.events        = 50000
  config.decay_card    = decay_card_conexc
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name  = "RscanPipPimPi0"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 2.0]})  # default only; the real c.m. energy is set per run
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:beam_energy, "the c.m. energy is taken per run from the measured CM energy
            (beam-energy conditions DB) instead of a fixed ECMS value, since the analysis
            uses the 2.00-3.08 GeV R-scan data points")

event_selection = Selection.new
event_selection.select_track {    # charged track selection
                  cos_theta   0.93   # |cos(theta)| < 0.93, theta = polar angle
                  Vz          10.0   # |Vz| < 10 cm
                  Vr          1.0    # Vr < 1 cm in the transverse plane
                  nChrp       ">=1"  # at least one positive track
                  nChrn       ">=1"  # at least one negative track
                  nNet        "==0"  # net charge summed over all tracks = 0
                }
               .select_photon {   # photon selection
                  tdc_emc_start     0     # TDC start time
                  tdc_emc_end       14    # TDC end time
                  energyThreshold_b 0.025 # barrel energy threshold: 25 MeV
                  energyThreshold_e 0.050 # endcap energy threshold: 50 MeV
                  angle_to_track    10.0  # photon-track angle > 10 degrees
                  nGam              ">=2" # at least two good photons
                }
               .pid(method: :probability) {  # pion identification (hadrons only)
                  prob_cut 0.001
                  identify :pion, against: [:kaon, :proton]  # pi+ and pi- vs K and p (charge-conj. shorthand)
                  npip ">=1"
                  npim ">=1"
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {  # 1-C mass-constrained fit: pi0 -> gamma gamma
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 ">=1"    # at least one pi0 candidate
                }
               .kinematic_fit([:pip, :pim, :pi0]) {   # 4C fit of the final state pi+ pi- pi0
                  nominal               # nominal fit: its corrected four-momenta are used
                  constrain_four_momentum
                  chi2_cut 50
                }

# Generate the BOSS algorithm for the ConExc signal process and run it on every
# R-scan data point together with its dedicated exclusive MC sample.
algorithm.with_decay_card(decay_card_conexc).apply(event_selection)
root_files = algorithm.execute_on(rscan_data + exMC_signal)