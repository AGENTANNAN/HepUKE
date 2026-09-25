# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# J/psi (3.097 GeV) real data and corresponding inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal process (EvtGen format, EvtGen particle names):
# J/psi -> gamma eta', eta' -> eta pi0 pi0, eta -> gamma gamma, pi0 -> gamma gamma
# (all-neutral 7-photon final state)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'        PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi0 pi0       PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma       PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma       PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive signal MC for the same decay chain, tied to the J/psi dataset
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gamma_etap_eta_2pi0"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name  = "JpsiToGammaEtaPrime"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})

# Full selection chain, applied identically to data, inclusive MC and exclusive MC
event_selection = Selection.new
event_selection.select_track {          # all-neutral final state: reject any charged track
                  nChrp  "==0"          # zero positive charged tracks
                  nChrn  "==0"          # zero negative charged tracks
                  nNet   "==0"          # net charge zero
                }
               .select_photon {          # EMC barrel |cos(theta)|<0.80 / endcap 0.86<|cos(theta)|<0.92 handled automatically
                  tdc_emc_start    0     # TDC start time
                  tdc_emc_end      14    # TDC end time
                  angle_to_track   10.0  # min angle to nearest charged track (degrees)
                  energyThreshold_b 0.025 # EMC barrel energy threshold (25 MeV)
                  energyThreshold_e 0.050 # EMC endcap energy threshold (50 MeV)
                  nGam   ">=7"           # at least seven photons
                }
                # 1C Kalman fit: reconstruct one eta -> gamma gamma (chi2 < 25)
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta) # nominal eta mass constraint
                  chi2_cut 25
                  neta ">=1"             # at least one eta candidate
                }
                # additional 1C Kalman fits: reconstruct at least two pi0 -> gamma gamma (chi2 < 25)
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0) # nominal pi0 mass constraint
                  chi2_cut 25
                  npi0 ">=2"             # at least two pi0 candidates
                }
                # 8C kinematic fit on gamma eta pi0 pi0: 4C four-momentum conservation plus the
                # eta, both pi0 (via the 1C Kalman fits above) and the eta' mass constraints.
                # If several combinations survive, the smallest-chi2 one is kept automatically.
               .kinematic_fit([:gamma, :eta, :pi0, :pi0]) {
                  nominal                                                        # nominal fit: corrected four-momenta are saved
                  constrain_four_momentum                                        # 4C energy-momentum conservation
                  invariant_mass_of(:eta, :pi0, :pi0).constrain_to_nominal_mass_of(:etap) # eta' mass (1C)
                  chi2_cut 100                                                   # chi2 < 100
                }

# BOSS-side criteria that cannot be expressed with the current DSL constructs
algorithm
  .note(:radiative_photon, "the radiative photon from J/psi -> gamma eta' is taken as the most energetic photon in the event")
  .note(:pi0_cos_theta, "each pi0 candidate is required to have |cos(theta_pi0)| < 0.95 in the pi0 rest frame")
  .note(:photon_timing, "each photon is required to be within +/- 500 ns of the most energetic photon's timing")

# Attach the decay card and render the full selection into the BOSS algorithm
algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Run the same selection on data, inclusive MC and the 500k exclusive signal MC
root_files = algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])