# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
psipp_data  = DatasetManager.real_data.find("712_3773")     # psi(3770) real data (BOSS 7.1.2)
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

# Decay card for the signal process:
#   e+e- -> e+e- pi0 (two-photon fusion), generated with the ConExc model.
# The literal token `ConExc` makes the DSL switch to the no-KKMC simulation template,
# so no `Particle vpho` line is written here (the DSL injects it per energy point).
# `-1 vhdr 11 -11 111` selects the user-defined final state by PDG codes
# (e+, e-, pi0) for the two-photon channel.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.0 ConExc -1 vhdr 11 -11 111;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 5M-event exclusive MC of the single-tag signal
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psipp_ee_pi0_conexc"
  config.related_dataset = psipp_data            # tie the signal MC to the 3.773 GeV data
  config.events          = 5_000_000             # 5M events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default              # ConExc handles the ISR / sigma0(m) modelling
end

### Event selection (BOSS) ###
alg_name = "Pi0TFF"                                # pi0 transition form factor (single-tag two-photon)
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.773]})   # sqrt(s) = 3.773 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

# Full selection chain: tracks -> photons -> PID -> 4C kinematic fit
event_selection = Selection.new
event_selection.select_track {                     # charged track selection
                  cos_theta  0.93                  # |cos(theta)| < 0.93
                  Vz         10.0                  # |Vz| < 10 cm
                  Vr         1.0                   # Vxy < 1 cm
                  nChrp      ">=1"                 # at least one positively charged track
                  nChrn      ">=1"                 # at least one negatively charged track
                }
               .select_photon {                    # photon selection
                  tdc_emc_start     0              # EMC TDC start
                  tdc_emc_end       14             # EMC TDC end
                  energyThreshold_b 0.025          # 25 MeV in the barrel
                  energyThreshold_e 0.050          # 50 MeV in the endcap
                  angle_to_track    10.0           # >= 10 deg from the nearest charged track
                  nGam              ">=2"          # at least two photons (pi0 -> gamma gamma)
                }
               .pid(method: :probability) {        # PID: probability method
                  prob_cut   0.001                 # PID probability > 0.001
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6  # e+/e- (vs pi/K)
                  nlp   ">=1"                      # at least one e+
                  nlm   ">=1"                      # at least one e-
                }
               .kinematic_fit([:lp, :lm, :gamma, :gamma]) {   # nominal 4C fit to e+e- gamma gamma
                  nominal                          # nominal fit -> corrected four-momenta are saved
                  constrain_four_momentum          # 4C energy-momentum conservation
                  chi2_cut 200                     # loose chi^2 cut (tight cut applied in ROOT)
                }

# Inexpressible BOSS-side procedures / criteria (kept for the systematics stage)
my_algorithm
  .note(:pid_correction_method, "electron identification additionally requires E/p > 0.8 on the
    identified lepton tracks; this requirement is enforced on top of the PID probability cut
    (prob_cut 0.001, electrons separated against pions and kaons)")
  .note(:undetected_electron_angle, "the un-tagged (undetected) electron/positron is required to
    satisfy cos(theta_undetected) > 0.99, i.e. it escapes down the beam pipe as expected for the
    single-tag two-photon topology")
  .note(:helicity_angle_window, "background suppression requires |cos(theta_H)| < 0.8 for the pi0
    helicity angle evaluated in the two-photon centre-of-mass frame")
  .note(:photon_isolation_ratio, "background suppression requires the photon isolation ratio
    R_gamma <= 0.15, computed from the EMC shower energy inside/outside a cone around the photon")
  .note(:extra_energy_sum, "background suppression requires the total extra (unassociated) shower
    energy in the event to satisfy sum(E_extra) <= 0.17 GeV")
  .note(:pi0_candidate_selection, "when several gamma-gamma combinations satisfy the selection, the
    pi0 candidate is chosen as the one with minimum p_t* evaluated in the gamma-gamma rest frame")

# Generate the complete algorithm for the process defined in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the ConExc exclusive signal MC
root_files = my_algorithm.execute_on([psipp_data, psipp_incMC, exMC_signal])