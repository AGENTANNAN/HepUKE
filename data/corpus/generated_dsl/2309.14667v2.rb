# =============================================================================
# BOSS-side DSL for J/psi -> Xi- anti-Xi+  (Xi- -> Lambda pi-, Lambda -> p pi-;
# anti-Xi+ -> anti-Lambda pi+, anti-Lambda -> anti-n pi0, pi0 -> gamma gamma)
# =============================================================================

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi real data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # Corresponding inclusive MC

# Decay card for the full signal chain (EvtGen syntax)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Xi- anti-Xi+      PHSP;
    Enddecay

    Decay Xi-
    1.0000 Lambda0 pi-       PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 anti-Lambda0 pi+  PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-            HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-n- pi0       PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma       PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 2M events for the full decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_XimXipbar"
  config.related_dataset = jpsi_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "XiXibar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta 0.93        # |cos(theta)| < 0.93
                  Vz        10.0        # |Vz| < 10 cm
                  Vr        1.0         # Vr < 1 cm in the transverse plane
                  nChrp     ">=1"       # at least one positive track
                  nChrn     ">=3"       # at least three negative tracks
                }
               .select_photon {         # Photon selection
                  tdc_emc_start     0     # TDC window 0-14
                  tdc_emc_end       14
                  energyThreshold_b 0.025 # 25 MeV barrel threshold
                  energyThreshold_e 0.050 # 50 MeV endcap threshold
                  angle_to_track    20.0  # opening angle to any charged track > 20 deg
                  nGam              ">=2" # at least two photons
                }
               .pid(method: :probability) {   # Probability PID
                  prob_cut 0.001
                  identify :proton, against: [:pion, :kaon]  # p+/pbar vs pi/K
                }
               .remove([:prp <= :chrgp, :prm <= :chrgn])  # drop PID'd (anti)protons from generic charged lists
               .remove(:prp) { condition "three_momentum_of(:prp) < 0.32" } # only p>0.32 GeV/c taken as proton
               .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks treated as pi+/pi-
               .remove(:pip) { condition "three_momentum_of(:pip) > 0.30" } # pi+ must have p<0.30 GeV/c
               .remove(:pim) { condition "three_momentum_of(:pim) > 0.30" } # pi- must have p<0.30 GeV/c
               .secondary_vertex_fit([:prp, :pim]) {   # Reconstruct Lambda -> p pi-
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {   # Reconstruct pi0 -> gamma gamma
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0     ">=1"
                }
               # Nominal kinematic fit: Lambda pi- pi+ pi0 with anti-neutron missing
               .kinematic_fit([:Lambda, :pim, :pip, :pi0]) {
                  nominal
                  miss_track_of :n_bar        # anti-neutron treated as a missing particle
                  constrain_four_momentum     # 4C energy-momentum constraint
                  chi2_cut 200                # loose BOSS-level chi2; tight cut applied in ROOT
                }

# Inexpressible BOSS-side procedures captured for the downstream uncertainty step
my_algorithm
  .note(:mass_windows, "Pre-fit mass windows applied to the reconstructed candidates: " \
        "|M(p pi-) - m_Lambda| < 11 MeV, |M(p pi- pi-) - m_Xi| < 11 MeV; the Lambda and " \
        "Xi- candidates are additionally required to have positive decay lengths and " \
        "|cos(theta_Xi)| < 0.84 in the CM frame.")
  .note(:missing_particle_mass, "In the nominal 4C fit the anti-neutron is treated as a " \
        "missing particle of unknown (unconstrained) mass; the fitted anti-neutron " \
        "four-momentum is used to form M(anti-n gamma gamma pi+).")
  .note(:bdt_photon_classifier, "A BDT-based photon classifier is applied to the selected " \
        "photons, tuned for ~90% signal efficiency and ~55% background rejection.")
  .note(:background_veto, "Photons with an opening angle < 15 degrees relative to the " \
        "anti-Xi+ direction are vetoed to suppress EMC showers originating from the anti-neutron.")
  .note(:blinding, "The angular analysis is performed blinded: the signal region in data is " \
        "hidden until the analysis procedure is frozen.")
  .with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])