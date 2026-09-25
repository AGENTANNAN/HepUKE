### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi (3.097 GeV) real data, round02 (year 2009)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC sample

# Decay card for the signal process J/psi -> gamma A0, A0 -> mu+ mu- (EvtGen format)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma A0          PHSP;
    Enddecay

    Decay A0
    1.0000 mu+ mu-           PHSP;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive MC sample for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gammaA0_mumu"
  config.related_dataset = jpsi_data           # associated real dataset
  config.events          = 200000              # 200k events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "GammaA0ToMuMu"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})        # ECMS = 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})
            # BOSS-side procedures that cannot be expressed in the formal DSL
            .note(:mass_dependent_photon_isolation,
                  "photon isolation from the nearest extrapolated charged track is mass dependent: the
                   nominal requirement is > 10 degrees (applied in select_photon), and is tightened to
                   > 20 degrees for candidates with m(A0) <= 0.3 GeV/c2, evaluated from the reconstructed
                   A0 mass; not expressible as a fixed photon cut in the DSL")
            .note(:electron_veto,
                  "both identified lepton tracks are required to satisfy E_cal/p < 0.9 (electron veto) to
                   suppress the J/psi -> gamma A0, A0 -> e+e- channel; applied on top of the
                   high-momentum lepton identification")
            .note(:muon_identification_criteria,
                  "at least one of the two tracks must be identified as a muon under the muon mass
                   hypothesis: 0.1 < E_cal < 0.3 GeV, |delta_t_TOF| < 0.26 ns, and MUC penetration depth
                   > (-40 + 70*p/GeV) cm for 0.5 <= p <= 1.1 GeV/c, or > 40 cm for p > 1.1 GeV/c")

# Build the event selection chain
event_selection = Selection.new
event_selection.select_track {                     # Charged track selection
                  cos_theta   0.93                 # |cos(theta)| < 0.93
                  Vz          10.0                 # |Vz| < 10 cm
                  Vr          1.0                  # Vr < 1 cm in the transverse plane
                  nChrp       "==1"                # exactly one positively charged track
                  nChrn       "==1"                # exactly one negatively charged track
                  nNet        "==0"                # net charge zero
                }
               .select_photon {                    # Photon selection
                  tdc_emc_start     0              # EMC TDC start
                  tdc_emc_end       14             # EMC TDC end
                  angle_to_track    10.0           # angle to the nearest extrapolated track > 10 degrees
                  energyThreshold_b 0.025          # E > 25 MeV in the barrel
                  energyThreshold_e 0.050          # E > 50 MeV in the endcap
                  nGam              ">=1"          # at least one photon
                }
               .pid(method: :probability) {        # PID: probability method
                  prob_cut   0.001                 # PID probability > 0.001
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6  # p>1.0 GeV/c -> lepton; EMC energy>0.6 GeV -> electron, else muon
                  nlp        ">=1"                 # at least one positive lepton
                  nlm        ">=1"                 # at least one negative lepton
                }
               .secondary_vertex_fit([:lp, :lm]) { # Fit the two muon candidates to a common vertex to form the A0
                  build_virtual_particle(:A0).by_minimizing_verfit_chi2   # A0 mass is not fixed; use the vertex chi2
                  remove_used_particle_from_candidate_list                # prevent track reuse
                }
               .kinematic_fit([:gamma, :A0]) {     # Nominal 4C kinematic fit on gamma mu+ mu-
                  nominal                          # nominal fit (corrected four-momenta used downstream)
                  constrain_four_momentum          # 4C energy-momentum constraint
                  chi2_cut 200                     # loose chi2 cut; the minimum-chi2 candidate is chosen automatically
                }

# Generate the complete algorithm for the process defined in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC and signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])