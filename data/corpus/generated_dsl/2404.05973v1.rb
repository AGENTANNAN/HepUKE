# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
data_4180  = DatasetManager.real_data.find("703_4180")      # 4.178 GeV energy-point real data
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")   # Corresponding inclusive MC sample

# Decay card for the signal process (EvtGen format, EvtGen particle names).
# No intermediate resonance is specified, so psi(4260) is used as the KKMC top mother.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s-                 PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+                 PHSP;
    Enddecay

    Decay D_s+
    1.000 pi+ phi                    PHSP;
    Enddecay

    Decay phi
    1.000 e+ e-                      VLL;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi-                  PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC of the full D_s*+ D_s- decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4180_DsstarDs"
  config.related_dataset = data_4180
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "DsstarDs"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.178]})   # 4.178 GeV center-of-mass energy

event_selection = Selection.new
  .select_track {              # Charged track selection
    cos_theta 0.93             # |cos(theta)| < 0.93
    Vz        10.0             # |Vz| < 10 cm
    Vr        1.0              # Vr < 1 cm in the transverse plane
    nChrp     ">=2"            # At least two positively charged tracks
    nChrn     ">=1"            # At least one negatively charged track
  }
  .select_photon {             # Photon selection
    tdc_emc_start     0        # TDC start
    tdc_emc_end       14       # TDC end
    angle_to_track    10.0     # At least 10 degrees from any charged track
    energyThreshold_b 0.025    # 25 MeV in the barrel
    energyThreshold_e 0.050    # 50 MeV in the endcap
    nGam              ">=1"    # At least one photon
  }
  .pid(method: :probability) { # PID by the probability method
    prob_cut 0.001             # L'(e) > 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6  # p > 1.0 -> lepton; EMC energy > 0.6 -> electron
    identify :pion, against: [:kaon, :proton]   # pi+/pi- separated from K and p
    identify :kaon, against: [:pion, :proton]   # K+/K- separated from pi and p
  }
  # Final (nominal) kinematic fit; loose chi2 cut, tight optimisation done in ROOT
  .kinematic_fit([:gamma, :pip, :ep, :em, :kp, :km, :pim]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:ep, :em).within(0.98, 1.04)         # e+e- (phi) window
    invariant_mass_of(:pip, :ep, :em).within(1.88, 2.02)   # D_s+ candidate mass window
    chi2_cut 200
  }

my_algorithm
  .note(:pid_correction_method, "electron PID additionally requires L'(e)/(L'(e)+L'(pi)+L'(K)) > 0.8 and E/p > 0.8 (0.7) for track momentum above (below) 0.4 GeV/c; only the L'(e) > 0.001 cut and the high-momentum lepton branches are expressible in the DSL")
  .note(:background_veto, "gamma -> e+e- conversions vetoed by discarding e+e- pairs whose vertex lies 2.0-8.0 cm from the IP, with photon recovery in a 5 degree cone around the electron tracks")
  .note(:transition_photon_selection, "when multiple photons are present, the D_s*+ -> gamma D_s+ transition photon is chosen as the one whose D_s+ gamma recoil mass is closest to m(D_s*+)")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute the algorithm on real data, inclusive MC and exclusive signal MC
root_files = my_algorithm.execute_on([data_4180, incMC_4180, exMC_signal])