# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(3686) real data at sqrt(s) = 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

# Decay card for the signal process (EvtGen format, EvtGen particle names):
#   psi(3686) -> phi eta eta',  phi -> K+ K-,  eta -> gamma gamma,  eta' -> gamma pi+ pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000  phi  eta  eta'    PHSP;
    Enddecay

    Decay phi
    1.000  K+  K-    VSS;
    Enddecay

    Decay eta
    1.000  gamma  gamma    PHSP;
    Enddecay

    Decay eta'
    1.000  gamma  pi+  pi-    PHSP;
    Enddecay

    End
DECAYCARD

# Create the exclusive signal MC sample: 1,000,000 events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_phi_eta_etap"
  config.related_dataset = psip_data                # associated real dataset (better simulation)
  config.events          = 1000000                  # 1,000,000 signal MC events
  config.decay_card      = decay_card_signal        # decay card defining the signal chain
  config.cross_section   = :default                 # default cross-section
end

### Event selection (BOSS) ###
alg_name     = "PhiEtaEtap"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})      # CMS energy = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"}) # alias for long type names

# Full event-selection chain (identical for data, inclusive MC and exclusive signal MC)
event_selection = Selection.new
  .select_track {                       # Charged track selection
      cos_theta 0.93                    # |cos(theta)| < 0.93
      Vz        10.0                    # |Vz| < 10 cm
      Vr        1.0                     # Vr < 1 cm in the transverse plane
      nChrp     "==2"                   # exactly two positive tracks
      nChrn     "==2"                   # exactly two negative tracks
      nNet      "==0"                   # net charge zero
  }
  .select_photon {                      # Photon selection
      tdc_emc_start     0               # EMC TDC start time
      tdc_emc_end       14              # EMC TDC end time
      energyThreshold_b 0.025           # barrel energy > 25 MeV
      energyThreshold_e 0.050           # endcap energy > 50 MeV
      angle_to_track    10.0            # at least 10 degrees from any charged track
      nGam              ">=3"           # at least three photons
  }
  .pid(method: :probability) {          # PID with the probability method
      prob_cut 0.001                    # probability > 0.001
      identify :kaon, against: [:pion, :proton]   # K+/K- separated from pi and p
      identify :pion, against: [:kaon, :proton]   # pi+/pi- separated from K and p
      nkp  "==1"                        # exactly one K+
      nkm  "==1"                        # exactly one K-
      npip "==1"                        # exactly one pi+
      npim "==1"                        # exactly one pi-
  }
  # Pre-selection 4C fit to K+ K- pi+ pi- gamma gamma gamma
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma]) {
      constrain_four_momentum           # 4-momentum conservation
      chi2_cut 50                       # chi2 < 50
  }
  # Competing 4-photon hypothesis (K+ K- pi+ pi- gamma gamma gamma gamma):
  # no chi2_cut and no nominal -> its chi2 is stored and the veto is applied in ROOT
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma, :gamma]) {
      constrain_four_momentum
  }
  # Competing 2-photon hypothesis (K+ K- pi+ pi- gamma gamma):
  # no chi2_cut and no nominal -> its chi2 is stored and the veto is applied in ROOT
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {
      constrain_four_momentum
  }
  # Nominal 5C fit: 4-momentum conservation + eta -> gamma gamma mass constraint
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma]) {
      nominal                                       # nominal fit (corrected four-momenta are kept)
      constrain_four_momentum                       # 4C energy-momentum constraint
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)  # 5th constraint
      chi2_cut 200                                  # loose chi2 cut; tight cut optimised in ROOT
  }

# BOSS-side procedures that have no formal DSL construct
my_algorithm
  .note(:helix_correction, "helix-parameter corrections are applied to the charged
    tracks of simulated events (exclusive signal MC and inclusive MC) before the
    pre-selection 4C kinematic fit; the resulting efficiency difference is evaluated
    by re-running the BOSS selection with and without the correction")
  .note(:background_veto, "three mass-window vetoes applied after the nominal 5C fit:
    |M(K+K- pi+pi- gamma gamma) - M(chi_cJ)| <= 0.010 GeV (psi(3686) -> gamma chi_cJ),
    |M(gamma gamma gamma K+K-) - M(J/psi)| <= 0.030 GeV and
    |M(gamma pi+pi- K+K-) - M(J/psi)| <= 0.030 GeV (J/psi and eta J/psi backgrounds);
    the masses are computed from the nominal-fit four-momenta, so the windows are
    applied at the ROOT level")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute the same selection chain on data, inclusive MC and exclusive signal MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])