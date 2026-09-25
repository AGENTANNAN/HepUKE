# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# J/psi (3.097 GeV) real data and corresponding inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal process (EvtGen format):
#   J/psi -> gamma eta_c, eta_c -> omega phi,
#   omega -> pi+ pi- pi0, phi -> K+ K-, pi0 -> gamma gamma
# Final state: 3 gamma pi+ pi- K+ K-
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.00000 gamma eta_c            PHSP;
    Enddecay

    Decay eta_c
    1.00000 omega phi              PHSP;
    Enddecay

    Decay omega
    1.00000 pi+ pi- pi0            OMEGA_DALITZ;
    Enddecay

    Decay phi
    1.00000 K+ K-                  VSS;
    Enddecay

    Decay pi0
    1.00000 gamma gamma            PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal channel: 1,000,000 events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gamma_etac_omega_phi"
  config.related_dataset = jpsi_data          # associated real dataset
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "GammaEtaCOmegaPhi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})       # sqrt(s) = 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {
      cos_theta 0.93     # |cos(theta)| < 0.93
      Vz        10.0     # |Vz| < 10 cm along the beam axis
      Vr        1.0      # Vr < 1 cm in the transverse plane
      nChrp     "==2"    # exactly two positive charged tracks
      nChrn     "==2"    # exactly two negative charged tracks
      nNet      "==0"    # net charge zero
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025   # >= 25 MeV in the EMC barrel
      energyThreshold_e 0.050   # >= 50 MeV in the EMC endcap
      nGam              ">=3"   # at least three photons
  }
  .pid(method: :probability) {           # probability PID method
      prob_cut 0.001                     # PID probability > 0.001
      identify :kaon, against: [:pion, :proton]   # K+ and K-
      identify :pion, against: [:kaon, :proton]   # pi+ and pi-
      nkp  "==1"                         # exactly one K+
      nkm  "==1"                         # exactly one K-
      npip "==1"                         # exactly one pi+
      npim "==1"                         # exactly one pi-
  }
  # Reconstruct pi0 -> gamma gamma with a 1C Kalman fit: constrain the two-photon
  # invariant mass to the nominal pi0 mass, chi2 < 25
  .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"                         # at least one pi0 candidate
  }
  # omega mass window: |M(pi+ pi- pi0) - M_omega| < 40 MeV/c^2
  .invariant_mass_of(:pip, :pim, :pi0).within(0.743, 0.823)
  # phi mass window: |M(K+ K-) - M_phi| < 15 MeV/c^2
  .invariant_mass_of(:kp, :km).within(1.004, 1.034)
  # eta' veto: reject 943 < M(pi+ pi- 3gamma) < 969 MeV/c^2
  .invariant_mass_of(:pip, :pim, :gamma, :gamma, :gamma).out_of(0.943, 0.969)
  # 5C kinematic fit: 4C four-momentum conservation + 1C pi0 mass constraint.
  # Uses the three individual photons so that the pi0 mass constraint is applied
  # inside the fit (5 constraints in total).
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma]) {
      nominal                            # nominal fit (only this one is saved)
      vertex_fit([0, 1, 2, 3])           # primary vertex fit on K+ K- pi+ pi-
      constrain_four_momentum            # 4C: total four-momentum -> CMS energy
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 200                       # loose BOSS cut; tight chi2 < 20 in ROOT
  }

# BOSS-side procedure that cannot be expressed in the DSL:
# the pi0 decay-angle cut |cos(theta_decay)| < 0.95
my_algorithm
  .note(:pi0_decay_angle_cut, "pi0 -> gamma gamma candidates are required to " \
        "satisfy |cos(theta_decay)| < 0.95, where theta_decay is the angle of the " \
        "decay photons in the pi0 rest frame; applied in the BOSS selector to " \
        "suppress combinatorial background")
  .with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])