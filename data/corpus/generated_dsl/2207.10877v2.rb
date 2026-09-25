# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# psi(3686) (3.686 GeV) real data and its inclusive MC  (709_3686, 448.1x10^6 events)
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for psi(3686) -> Lambda_c+ anti-Sigma-  (+ c.c.)
#   Lambda_c+   -> p K- pi+
#   anti-Sigma- -> anti-p- pi0 , pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Lambda_c+ anti-Sigma- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 p+ K- pi+ PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000 anti-p- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive signal MC events, generated with the same decay card
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_LcpSigmam"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "LcSigma"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})          # sqrt(s) = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})
            # BOSS-side procedure that cannot be expressed in the DSL:
            .note(:helix_correction, "For Monte Carlo the charged tracks are corrected with the helix-parameter
              correction before the 5C kinematic fit; all downstream track quantities use the corrected helix
              parameters (the efficiency difference is estimated by re-running the BOSS selection).")
            .note(:loose_vertex_cuts, "The two daughters of the long-lived anti-Sigma- (anti-p- and pi0) are
              selected with looser Vr/Vz vertex cuts than the global charged-track selection in the generated
              code, so that the displaced anti-Sigma- decay is retained.")

event_selection = Selection.new
event_selection.select_track {          # charged-track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz 10.0                             # |Vz| < 10 cm
    Vr 1.0                              # Vr < 1 cm
    nChrp ">=2 && nChrp <=10"           # 2-10 positive tracks
    nChrn ">=2 && nChrn <=10"           # 2-10 negative tracks
    nTot  ">=4 && nTot <=10"            # 4-10 charged tracks in total
    }
  .select_photon {                      # photon selection
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0                 # more than 10 degrees from the nearest charged track
    energyThreshold_b 0.025             # E > 25 MeV in the barrel  (|cos(theta)| < 0.80)
    energyThreshold_e 0.050             # E > 50 MeV in the endcap (0.86 < |cos(theta)| < 0.92)
    nGam ">=2"                          # at least two photons for pi0 -> gamma gamma
    }
  .pid(method: :probability) {          # PID by highest-probability assignment (probability method)
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]     # p+ and anti-p-
    identify :kaon,   against: [:pion, :proton]   # K+ and K-
    identify :pion,   against: [:kaon, :proton]   # pi+ and pi-
    nprp "==1"                           # exactly one p
    nprm "==1"                           # exactly one anti-p-
    nkm  "==1"                           # exactly one K-
    npip "==1"                           # exactly one pi+
    }
  .kalman_kinematic_fit([:gamma, :gamma]) {       # 1C mass-constrained fit: pi0 -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"                                    # at least one pi0 candidate
    }
  .kinematic_fit([:prp, :km, :pip, :prm, :pi0]) { # 5C fit: 4-momentum conservation + pi0 mass constraint
    invariant_mass_of(:pip, :prm).out_of(1.090, 1.130)   # veto M(pi+ anti-p-) : anti-Lambda region
    invariant_mass_of(:km, :pip).out_of(0.756, 1.036)    # veto M(K- pi+) : K*/phi region
    invariant_mass_of(:prm, :pi0).between(1.150, 1.230)  # Sigma- mass window
    constrain_four_momentum                     # the pi0 mass constraint is carried by the mass-constrained pi0
    chi2_cut 60
    nominal                                     # nominal fit; the best candidate is the smallest 5C chi2 (default)
    }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])