### Dataset preparation ###
# J/psi(3097) real data and its inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the full J/psi -> Lambda anti-Lambda signal chain
# (all sub-decays: Lambda -> p pi-, anti-Lambda -> anti-p pi+, Sigma+ -> p pi0, pi0 -> gamma gamma)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Lambda0 anti-Lambda0        PHSP;
    Enddecay

    Decay Lambda0
    1.000  p+  pi-                     HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000  anti-p-  pi+                HypWK;
    Enddecay

    Decay Sigma+
    1.000  p+  pi0                     PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                PHSP;
    Enddecay

    End
DECAYCARD

# 1M exclusive signal MC events for the full J/psi -> Lambda anti-Lambda decay card
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_lambdabar_sigma_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "LambdaBarPi0Proton"        # single tag (anti-Lambda) + double tag (Sigma+ -> p pi0)
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # J/psi CMS energy (GeV)

event_selection = Selection.new
event_selection.select_track {         # Charged track selection
                  cos_theta 0.93       # |cos(theta)| < 0.93
                  Vz        10.0       # |Vz| < 10 cm
                  Vr        1.0        # Vr < 1 cm
                  nChrp     ">=1"      # at least one positive track
                  nChrn     ">=1"      # at least one negative track
                }
               .select_photon {        # Photon selection
                  tdc_emc_start     0      # EMC TDC start
                  tdc_emc_end       14     # EMC TDC end
                  energyThreshold_b 0.025  # 25 MeV (barrel)
                  energyThreshold_e 0.050  # 50 MeV (endcap)
                  angle_to_track    10.0   # >= 10 deg separation from any charged track
                  nGam              ">=2"  # at least two photons
                }
               .pid(method: :probability) {   # PID by probability method
                  prob_cut 0.001               # probability > 0.001
                  identify :proton, against: [:pion, :kaon]   # p+ and anti-p- vs pi, K
                }
               .remove([:prp <= :chrgp, :prm <= :chrgn])   # keep non-proton tracks for pions
               .assign({:chrgp => :pip, :chrgn => :pim})   # remaining +/- tracks treated as pions
               .secondary_vertex_fit([:prm, :pip]) {       # Reconstruct anti-Lambda -> anti-p pi+
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
               }
               .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C Kalman fit: pi0 -> gamma gamma
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 ">=1"                                # at least one pi0 candidate
               }
               .kinematic_fit([:Lambda_bar, :pi0, :prp]) {  # 4C kinematic fit on anti-Lambda pi0 p
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
               }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])