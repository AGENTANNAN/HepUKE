# Core DSL classes are loaded automatically at execution.

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi real data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # Matching inclusive MC sample

# Decay card: J/psi -> gamma eta_c, eta_c -> Xi0 anti-Xi0,
#             Xi0 -> Lambda pi0, Lambda -> p pi-  (and charge conjugate)
# eta_c / Xi0 production modelled in phase space, Lambda via hyperonic weak decay.
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0000 Xi0 anti-Xi0 PHSP;
  Enddecay

  Decay Xi0
  1.0000 Lambda0 pi0 PHSP;
  Enddecay

  Decay anti-Xi0
  1.0000 anti-Lambda0 pi0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 500k-event exclusive MC with the full decay chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etac_Xi0Xibar0"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "Xi0Xibar0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})       # sqrt(s) = 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})
            # BOSS-side procedures that have no dedicated DSL construct:
            .note(:helix_correction,
                  "helix-parameter correction applied to all charged tracks before the 4C kinematic fit")
            .note(:pid_correction_method,
                  "momentum-dependent p/pi assignment: charged tracks with p > 0.247 GeV/c " \
                  "treated as protons and p < 0.247 GeV/c as pions, applied in addition to the probability PID")
            .note(:lambda_mass_window,
                  "Lambda/anti-Lambda secondary-vertex-fit candidates required to satisfy " \
                  "vertex-fit chi2 < 200 and vertex mass within [1.111, 1.120] GeV/c^2")

event_selection = Selection.new
event_selection.select_track {                            # Charged track selection
                  cos_theta 0.93        # |cos(theta)| < 0.93
                  nChrp ">=2"           # at least 2 positive tracks
                  nChrn ">=2"           # at least 2 negative tracks
                  nTot  ">=4"           # at least 4 charged tracks in total
                  nNet  "==0"           # net charge zero
                }
                .select_photon {                          # Photon selection
                  tdc_emc_start 0            # EMC TDC window start
                  tdc_emc_end 14             # EMC TDC window end
                  energyThreshold_b 0.025    # 25 MeV in the barrel
                  energyThreshold_e 0.050    # 50 MeV in the endcap
                  angle_to_track 10.0        # at least 10 deg from any charged track
                  nGam ">=5"                 # at least 5 photons
                }
                .select_isolated_photon {                 # Isolated photon selection
                  angle_to_prp_track 20.0    # at least 20 deg from (anti)proton tracks
                  angle_to_prm_track 20.0
                  nGam ">=5"                 # at least 5 isolated photons
                }
                .pid(method: :probability) {              # PID (probability method)
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]  # protons against kaons and pions
                  identify :pion,   against: [:kaon]         # pions against kaons
                }
                .kalman_kinematic_fit([:gamma, :gamma]) { # Form pi0 from gamma-gamma pairs (1C mass constraint)
                  invariant_mass_of(:gamma, :gamma).within(0.098, 0.165)  # M(gamma gamma) window
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 200                            # chi2 < 200
                  npi0 ">=2"                              # at least two pi0 candidates
                }
                .secondary_vertex_fit([:prp, :pim]) {     # Lambda -> p pi-
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
                .secondary_vertex_fit([:prm, :pip]) {     # anti-Lambda -> anti-p pi+
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
                # Nominal 4C kinematic fit on gamma Lambda anti-Lambda pi0 pi0
                .kinematic_fit([:gamma, :Lambda, :Lambda_bar, :pi0, :pi0]) {
                  nominal
                  constrain_four_momentum                 # 4C energy-momentum constraint
                  invariant_mass_of(:Lambda, :pi0).within(1.296, 1.331)      # Xi0 -> Lambda pi0 window
                  invariant_mass_of(:Lambda_bar, :pi0).within(1.296, 1.331)  # anti-Xi0 window
                  chi2_cut 200                            # loose chi2 cut (tight cut in ROOT)
                }
                # Competing 4C hypothesis WITHOUT the radiative photon; chi2 stored for background suppression
                .kinematic_fit([:Lambda, :Lambda_bar, :pi0, :pi0]) {
                  constrain_four_momentum
                }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])