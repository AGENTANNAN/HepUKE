# Dataset preparation
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Signal decay card: J/psi -> gamma X(2370), X(2370) -> K*(892)0 anti-K0 + c.c.
# with K*(892)0 -> K_S0 pi0 (and c.c. via anti-K*(892)0 -> K_S0 pi0), K_S0 -> pi+ pi-, pi0 -> gamma gamma
# Alias X(2370) as pseudoscalar dummy resonance (mass ~ 2.359 GeV) modelled via PHSP.
decay_card_signal = <<~DECAYCARD
  Alias      Xhad    eta(2S)

  Decay J/psi
  1.0000 gamma Xhad                          PHSP;
  Enddecay

  Decay Xhad
  0.5000 K*0 anti-K0                         PHSP;
  0.5000 anti-K*0 K0                         PHSP;
  Enddecay

  Decay K*0
  1.0000 K_S0 pi0                            VSS;
  Enddecay

  Decay anti-K*0
  1.0000 K_S0 pi0                            VSS;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-                             PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                         PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_X2370_KstarKbar_KsKsPi0"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Event selection for J/psi -> gamma K_S0 K_S0 pi0
alg_name = "JpsiGammaKsKsPi0"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz    100.0
                  Vr    10.0
                  nChrp ">=2"     # two pi+ from two K_S0 decays
                  nChrn ">=2"     # two pi- from two K_S0 decays
                  nNet  "==0"
                }
                .select_photon {
                  tdc_emc_start 0
                  tdc_emc_end 14
                  angle_to_track 10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam ">=3"      # 1 radiative gamma + 2 gammas from pi0
                }
                .assign({ :chrgp => :pip, :chrgn => :pim })
                .secondary_vertex_fit([:pip, :pim]) {
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
                .secondary_vertex_fit([:pip, :pim]) {
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
                .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 ">=1"
                }
                .kinematic_fit([:gamma, :K_S0, :K_S0, :pi0]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }

alg.note(:kstar_selection,
         "Downstream selection: require at least one K_S0-pi0 combination with " \
         "|M(K_S0 pi0) - m_K*(892)0| <= 50 MeV/c^2 to select the K*(892)0 -> K_S0 pi0 candidate.")
   .note(:analysis_selection_reference,
         "Full event-selection criteria follow Ref. [28] (arXiv:2605.26495) for " \
         "J/psi -> gamma K_S0 K_S0 pi0; the background contributions are estimated " \
         "to be negligible after selection.")

alg.with_decay_card(decay_card_signal).apply(event_selection)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
