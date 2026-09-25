# Dataset preparation
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'                       PHSP;
  Enddecay

  Decay eta'
  1.0000 eta pi+ pi-                      ETAPRIME_DALITZ;
  Enddecay

  Decay eta
  1.0000 gamma gamma                      PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etap_etapipi"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Event selection (BOSS): J/psi -> gamma eta', eta' -> eta pi+ pi-, eta -> gamma gamma
alg_name = "JpsiGammaEtapEtaPiPi"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz        10.0
                  Vr        1.0
                  nChrp     "==1"
                  nChrn     "==1"
                  nNet      "==0"
                }
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    20.0     # >= 20 deg isolation from charged tracks
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam              ">=3"    # more than two photons
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :pion, against: [:kaon]
                  npip ">=1"                 # at least one identified pion
                }
               .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200               # keep loose; paper uses < 200
                }

algorithm
  .note(:best_gamma_combination,
        "Best gamma gamma gamma pi+ pi- combination chosen by the minimum " \
        "chi^2 of the 4C kinematic fit over all photon candidates.")
  .note(:eta_selection,
        "eta candidate chosen from the gamma gamma pair with invariant mass " \
        "closest to nominal m_eta; eta signal region 0.518 - 0.578 GeV/c^2 " \
        "applied at ROOT level.")

algorithm.with_decay_card(decay_card_signal).apply(event_selection)
algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
