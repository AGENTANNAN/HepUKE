### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # corresponding inclusive MC

# Decay card: J/psi -> gamma eta', eta' -> e+ e- omega, omega -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000000 gamma eta' HELAMP 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay eta'
    1.000000 e+ e- omega PHSP;
    Enddecay

    Decay omega
    1.000000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 5M-event exclusive signal MC
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "jpsi_gamma_etap_ee_omega"
  config.related_dataset = jpsi_data
  config.events         = 5_000_000
  config.decay_card     = decay_card_signal
  config.cross_section  = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiGammaEtapEEOmega"
my_alg = Algorithm.new(alg_name)
my_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
  .select_track {                       # 4 charged tracks: 2 positive, 2 negative, net charge 0
    cos_theta 0.93
    Vr        2.0
    Vz        20.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {                      # >= 3 photons, > 10 deg from any charged track
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"
  }
  .pid(method: :probability) {          # at least 1 pi+, 1 pi-, 1 e+, 1 e-
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]    # pi vs K separation (leptons handled above)
    npip ">=1"
    npim ">=1"
    nlp  ">=1"
    nlm  ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 -> gamma gamma (1C mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0  ">=1"
  }
  .kinematic_fit([:gamma, :lp, :lm, :pip, :pim, :pi0]) {   # nominal 4C fit on gamma e+ e- pi+ pi- pi0
    nominal
    invariant_mass_of(:pip, :pim, :pi0, :lp, :lm).within(0.90, 1.02)   # eta' candidate mass window
    constrain_four_momentum
    chi2_cut 70
  }
  .assign({:lp => :pip, :lm => :pim})   # reinterpret leptons as pions for the competing hypothesis
  .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pip, :pim, :pim]) {   # competing J/psi -> pi+ pi- pi+ pi- gamma gamma gamma;
    constrain_four_momentum                                            # no chi2_cut / no nominal -> chi2 stored for ROOT-level e/pi mis-ID veto
  }

# BOSS-side procedures that cannot be expressed in the current DSL
my_alg.note(:leading_photon_energy,
            "the most energetic photon in the event is required to have E > 1 GeV")
     .note(:pi0_mass_window,
           "require |M(gamma gamma) - m_pi0| < 0.015 GeV/c^2, evaluated for the two lowest-energy photons (the specific photon pairing is not expressible in the DSL)")
     .note(:conversion_veto,
           "photon-conversion veto: retain photons with conversion transverse radius Rxy < 2 cm; otherwise reject if cos(theta_eg) > 0 and -1.0 < dxy < 0.8 cm")
     .note(:background_veto,
           "suppress pi0 -> gamma e+ e- contamination, and veto eta -> gamma e+ e- : reject events with 0.52 < M(gamma e+ e-) < 0.58 GeV/c^2 and M(gamma e+ e-) < 0.5 GeV/c^2")

# Generate the algorithm for the signal process and execute
my_alg.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])