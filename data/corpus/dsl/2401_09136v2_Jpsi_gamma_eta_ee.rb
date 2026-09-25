# Paper: arXiv:2401.09136v2
# J/psi -> gamma eta/eta', eta/eta' -> gamma e+ e-
# 10 billion J/psi events at sqrt(s)=3.097 GeV
# Ordinary analysis with 4C kinematic fit

### Dataset preparation ###
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for J/psi -> gamma eta, eta -> gamma e+ e-
decay_card_eta = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta PHSP;
    Enddecay
    Decay eta
    1.000 gamma e+ e- PHSP;
    Enddecay
    End
DECAYCARD

# Decay card for J/psi -> gamma eta', eta' -> gamma e+ e-
decay_card_etap = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta' PHSP;
    Enddecay
    Decay eta'
    1.000 gamma e+ e- PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC samples
exMC_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_jpsi_gamma_eta_ee"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_eta
  config.cross_section = :default
end

exMC_etap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_jpsi_gamma_etap_ee"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_etap
  config.cross_section = :default
end

### Event selection ###

# Two separate Algorithms per Rule T1: different PID cuts for eta vs eta'
# They share the same selection chain (same final state gamma gamma e+ e-)

# -- eta mode --
alg_eta = Algorithm.new("JpsiGammaEtaEE")
alg_eta.set_header(["JpsiGammaEtaEEAlg/JpsiGammaEtaEE.h"])
        .set_constant({"ECMS" => [:double, 3.097]})

# -- eta' mode --
alg_etap = Algorithm.new("JpsiGammaEtapEE")
alg_etap.set_header(["JpsiGammaEtapEEAlg/JpsiGammaEtapEE.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

# Common selection chain
event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==1"
    nChrn "==1"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=2"
  end
  .pid do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  end
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) do
    constrain_four_momentum
    chi2_cut 100
    nominal
  end

# Notes for both modes (inexpressible BOSS-side procedures)
# Custom generator note
alg_eta.note(:custom_generator, "Custom generator incorporating theoretical amplitudes for exclusive eta/eta' -> gamma e+ e- decays; PHSP used as placeholder in decay card")
alg_etap.note(:custom_generator, "Custom generator incorporating theoretical amplitudes for exclusive eta/eta' -> gamma e+ e- decays; PHSP used as placeholder in decay card")

# Electron PID likelihood ratio cuts (inexpressible in DSL PID block)
alg_eta.note(:electron_pid_cut, "Electron PID: likelihood ratio L(e)/(L(e)+L(pi)) > 0.5 for eta mode (applied as post-PID cut in analysis code)")
alg_etap.note(:electron_pid_cut, "Electron PID: likelihood ratio L(e)/(L(e)+L(pi)) > 0.95 for eta' mode (applied as post-PID cut in analysis code)")

# Photon conversion veto
alg_eta.note(:photon_conversion_veto, "Photon conversion veto: reject events with R_xy > 2 cm when cos(theta_eg) > 0 and |Delta_xy| < 0.8 cm; R_xy is conversion length, Delta_xy is distance between track circle intersections")
alg_etap.note(:photon_conversion_veto, "Photon conversion veto: reject events with R_xy > 2 cm when cos(theta_eg) > 0 and |Delta_xy| < 0.8 cm")

# QED background suppression
alg_eta.note(:qed_background_suppression, "Low-energy photon E > 0.15 GeV; angle between low-E photon and e+/e- > 10 deg; angle between radiative photon and e+/e- > 20 deg")
alg_etap.note(:qed_background_suppression, "Low-energy photon E > 0.15 GeV; angle between low-E photon and e+/e- > 10 deg; angle between radiative photon and e+/e- > 20 deg")

# eta/eta' mass window veto on gamma gamma
alg_eta.note(:eta_etap_veto, "|M(gamma gamma) - 0.547| > 0.03 GeV/c^2 and |M(gamma gamma) - 0.958| > 0.03 GeV/c^2 to suppress J/psi -> e+e- eta/eta' with eta/eta' -> gamma gamma")
alg_etap.note(:eta_etap_veto, "|M(gamma gamma) - 0.547| > 0.03 GeV/c^2 and |M(gamma gamma) - 0.958| > 0.03 GeV/c^2 to suppress J/psi -> e+e- eta/eta' with eta/eta' -> gamma gamma")

# Radiative photon identification
alg_eta.note(:radiative_photon_id, "Photon with maximum energy identified as radiative photon from J/psi decay; low-energy photon from eta/eta' Dalitz decay")
alg_etap.note(:radiative_photon_id, "Photon with maximum energy identified as radiative photon from J/psi decay; low-energy photon from eta/eta' Dalitz decay")

# Signal regions and sidebands (post-fit, applied in ROOT analysis)
alg_eta.note(:signal_region, "eta signal region: M(gamma e+ e-) in [0.52, 0.56] GeV/c^2; sidebands [0.50,0.51] U [0.59,0.60] GeV/c^2")
alg_etap.note(:signal_region, "eta' signal region: M(gamma e+ e-) in [0.94, 0.98] GeV/c^2; sidebands [0.88,0.90] U [1.00,1.02] GeV/c^2")

# Primary vertex fit
alg_eta.note(:primary_vertex_fit, "Candidate events required to pass primary vertex fit before selection")
alg_etap.note(:primary_vertex_fit, "Candidate events required to pass primary vertex fit before selection")

alg_eta.with_decay_card(decay_card_eta).apply(event_selection)
alg_etap.with_decay_card(decay_card_etap).apply(event_selection)

alg_eta.execute_on([jpsi_data, jpsi_incMC, exMC_eta])
alg_etap.execute_on([jpsi_data, jpsi_incMC, exMC_etap])