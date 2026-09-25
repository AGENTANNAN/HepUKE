# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi real data at 3.097 GeV (1310.6e6 events)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097") # Corresponding inclusive MC sample

# ---- Decay cards (EvtGen): J/psi -> gamma eta', one per exclusive eta' decay mode ----

# Mode 1: eta' -> gamma pi+ pi-
decay_card_gammapipi = <<~DECAYCARD
  Decay J/psi
  1.0 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.0 gamma pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode 2: eta' -> eta pi+ pi-  (eta -> gamma gamma)
decay_card_etapipi = <<~DECAYCARD
  Decay J/psi
  1.0 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.0 eta pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode 3: eta' -> eta pi0 pi0  (eta, pi0 -> gamma gamma)
decay_card_etapi0pi0 = <<~DECAYCARD
  Decay J/psi
  1.0 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.0 eta pi0 pi0 PHSP;
  Enddecay

  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode 4: eta' -> gamma omega  (omega -> pi+ pi- pi0, pi0 -> gamma gamma)
decay_card_gammaomega = <<~DECAYCARD
  Decay J/psi
  1.0 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.0 gamma omega PHSP;
  Enddecay

  Decay omega
  1.0 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode 5: eta' -> gamma gamma
decay_card_gammagamma = <<~DECAYCARD
  Decay J/psi
  1.0 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- Exclusive MC samples (100k events each) ----
exMC_gammapipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaetaprime_gammapipi"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_gammapipi
  config.cross_section   = :default
end

exMC_etapipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaetaprime_etapipi"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_etapipi
  config.cross_section   = :default
end

exMC_etapi0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaetaprime_etapi0pi0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_etapi0pi0
  config.cross_section   = :default
end

exMC_gammaomega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaetaprime_gammaomega"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_gammaomega
  config.cross_section   = :default
end

exMC_gammagamma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaetaprime_gammagamma"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_gammagamma
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# =====================================================================
# Mode 1: J/psi -> gamma eta', eta' -> gamma pi+ pi-   (final: gamma gamma pi+ pi-)
# =====================================================================
alg_gammapipi = Algorithm.new("JpsiEtaPrimeGammaPiPi")
alg_gammapipi.set_header(["JpsiEtaPrimeGammaPiPiAlg/JpsiEtaPrimeGammaPiPi.h"])
             .set_constant({"ECMS" => [:double, 3.097]})
             .set_alias({"std::vector<double>" => "Vdouble"})
             # Inexpressible BOSS-side procedure: inclusive eta' tag via photon conversion
             .note(:inclusive_etaprime_conversion_tag,
                   "Inclusive eta' tag uses the photon-conversion method: the radiative photon from " \
                   "J/psi -> gamma eta' converts to e+e- with a pi0 conversion veto; the highest-energy " \
                   "photon (E ~ 1.4 GeV) is taken as the radiative candidate. Normalized with " \
                   "BF(J/psi -> gamma eta') = (5.27 +/- 0.03 +/- 0.05)e-3 and an efficiency correction of " \
                   "1.0085 +/- 0.0050.")

sel_gammapipi = Selection.new
  .select_track {                       # Charged track selection
      cos_theta 0.93
      Vz        100.0
      Vr        10.0
      nChrp     "==1"                   # exactly one pi+
      nChrn     "==1"                   # exactly one pi-
      nNet      "==0"                   # net charge zero
  }
  .select_photon {                      # Photon selection
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=2"           # at least 2 photons
  }
  .pid(method: :probability) {          # PID: pi/K separation
      prob_cut 0.001
      identify :pion, against: [:kaon]
      npip "==1"
      npim "==1"
  }
  .kinematic_fit([:gamma, :gamma, :pip, :pim]) {   # 4C fit: gamma gamma pi+ pi-
      nominal
      constrain_four_momentum
      chi2_cut 200                      # loose BOSS cut; tight cut in ROOT
  }

alg_gammapipi.with_decay_card(decay_card_gammapipi).apply(sel_gammapipi)

# =====================================================================
# Mode 2: J/psi -> gamma eta', eta' -> eta pi+ pi-, eta -> gamma gamma
#         (final: gamma eta(->gamma gamma) pi+ pi-)
# =====================================================================
alg_etapipi = Algorithm.new("JpsiEtaPrimeEtaPiPi")
alg_etapipi.set_header(["JpsiEtaPrimeEtaPiPiAlg/JpsiEtaPrimeEtaPiPi.h"])
           .set_constant({"ECMS" => [:double, 3.097]})
           .set_alias({"std::vector<double>" => "Vdouble"})

sel_etapipi = Selection.new
  .select_track {
      cos_theta 0.93
      Vz        100.0
      Vr        10.0
      nChrp     "==1"
      nChrn     "==1"
      nNet      "==0"
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=3"           # at least 3 photons
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon]
      npip "==1"
      npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {        # reconstruct eta from gamma gamma (1-C)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 25
      neta ">=1"                        # at least 1 eta candidate
  }
  .kinematic_fit([:gamma, :eta, :pip, :pim]) {     # 4C fit: gamma eta pi+ pi-
      nominal
      constrain_four_momentum
      chi2_cut 200
  }

alg_etapipi.with_decay_card(decay_card_etapipi).apply(sel_etapipi)

# =====================================================================
# Mode 3: J/psi -> gamma eta', eta' -> eta pi0 pi0, eta/pi0 -> gamma gamma
#         (final: gamma eta(->gamma gamma) pi0 pi0)
# =====================================================================
alg_etapi0pi0 = Algorithm.new("JpsiEtaPrimeEtaPi0Pi0")
alg_etapi0pi0.set_header(["JpsiEtaPrimeEtaPi0Pi0Alg/JpsiEtaPrimeEtaPi0Pi0.h"])
             .set_constant({"ECMS" => [:double, 3.097]})
             .set_alias({"std::vector<double>" => "Vdouble"})

sel_etapi0pi0 = Selection.new
  .select_track {
      cos_theta 0.93
      Vz        100.0
      Vr        10.0
      nChrp     "==0"                   # neutral final state: no charged tracks
      nChrn     "==0"
      nNet      "==0"
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=5"           # at least 5 photons
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {        # reconstruct eta from gamma gamma (1-C)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 25
      neta ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {        # reconstruct pi0 from gamma gamma (1-C)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=2"                        # at least 2 pi0 candidates
  }
  .kinematic_fit([:gamma, :eta, :pi0, :pi0]) {     # 4C fit: gamma eta pi0 pi0
      nominal
      constrain_four_momentum
      chi2_cut 200
  }

alg_etapi0pi0.with_decay_card(decay_card_etapi0pi0).apply(sel_etapi0pi0)

# =====================================================================
# Mode 4: J/psi -> gamma eta', eta' -> gamma omega,
#         omega -> pi+ pi- pi0, pi0 -> gamma gamma
#         (final: gamma gamma pi+ pi- pi0)
# =====================================================================
alg_gammaomega = Algorithm.new("JpsiEtaPrimeGammaOmega")
alg_gammaomega.set_header(["JpsiEtaPrimeGammaOmegaAlg/JpsiEtaPrimeGammaOmega.h"])
              .set_constant({"ECMS" => [:double, 3.097]})
              .set_alias({"std::vector<double>" => "Vdouble"})

sel_gammaomega = Selection.new
  .select_track {
      cos_theta 0.93
      Vz        100.0
      Vr        10.0
      nChrp     "==1"
      nChrn     "==1"
      nNet      "==0"
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=3"           # at least 3 photons
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon]
      npip "==1"
      npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {        # reconstruct pi0 from gamma gamma (1-C)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"
  }
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :pi0]) {   # 4C fit: gamma gamma pi+ pi- pi0
      nominal
      constrain_four_momentum
      chi2_cut 200
  }

alg_gammaomega.with_decay_card(decay_card_gammaomega).apply(sel_gammaomega)

# =====================================================================
# Mode 5: J/psi -> gamma eta', eta' -> gamma gamma   (final: gamma gamma gamma)
# =====================================================================
alg_gammagamma = Algorithm.new("JpsiEtaPrimeGammaGamma")
alg_gammagamma.set_header(["JpsiEtaPrimeGammaGammaAlg/JpsiEtaPrimeGammaGamma.h"])
              .set_constant({"ECMS" => [:double, 3.097]})
              .set_alias({"std::vector<double>" => "Vdouble"})

sel_gammagamma = Selection.new
  .select_track {
      cos_theta 0.93
      Vz        100.0
      Vr        10.0
      nChrp     "==0"                   # neutral final state: no charged tracks
      nChrn     "==0"
      nNet      "==0"
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=3"           # at least 3 photons
  }
  .kinematic_fit([:gamma, :gamma, :gamma]) {       # 4C fit: gamma gamma gamma
      nominal
      constrain_four_momentum
      chi2_cut 200
  }

alg_gammagamma.with_decay_card(decay_card_gammagamma).apply(sel_gammagamma)

### Execute on datasets ###
root_files_gammapipi  = alg_gammapipi.execute_on([jpsi_data, jpsi_incMC, exMC_gammapipi])
root_files_etapipi    = alg_etapipi.execute_on([jpsi_data, jpsi_incMC, exMC_etapipi])
root_files_etapi0pi0  = alg_etapi0pi0.execute_on([jpsi_data, jpsi_incMC, exMC_etapi0pi0])
root_files_gammaomega = alg_gammaomega.execute_on([jpsi_data, jpsi_incMC, exMC_gammaomega])
root_files_gammagamma = alg_gammagamma.execute_on([jpsi_data, jpsi_incMC, exMC_gammagamma])