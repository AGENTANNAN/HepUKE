# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
# Real data samples and their matching inclusive MC at the two resonances
psi2S_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) @ 3.686 GeV
psi2S_incMC = DatasetManager.inclusive_mc.find("709_3686")  # matching inclusive MC
jpsi_data   = DatasetManager.real_data.find("708_3097")     # J/psi @ 3.097 GeV
jpsi_incMC  = DatasetManager.inclusive_mc.find("708_3097")  # matching inclusive MC

### Decay cards (EvtGen format) ###
# Mode I : eta' -> gamma pi+ pi-
decay_card_psi2S_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.0000 p+ anti-p- eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_jpsi_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000 p+ anti-p- eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode II : eta' -> eta pi+ pi-, eta -> gamma gamma
decay_card_psi2S_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.0000 p+ anti-p- eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_jpsi_modeII = <<~DECAYCARD
    Decay J/psi
    1.0000 p+ anti-p- eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive signal MC (100k events per mode/resonance combination) ###
exMC_psi2S_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psi2S_ppbar_etap_gpipi"
  config.related_dataset = psi2S_data
  config.events          = 100000
  config.decay_card      = decay_card_psi2S_modeI
  config.cross_section   = :default
end

exMC_jpsi_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_ppbar_etap_gpipi"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_jpsi_modeI
  config.cross_section   = :default
end

exMC_psi2S_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psi2S_ppbar_etap_etapipi"
  config.related_dataset = psi2S_data
  config.events          = 100000
  config.decay_card      = decay_card_psi2S_modeII
  config.cross_section   = :default
end

exMC_jpsi_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_ppbar_etap_etapipi"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_jpsi_modeII
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# --- Mode I : p pbar gamma pi+ pi-  (eta' -> gamma pi+ pi-) ---
selection_modeI = Selection.new
  .select_track {
    cos_theta 0.93      # |cos(theta)| < 0.93
    Vz        10.0      # |Vz| < 10 cm along the beam axis
    Vr        1.0       # Vr < 1 cm in the transverse plane
    nChrp     "==2"     # exactly two positive tracks
    nChrn     "==2"     # exactly two negative tracks
    nNet      "==0"     # net charge zero
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0    # >= 10 deg away from any charged track
    energyThreshold_b 0.025   # E > 25 MeV (barrel)
    energyThreshold_e 0.050   # E > 50 MeV (endcap)
    nGam              ">=1"   # at least one photon
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion,   against: [:kaon, :proton]   # pi+/pi- vs K and p
    identify :proton, against: [:pion, :kaon]     # p/pbar vs pi and K
    npip "==1"
    npim "==1"
    nprp "==1"
    nprm "==1"
  }
  .kinematic_fit([:prp, :prm, :pip, :pim, :gamma]) {   # 4C fit to p pbar pi+ pi- gamma
    nominal
    constrain_four_momentum
    chi2_cut 200                                       # loose cut; tight cut determined in ROOT
  }

# --- Mode II : p pbar eta pi+ pi-  (eta' -> eta pi+ pi-, eta -> gamma gamma) ---
selection_modeII = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"   # at least two photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion,   against: [:kaon, :proton]
    identify :proton, against: [:pion, :kaon]
    npip "==1"
    npim "==1"
    nprp "==1"
    nprm "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct eta from two photons (1-C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:prp, :prm, :pip, :pim, :eta]) {   # 4C fit to p pbar pi+ pi- eta
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

### Algorithms ###
# psi(3686) Mode I
alg_psi2S_modeI = Algorithm.new("Psi2SppbarEtapToGPiPi")
alg_psi2S_modeI.set_header(["Psi2SppbarEtapToGPiPiAlg/Psi2SppbarEtapToGPiPi.h"])
               .set_constant({"ECMS" => [:double, 3.686]})
               .set_alias({"std::vector<double>" => "Vdouble"})
               .note(:background_veto,
                 "At psi(3686) energy, veto psi(3686) -> gamma chi_cJ, psi(3686) -> pi+ pi- J/psi and psi(3686) -> eta J/psi backgrounds via recoil-mass windows; window boundaries fixed from MC study and applied in the ROOT analysis after the 4C fit.")
alg_psi2S_modeI.with_decay_card(decay_card_psi2S_modeI).apply(selection_modeI.dup)

# J/psi Mode I
alg_jpsi_modeI = Algorithm.new("JpsippbarEtapToGPiPi")
alg_jpsi_modeI.set_header(["JpsippbarEtapToGPiPiAlg/JpsippbarEtapToGPiPi.h"])
              .set_constant({"ECMS" => [:double, 3.097]})
              .set_alias({"std::vector<double>" => "Vdouble"})
alg_jpsi_modeI.with_decay_card(decay_card_jpsi_modeI).apply(selection_modeI.dup)

# psi(3686) Mode II
alg_psi2S_modeII = Algorithm.new("Psi2SppbarEtapToEtaPiPi")
alg_psi2S_modeII.set_header(["Psi2SppbarEtapToEtaPiPiAlg/Psi2SppbarEtapToEtaPiPi.h"])
                .set_constant({"ECMS" => [:double, 3.686]})
                .set_alias({"std::vector<double>" => "Vdouble"})
                .note(:background_veto,
                  "At psi(3686) energy, veto psi(3686) -> gamma chi_cJ, psi(3686) -> pi+ pi- J/psi and psi(3686) -> eta J/psi backgrounds via recoil-mass windows; window boundaries fixed from MC study and applied in the ROOT analysis after the 4C fit.")
alg_psi2S_modeII.with_decay_card(decay_card_psi2S_modeII).apply(selection_modeII.dup)

# J/psi Mode II
alg_jpsi_modeII = Algorithm.new("JpsippbarEtapToEtaPiPi")
alg_jpsi_modeII.set_header(["JpsippbarEtapToEtaPiPiAlg/JpsippbarEtapToEtaPiPi.h"])
               .set_constant({"ECMS" => [:double, 3.097]})
               .set_alias({"std::vector<double>" => "Vdouble"})
alg_jpsi_modeII.with_decay_card(decay_card_jpsi_modeII).apply(selection_modeII.dup)

### Execution on real data, matching inclusive MC and exclusive signal MC ###
root_files = []
root_files += alg_psi2S_modeI.execute_on([psi2S_data, psi2S_incMC, exMC_psi2S_modeI])
root_files += alg_jpsi_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_modeI])
root_files += alg_psi2S_modeII.execute_on([psi2S_data, psi2S_incMC, exMC_psi2S_modeII])
root_files += alg_jpsi_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_modeII])