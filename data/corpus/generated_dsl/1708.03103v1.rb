# ============================================================
# psi(3686) radiative decays: psi(3686) -> gamma eta', gamma eta, gamma pi0
# Five independent signal modes -> five Algorithm objects (Rule T1)
# ============================================================

### Datasets ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC

### Decay cards (EvtGen syntax) ###
# Mode I: psi(2S) -> gamma eta', eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_modeI = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.000 pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: psi(2S) -> gamma eta', eta' -> pi0 pi0 eta, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.000 pi0 pi0 eta PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode III: psi(2S) -> gamma eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_modeIII = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta PHSP;
  Enddecay

  Decay eta
  1.000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode IV: psi(2S) -> gamma eta, eta -> pi0 pi0 pi0
decay_card_modeIV = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma eta PHSP;
  Enddecay

  Decay eta
  1.000 pi0 pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode V: psi(2S) -> gamma pi0, pi0 -> gamma gamma
decay_card_modeV = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC samples (200k events each, five modes) ###
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_etap_pipim_eta"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_etap_pi0pi0_eta"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

exMC_modeIII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_eta_pipim_pi0"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_modeIII
  config.cross_section   = :default
end

exMC_modeIV = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_eta_pi0pi0pi0"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_modeIV
  config.cross_section   = :default
end

exMC_modeV = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_pi0"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_modeV
  config.cross_section   = :default
end

# ============================================================
### Event selection (BOSS) ###
# ============================================================

# ------------------------------------------------------------
# Mode I: psi(2S) -> gamma eta', eta' -> pi+ pi- eta, eta -> gamma gamma
# charged final state  gamma pi+ pi- gamma gamma
# ------------------------------------------------------------
alg_modeI = Algorithm.new("EtapPipimEta")
alg_modeI.set_header(["EtapPipimEtaAlg/EtapPipimEta.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_modeI = Selection.new
sel_modeI.select_track {                    # exactly one pi+ and one pi-
            cos_theta 0.93                  # |cos(theta)| < 0.93
            Vz        10.0                  # |Vz| < 10 cm
            Vr        1.0                   # Vr < 1 cm
            nChrp     "==1"                 # exactly one positive track
            nChrn     "==1"                 # exactly one negative track
            nNet      "==0"                 # net charge zero
          }
         .select_photon {                   # >= 3 photons; angle to nearest track > 10 deg (mode I only)
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025         # E > 25 MeV (barrel)
            energyThreshold_e 0.050         # E > 50 MeV (endcap)
            angle_to_track    10.0          # > 10 deg from any charged track (eta'->pi+pi-eta only)
            nGam              ">=3"
          }
         .pid(method: :probability) {       # pion PID, probability method
            prob_cut 0.001
            identify :pion, against: [:kaon, :proton]
          }
         .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {   # 4C fit: recoil gamma + 2 gamma (eta) + pi+ pi-
            nominal
            constrain_four_momentum
            chi2_cut 80
          }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

# ------------------------------------------------------------
# Mode III: psi(2S) -> gamma eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma
# charged final state  gamma pi+ pi- gamma gamma
# ------------------------------------------------------------
alg_modeIII = Algorithm.new("EtaPipimPi0")
alg_modeIII.set_header(["EtaPipimPi0Alg/EtaPipimPi0.h"])
           .set_constant({"ECMS" => [:double, 3.686]})

sel_modeIII = Selection.new
sel_modeIII.select_track {                  # exactly one pi+ and one pi-
             cos_theta 0.93
             Vz        10.0
             Vr        1.0
             nChrp     "==1"
             nChrn     "==1"
             nNet      "==0"
           }
           .select_photon {                 # >= 3 photons (no extra angle requirement for this channel)
             tdc_emc_start     0
             tdc_emc_end       14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             nGam              ">=3"
           }
           .pid(method: :probability) {     # pion PID, probability method
             prob_cut 0.001
             identify :pion, against: [:kaon, :proton]
           }
           .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {   # 4C fit: recoil gamma + 2 gamma (pi0) + pi+ pi-
             nominal
             constrain_four_momentum
             chi2_cut 80
           }

alg_modeIII.with_decay_card(decay_card_modeIII).apply(sel_modeIII)

# ------------------------------------------------------------
# Mode II: psi(2S) -> gamma eta', eta' -> pi0 pi0 eta, eta -> gamma gamma
# neutral final state  7 gamma ; charged tracks vetoed
# ------------------------------------------------------------
alg_modeII = Algorithm.new("EtapPi0Pi0Eta")
alg_modeII.set_header(["EtapPi0Pi0EtaAlg/EtapPi0Pi0Eta.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_modeII = Selection.new
sel_modeII.select_track {                   # veto all charged tracks
            nChrp "==0"
            nChrn "==0"
          }
          .select_photon {                  # 7 photons
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam              ">=7"
          }
          .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) {   # 4C fit over 7 gamma
            nominal
            constrain_four_momentum
            chi2_cut 80
          }

# inexpressible: custom photon-combination figure of merit (chi2_M) retained for eta'->pi0pi0eta
alg_modeII.note(:custom_photon_combination,
  "retain only the photon combination minimising chi2_M = sum over pi0 candidates of " \
  "(M_gg - M_pi0)^2 / sigma_pi0^2 + (M_gg - M_eta)^2 / sigma_eta^2, with sigma_pi0 = 4.8 MeV, " \
  "sigma_eta = 8.7 MeV; the eta'->pi0pi0eta channel keeps the combination with the smallest chi2_M")

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

# ------------------------------------------------------------
# Mode IV: psi(2S) -> gamma eta, eta -> pi0 pi0 pi0
# neutral final state  7 gamma ; charged tracks vetoed
# ------------------------------------------------------------
alg_modeIV = Algorithm.new("EtaPi0Pi0Pi0")
alg_modeIV.set_header(["EtaPi0Pi0Pi0Alg/EtaPi0Pi0Pi0.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_modeIV = Selection.new
sel_modeIV.select_track {                   # veto all charged tracks
            nChrp "==0"
            nChrn "==0"
          }
          .select_photon {                  # 7 photons
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam              ">=7"
          }
          .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) {   # 4C fit over 7 gamma
            nominal
            constrain_four_momentum
            chi2_cut 80
          }

alg_modeIV.with_decay_card(decay_card_modeIV).apply(sel_modeIV)

# ------------------------------------------------------------
# Mode V: psi(2S) -> gamma pi0, pi0 -> gamma gamma
# neutral final state  3 gamma (barrel only) ; charged tracks vetoed
# ------------------------------------------------------------
alg_modeV = Algorithm.new("GamPi0")
alg_modeV.set_header(["GamPi0Alg/GamPi0.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_modeV = Selection.new
sel_modeV.select_track {                    # veto all charged tracks
            nChrp "==0"
            nChrn "==0"
          }
          .select_photon {                  # exactly 3 photons
            tdc_emc_start     0
            tdc_emc_end       14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam              "==3"
          }
          .kinematic_fit([:gamma, :gamma, :gamma]) {   # 4C fit over 3 gamma (tighter chi2)
            nominal
            constrain_four_momentum
            chi2_cut 40
          }

# inexpressible suppression criteria for the gamma pi0 channel
alg_modeV
  .note(:barrel_photon_requirement,
    "the three accepted photons are required to lie in the EMC barrel region, |cos(theta)| < 0.8")
  .note(:pi0_helicity_cut,
    "additional suppression: require |cos(theta_hel)| < 0.7 of the pi0 candidate")
  .note(:background_veto,
    "veto converted photons from e+e- -> gamma gamma (gamma_ISR): reject events with fewer than 8 " \
    "MDC hits between the interaction point and the two EMC shower positions")

alg_modeV.with_decay_card(decay_card_modeV).apply(sel_modeV)

# ============================================================
### Execute on datasets ###
# ============================================================
root_files_modeI   = alg_modeI.execute_on([psip_data, psip_incMC, exMC_modeI])
root_files_modeII  = alg_modeII.execute_on([psip_data, psip_incMC, exMC_modeII])
root_files_modeIII = alg_modeIII.execute_on([psip_data, psip_incMC, exMC_modeIII])
root_files_modeIV  = alg_modeIV.execute_on([psip_data, psip_incMC, exMC_modeIV])
root_files_modeV   = alg_modeV.execute_on([psip_data, psip_incMC, exMC_modeV])