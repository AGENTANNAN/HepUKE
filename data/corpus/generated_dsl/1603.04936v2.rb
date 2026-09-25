# ============================================================================
# Search for h_c -> gamma eta' and h_c -> gamma eta via psi(2S) -> pi0 h_c
# BOSS part: dataset preparation + event selection up to the final kinematic fit
# ============================================================================

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(2S) data (2009 + 2012, ~4.48e8 events)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # matching inclusive MC
cont_data  = DatasetManager.real_data.find("709_3650")     # continuum, 44 pb^-1 at sqrt(s) = 3.65 GeV

# ------------------------------- Decay cards -------------------------------
# Mode 1: psi' -> pi0 h_c ; h_c -> gamma eta' ; eta' -> pi+ pi- eta ; eta -> gamma gamma
decay_card_mode1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 2: psi' -> pi0 h_c ; h_c -> gamma eta' ; eta' -> gamma pi+ pi-
decay_card_mode2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 3: psi' -> pi0 h_c ; h_c -> gamma eta ; eta -> gamma gamma (all-neutral)
decay_card_mode3 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 4: psi' -> pi0 h_c ; h_c -> gamma eta ; eta -> pi+ pi- pi0 ; pi0 -> gamma gamma
decay_card_mode4 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta PHSP;
    Enddecay

    Decay eta
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --------------------------- Exclusive MC samples --------------------------
# 200k events for each of the four signal chains
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_etap_pipimeta_etagg"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_mode1
  config.cross_section   = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_etap_gampipi"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_mode2
  config.cross_section   = :default
end

exMC_mode3 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_eta_gg"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_mode3
  config.cross_section   = :default
end

exMC_mode4 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0hc_eta_pipipi0"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = decay_card_mode4
  config.cross_section   = :default
end

### Event selection ###
# Common kinematics: |cos(theta)| < 0.93, |Vz| < 10 cm, Vr < 1 cm;
# photons: EMC timing 0-700 ns, >= 10 deg from nearest charged track,
# E > 25 MeV (barrel) / > 50 MeV (endcap); PID prob > 0.001 separating pi vs K, p.

# ---------------------------------------------------------------------------
# Mode 1: pi+ pi- eta gamma gamma   ->  6C fit (4C + m(pi0) + m(eta)), chi2 < 120
# ---------------------------------------------------------------------------
alg_name_mode1 = "HcEtapPipPimEtaGG"
alg_mode1 = Algorithm.new(alg_name_mode1)
alg_mode1.set_header(["#{alg_name_mode1}Alg/#{alg_name_mode1}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_mode1 = Selection.new
  .select_track {                     # charged track selection
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==1"                 # exactly one positive track
      nChrn     "==1"                 # exactly one negative track
      nNet      "==0"                 # net charge zero
  }
  .select_photon {                    # photon selection
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=5"         # pi0->gg (2) + h_c->g (1) + eta->gg (2)
  }
  .pid(method: :probability) {        # pion identification
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]
      npip "==1"
      npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # pi0 -> gamma gamma (1C, chi2 < 25)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # eta -> gamma gamma (1C, chi2 < 25)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 25
      neta ">=1"
  }
  .kinematic_fit([:pip, :pim, :gamma, :eta, :pi0]) {  # 6C fit
      nominal                      # nominal fit (corrected 4-momenta saved)
      vertex_fit([0, 1])           # constrain pi+ and pi- to a common vertex
      constrain_four_momentum      # 4C
      chi2_cut 120                 # loose cut; tight signal region applied in ROOT
      # all photon/particle combinations are looped over, smallest chi2 retained
  }

alg_mode1.with_decay_card(decay_card_mode1).apply(sel_mode1)
alg_mode1.execute_on([psip_data, psip_incMC, cont_data, exMC_mode1])

# ---------------------------------------------------------------------------
# Mode 2: gamma pi+ pi- eta'   ->  5C fit (4C + m(pi0)), chi2 < 50
# ---------------------------------------------------------------------------
alg_name_mode2 = "HcEtapGamPiPi"
alg_mode2 = Algorithm.new(alg_name_mode2)
alg_mode2.set_header(["#{alg_name_mode2}Alg/#{alg_name_mode2}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_mode2 = Selection.new
  .select_track {
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
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=4"         # pi0->gg (2) + h_c->g (1) + eta'->g (1)
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]
      npip "==1"
      npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # pi0 -> gamma gamma (1C, chi2 < 25)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"
  }
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :pi0]) {  # 5C fit
      nominal
      vertex_fit([0, 1])           # constrain pi+ and pi- to a common vertex
      constrain_four_momentum      # 4C
      chi2_cut 50
      # smallest-chi2 combination retained
  }

alg_mode2.with_decay_card(decay_card_mode2).apply(sel_mode2)
alg_mode2.execute_on([psip_data, psip_incMC, cont_data, exMC_mode2])

# ---------------------------------------------------------------------------
# Mode 3: all-neutral (eta -> gamma gamma)   ->  6C fit, chi2 < 200
# ---------------------------------------------------------------------------
alg_name_mode3 = "HcEtaGG"
alg_mode3 = Algorithm.new(alg_name_mode3)
alg_mode3.set_header(["#{alg_name_mode3}Alg/#{alg_name_mode3}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_mode3 = Selection.new
  .select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==0"                 # all-neutral mode: no charged tracks
      nChrn     "==0"
      nNet      "==0"
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=5"         # pi0->gg (2) + h_c->g (1) + eta->gg (2)
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # pi0 -> gamma gamma (1C, chi2 < 25)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # eta -> gamma gamma (1C, chi2 < 25)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 25
      neta ">=1"
  }
  .kinematic_fit([:gamma, :pi0, :eta]) {        # 6C fit (4C + m(pi0) + m(eta))
      nominal
      constrain_four_momentum
      chi2_cut 200
      # no charged tracks -> no vertex constraint; smallest-chi2 combination kept
  }

alg_mode3.with_decay_card(decay_card_mode3).apply(sel_mode3)
alg_mode3.execute_on([psip_data, psip_incMC, cont_data, exMC_mode3])

# ---------------------------------------------------------------------------
# Mode 4: eta -> pi+ pi- pi0   ->  6C fit (4C + m(pi0) + m(pi0)), chi2 < 120
# ---------------------------------------------------------------------------
alg_name_mode4 = "HcEtaPiPiPi0"
alg_mode4 = Algorithm.new(alg_name_mode4)
alg_mode4.set_header(["#{alg_name_mode4}Alg/#{alg_name_mode4}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_mode4 = Selection.new
  .select_track {
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
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=5"         # pi0->gg (2) + h_c->g (1) + pi0(eta)->gg (2)
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]
      npip "==1"
      npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {     # pi0 -> gamma gamma (two pi0 required)
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=2"                    # two pi0 candidates (psi'->pi0 and eta->pi0)
  }
  .kinematic_fit([:pip, :pim, :gamma, :pi0, :pi0]) {   # 6C fit (4C + 2 x m(pi0))
      nominal
      vertex_fit([0, 1])            # constrain pi+ and pi- to a common vertex
      constrain_four_momentum
      chi2_cut 120
      # smallest-chi2 combination retained
  }

alg_mode4.with_decay_card(decay_card_mode4).apply(sel_mode4)
alg_mode4.execute_on([psip_data, psip_incMC, cont_data, exMC_mode4])