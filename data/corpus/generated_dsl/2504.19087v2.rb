# ============================================================
#  ψ(3686) → γ χ_cJ  (J = 0,1,2),  χ_cJ → η η η'
#  Two independent chains by the η' decay mode:
#     Mode 1: η' → γ  π+π−
#     Mode 2: η' → η  π+π−
# ============================================================

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # matching inclusive MC

# ---------------- Decay card: mode 1 (η' → γ π+π−) ----------------
decay_card_mode1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 eta eta eta' PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta'
    1.000 gamma pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# ---------------- Decay card: mode 2 (η' → η π+π−) ----------------
decay_card_mode2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 eta eta eta' PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta'
    1.000 eta pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# ---------------- Exclusive signal MC (100k events per η' mode) ----------------
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_chic1_etap_gammapipi"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_mode1
  config.cross_section   = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_chic1_etap_etapipi"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_mode2
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ================= Mode 1: η' → γ π+π−  (nominal 6C fit) =================
alg_name_mode1 = "ChiC1EtaEtaEtapGamPiPi"
alg_mode1 = Algorithm.new(alg_name_mode1)
alg_mode1.set_header(["#{alg_name_mode1}Alg/#{alg_name_mode1}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

mode1_selection = Selection.new
  .select_track {                       # exactly two charged tracks
      cos_theta 0.93                    # |cosθ| < 0.93
      Vz        10.0                    # |Vz| < 10 cm
      Vr        1.0                     # Vr < 1 cm
      nChrp     "==1"                   # one positive track
      nChrn     "==1"                   # one negative track
      nNet      "==0"                   # net charge zero
  }
  .select_photon {                      # photon selection
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0            # minimum opening angle 10°
      energyThreshold_b 0.025           # ≥ 25 MeV in the barrel
      energyThreshold_e 0.050           # ≥ 50 MeV in the endcap
      nGam              ">=8"           # ≥ 8 photons in mode 1
  }
  .pid(method: :probability) {          # pion identification
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]   # π+ and π− vs K, p
      npip "==1"
      npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {                 # η → γγ, 1C fit to the η mass
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 200
      neta ">=3"                                            # ≥ 3 η candidates in mode 1
  }
  .kinematic_fit([:gamma, :eta, :eta, :eta, :pip, :pim]) {  # nominal 6C fit
      nominal
      constrain_four_momentum
      chi2_cut 200
  }
  .kinematic_fit([:eta, :eta, :eta, :pip, :pim]) {          # competing hyp. (one fewer photon)
      constrain_four_momentum
  }
  .kinematic_fit([:eta, :eta, :pip, :pim]) {                # competing hyp. (two fewer photons)
      constrain_four_momentum
  }

alg_mode1
  .note(:background_veto, "π0 veto: every γγ pair with |M_γγ − m_π0| < 15 MeV is rejected")
  .note(:etap_mass_window, "η' candidate required within 30 MeV of the nominal η' mass")
alg_mode1.with_decay_card(decay_card_mode1).apply(mode1_selection)
root_files_mode1 = alg_mode1.execute_on([psip_data, psip_incMC, exMC_mode1])

# ================= Mode 2: η' → η π+π−  (nominal 7C fit) =================
alg_name_mode2 = "ChiC1EtaEtaEtapEtaPiPi"
alg_mode2 = Algorithm.new(alg_name_mode2)
alg_mode2.set_header(["#{alg_name_mode2}Alg/#{alg_name_mode2}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

mode2_selection = Selection.new
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
      nGam              ">=9"           # ≥ 9 photons in mode 2
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]
      npip "==1"
      npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {                        # η → γγ, 1C fit to the η mass
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 200
      neta ">=4"                                                   # ≥ 4 η candidates in mode 2
  }
  .kinematic_fit([:gamma, :eta, :eta, :eta, :eta, :pip, :pim]) {   # nominal 7C fit
      nominal
      constrain_four_momentum
      chi2_cut 200
  }
  .kinematic_fit([:eta, :eta, :eta, :eta, :pip, :pim]) {           # competing hyp. (one fewer photon)
      constrain_four_momentum
  }
  .kinematic_fit([:eta, :eta, :eta, :pip, :pim]) {                 # competing hyp. (two fewer photons)
      constrain_four_momentum
  }

alg_mode2
  .note(:background_veto, "π0 veto: every γγ pair with |M_γγ − m_π0| < 15 MeV is rejected")
  .note(:etap_mass_window, "η' candidate required within 30 MeV of the nominal η' mass")
alg_mode2.with_decay_card(decay_card_mode2).apply(mode2_selection)
root_files_mode2 = alg_mode2.execute_on([psip_data, psip_incMC, exMC_mode2])