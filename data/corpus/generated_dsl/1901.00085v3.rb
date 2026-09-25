# ============================================================
# Dataset preparation — J/psi(3097)
# ============================================================
jpsi_data  = DatasetManager.real_data.find("708_3097")      # 1.3e9-event J/psi real data set
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

# Decay card — Mode I: J/psi -> phi eta eta', phi -> K+K-, eta -> gamma gamma, eta' -> gamma pi+pi-
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta etap        PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-               VSS;
    Enddecay

    Decay eta
    1.0000 gamma gamma         PHSP;
    Enddecay

    Decay etap
    1.0000 gamma pi+ pi-       PHSP;
    Enddecay

    End
DECAYCARD

# Decay card — Mode II: J/psi -> phi eta eta', phi -> K+K-, eta' -> eta pi+pi-, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta etap        PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K-               VSS;
    Enddecay

    Decay etap
    1.0000 eta pi+ pi-         PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma         PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC samples, one per mode
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_phi_eta_etap_modeI"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_phi_eta_etap_modeII"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# ============================================================
# Event selection (BOSS) — Mode I: K+K-pi+pi-gamma gamma gamma
# ============================================================
alg_name_modeI = "JpsiPhiEtaEtapModeI"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})            # J/psi center-of-mass energy
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
  .select_track {                                              # Charged track selection
    cos_theta   0.93                                           # |cos(theta)| < 0.93
    Vz          10.0                                           # |Vz| < 10 cm
    Vr          1.0                                            # Vr < 1 cm
    nChrp       "==2"                                          # exactly 2 positive tracks
    nChrn       "==2"                                          # exactly 2 negative tracks
    nNet        "==0"                                          # net charge zero
  }
  .select_photon {                                             # Photon selection
    tdc_emc_start   0
    tdc_emc_end     14
    energyThreshold_b 0.025                                    # barrel E > 25 MeV
    energyThreshold_e 0.050                                    # endcap E > 50 MeV
    angle_to_track  10.0                                       # angle to nearest charged track > 10 deg
    nGam            ">=3"                                      # at least 3 photons (Mode I)
  }
  .pid(method: :probability) {                                 # PID: probability method
    prob_cut  0.001
    identify :kaon, against: [:pion, :proton]                  # K+ and K- vs pi/p
    identify :pion, against: [:kaon, :proton]                  # pi+ and pi- vs K/p
    nkp   "==1"                                                # exactly one K+
    nkm   "==1"                                                # exactly one K-
    npip  "==1"                                                # exactly one pi+
    npim  "==1"                                                # exactly one pi-
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma]) {  # nominal 4C fit, K+K-pi+pi-g g g
    nominal
    constrain_four_momentum
    chi2_cut 40                                                # chi2 < 40
    invariant_mass_of(:gamma, :gamma).out_of(0.12, 0.15)       # reject pi0-like gamma-gamma pairs
    invariant_mass_of(:gamma, :gamma).within(0.522, 0.573)     # eta window
    invariant_mass_of(:kp, :km).within(1.010, 1.030)           # phi window
    invariant_mass_of(:gamma, :pip, :pim).within(0.936, 0.979) # eta' window (gamma pi+ pi-)
    invariant_mass_of(:pip, :pim).smaller_than(0.87)           # M(pi+pi-) < 0.87 GeV
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

# ============================================================
# Event selection (BOSS) — Mode II: K+K-pi+pi-gamma gamma gamma gamma
# ============================================================
alg_name_modeII = "JpsiPhiEtaEtapModeII"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
  .select_track {                                              # Charged track selection
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nChrp       "==2"
    nChrn       "==2"
    nNet        "==0"
  }
  .select_photon {                                             # Photon selection
    tdc_emc_start   0
    tdc_emc_end     14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track  10.0
    nGam            ">=4"                                      # at least 4 photons (Mode II)
  }
  .pid(method: :probability) {                                 # PID: probability method
    prob_cut  0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp   "==1"
    nkm   "==1"
    npip  "==1"
    npim  "==1"
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma, :gamma]) {  # nominal 4C fit
    nominal
    constrain_four_momentum
    chi2_cut 80                                                # chi2 < 80
    invariant_mass_of(:gamma, :gamma).out_of(0.12, 0.15)       # reject pi0-like gamma-gamma pairs
    invariant_mass_of(:gamma, :gamma).within(0.509, 0.586)     # eta window
    invariant_mass_of(:kp, :km).within(1.010, 1.030)           # phi window
    invariant_mass_of(:gamma, :gamma, :pip, :pim).within(0.920, 0.995)  # eta' window (eta pi+ pi-)
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

# ============================================================
# Execution
# ============================================================
root_files_modeI  = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])
root_files_modeII = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])