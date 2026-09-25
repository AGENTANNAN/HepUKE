# ============================================================================
# BESIII J/psi analysis:
#   J/psi -> phi pi0 f0(980) ,  J/psi -> phi f1(1285) ,  J/psi -> phi eta'
#   phi -> K+ K- ,  f0(980) -> pi+ pi- / pi0 pi0 ,  f1 -> pi0 f0(980)
#   eta' -> pi+ pi- pi0 / 3 pi0
# Charged final state : K+ K- pi+ pi- pi0  (fitted as K+ K- pi+ pi- gamma gamma)
# Neutral final state : K+ K- pi0 pi0 pi0  (fitted as K+ K- + six photons)
# ============================================================================

### ------------------------------- Datasets -------------------------------- ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi 2009+2012 data (1.311e9 events)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # inclusive J/psi MC  (1.2e9 events)

### -------------------------- Decay cards (6 modes) ------------------------ ###

# Mode 1 : J/psi -> phi pi0 f0(980), f0 -> pi+ pi-   (charged final state)
decay_card_m1 = <<~DECAYCARD
  Decay J/psi
  1.000 phi pi0 f_0 PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay f_0
  1.000 pi+ pi- PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode 2 : J/psi -> phi pi0 f0(980), f0 -> pi0 pi0   (neutral final state)
decay_card_m2 = <<~DECAYCARD
  Decay J/psi
  1.000 phi pi0 f_0 PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay f_0
  1.000 pi0 pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode 3 : J/psi -> phi f1(1285), f1 -> pi0 f0(980), f0 -> pi+ pi-   (charged final state)
decay_card_m3 = <<~DECAYCARD
  Decay J/psi
  1.000 phi f_1 PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay f_1
  1.000 pi0 f_0 PHSP;
  Enddecay
  Decay f_0
  1.000 pi+ pi- PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode 4 : J/psi -> phi f1(1285), f1 -> pi0 f0(980), f0 -> pi0 pi0   (neutral final state)
decay_card_m4 = <<~DECAYCARD
  Decay J/psi
  1.000 phi f_1 PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay f_1
  1.000 pi0 f_0 PHSP;
  Enddecay
  Decay f_0
  1.000 pi0 pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode 5 : J/psi -> phi eta', eta' -> pi+ pi- pi0   (charged final state)
decay_card_m5 = <<~DECAYCARD
  Decay J/psi
  1.000 phi eta' PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay eta'
  1.000 pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode 6 : J/psi -> phi eta', eta' -> 3 pi0   (neutral final state)
decay_card_m6 = <<~DECAYCARD
  Decay J/psi
  1.000 phi eta' PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay eta'
  1.000 pi0 pi0 pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

### ------------------- Exclusive MC samples (200k events each) -------------- ###
exMC_m1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "xmc_jpsi_phiPi0f0_pipi"       # phi pi0 f0, f0 -> pi+ pi-
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_m1
  config.cross_section   = :default
end

exMC_m2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "xmc_jpsi_phiPi0f0_pi0pi0"     # phi pi0 f0, f0 -> pi0 pi0
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_m2
  config.cross_section   = :default
end

exMC_m3 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "xmc_jpsi_phif1_pipi"          # phi f1, f0 -> pi+ pi-
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_m3
  config.cross_section   = :default
end

exMC_m4 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "xmc_jpsi_phif1_pi0pi0"        # phi f1, f0 -> pi0 pi0
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_m4
  config.cross_section   = :default
end

exMC_m5 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "xmc_jpsi_phietap_pipipi0"     # phi eta', eta' -> pi+ pi- pi0
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_m5
  config.cross_section   = :default
end

exMC_m6 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "xmc_jpsi_phietap_3pi0"        # phi eta', eta' -> 3 pi0
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_m6
  config.cross_section   = :default
end

### ========================= Charged channel (K+ K- pi+ pi- pi0) ========== ###
# Modes 1, 3, 5 share the charged final state  ->  one common selection.
charged_selection = Selection.new
charged_selection
  .select_track {                       # charged track quality + multiplicity
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     "==2"                     # two positive tracks
    nChrn     "==2"                     # two negative tracks
    nNet      "==0"                     # net charge zero
  }
  .select_photon {                      # photon selection
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0              # >= 10 deg from any charged track
    energyThreshold_b 0.025             # 25 MeV (barrel)
    energyThreshold_e 0.050             # 50 MeV (endcap)
    nGam              ">=2"             # at least two photons (pi0 -> gamma gamma)
  }
  .pid(method: :probability) {          # PID (probability method), K/pi separation
    prob_cut 0.0                        # probability cut 0.0
    identify :kaon, against: [:pion]    # K+ and K-  vs  pions
    nkp ">=1"                           # at least one K+
    nkm ">=1"                           # at least one K-
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])       # remove kaons from charged-track samples
  .assign({:chrgp => :pip, :chrgn => :pim})     # remaining tracks assigned as pi+ / pi-
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {   # 5C fit
    nominal                                     # nominal fit
    constrain_four_momentum                     # 4C four-momentum conservation
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # + 1C pi0 mass
    chi2_cut 100                                # chi^2 < 100
  }

alg_charged = Algorithm.new("JPsiPhiCharged")
alg_charged.set_header(["JPsiPhiChargedAlg/JPsiPhiCharged.h"])
           .set_constant({"ECMS" => [:double, 3.097]})
           .note(:background_veto,
                 "K*0 veto: reject events with |M(K+ pi-) - M(K*0)| < 0.050 GeV/c^2 " \
                 "(equivalently the K- pi+ combination) to suppress K*0 -> K pi background; " \
                 "M(K*0)=0.8955 GeV/c^2, veto window [0.8455, 0.9455] GeV/c^2.")
           .note(:phi_selection,
                 "the K+ K- pair whose invariant mass is closest to the phi(1020) mass is " \
                 "selected as the phi candidate before the kinematic fit.")
# Decay card for the charged channel (defines the K+ K- pi+ pi- gamma gamma variables
# common to modes 1, 3, 5).
alg_charged.with_decay_card(decay_card_m1).apply(charged_selection)
# Charged channel runs on the three charged-mode exclusive samples
alg_charged.execute_on([jpsi_data, jpsi_incMC, exMC_m1, exMC_m3, exMC_m5])

### ========================= Neutral channel (K+ K- pi0 pi0 pi0) =========== ###
# Modes 2, 4, 6 share the neutral final state  ->  one common selection.
neutral_selection = Selection.new
neutral_selection
  .select_track {                       # charged track quality + multiplicity
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     "==1"                     # one positive track
    nChrn     "==1"                     # one negative track
    nNet      "==0"                     # net charge zero
  }
  .select_photon {                      # photon selection
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0              # >= 10 deg from any charged track
    energyThreshold_b 0.025             # 25 MeV (barrel)
    energyThreshold_e 0.050             # 50 MeV (endcap)
    nGam              ">=6"             # at least six photons (three pi0 -> gamma gamma)
  }
  .pid(method: :probability) {          # PID (probability method), K/pi separation
    prob_cut 0.0
    identify :kaon, against: [:pion]    # K+ and K-
    nkp ">=1"                           # at least one K+
    nkm ">=1"                           # at least one K-
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])       # remove kaons from charged-track samples
  .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) {  # 7C fit
    nominal                                     # nominal fit
    constrain_four_momentum                     # 4C four-momentum conservation
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # + 1C
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # + 1C
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # + 1C
    chi2_cut 90                                 # chi^2 < 90
  }

alg_neutral = Algorithm.new("JPsiPhiNeutral")
alg_neutral.set_header(["JPsiPhiNeutralAlg/JPsiPhiNeutral.h"])
           .set_constant({"ECMS" => [:double, 3.097]})
           .note(:phi_selection,
                 "the K+ K- pair whose invariant mass is closest to the phi(1020) mass is " \
                 "selected as the phi candidate before the kinematic fit.")
           .note(:pi0_pairing,
                 "the six photons are paired into three pi0 candidates by minimising the " \
                 "pi0 mass chi^2 before the 7C kinematic fit.")
# Decay card for the neutral channel (defines the K+ K- + six photon variables
# common to modes 2, 4, 6).
alg_neutral.with_decay_card(decay_card_m2).apply(neutral_selection)
# Neutral channel runs on the three neutral-mode exclusive samples
alg_neutral.execute_on([jpsi_data, jpsi_incMC, exMC_m2, exMC_m4, exMC_m6])