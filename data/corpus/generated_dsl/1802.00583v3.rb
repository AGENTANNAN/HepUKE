# ============================================================================
# BOSS-side (dataset preparation + event selection) DSL for
#   Mode A1: J/psi -> phi eta pi0, phi -> K+K-,  eta -> gamma gamma   (3.097 GeV)
#   Mode A2: J/psi -> phi eta pi0, phi -> K+K-,  eta -> pi+pi-pi0     (3.097 GeV)
#   Mode B : psi(3686) -> gamma chi_c1,          chi_c1 -> pi0 pi+pi-  (3.686 GeV)
# ============================================================================

### ---------------------------------------------------------------- Datasets ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # 3.097 GeV J/psi real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # J/psi inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")     # 3.686 GeV psi(3686) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # psi(3686) inclusive MC

### ------------------------------------------------------------ Decay cards ###
# Mode A1: J/psi -> phi eta pi0, phi -> K+K-, eta -> gamma gamma
decay_card_A1 = <<~DECAYCARD
  Decay J/psi
  1.000 phi eta pi0 PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode A2: J/psi -> phi eta pi0, phi -> K+K-, eta -> pi+ pi- pi0
decay_card_A2 = <<~DECAYCARD
  Decay J/psi
  1.000 phi eta pi0 PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
  Decay eta
  1.000 pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode B: psi(2S) -> gamma chi_c1, chi_c1 -> pi0 pi+ pi-
decay_card_B = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c1 PHSP;
  Enddecay
  Decay chi_c1
  1.000 pi0 pi+ pi- PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

### --------------------------------------------- Exclusive MC (100k events) ###
exMC_A1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_phi_eta_pi0_eta2g"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_A1
  config.cross_section   = :default
end

exMC_A2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_phi_eta_pi0_eta3pi"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_A2
  config.cross_section   = :default
end

exMC_B = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chi_c1"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_B
  config.cross_section   = :default
end

# ============================================================================
# Mode A1 : J/psi -> phi eta pi0, phi -> K+K-, eta -> gamma gamma
# ============================================================================
alg_name_A1 = "JpsiToPhiEtaPi0EtaGG"
alg_A1 = Algorithm.new(alg_name_A1)
alg_A1.set_header(["#{alg_name_A1}Alg/#{alg_name_A1}.h"])
      .set_constant({ "ECMS" => [:double, 3.097] })   # CM energy of the J/psi data

sel_A1 = Selection.new
  .select_track {                 # charged-track quality
    cos_theta 0.93                # |cos theta| < 0.93
    Vz        10.0                # |Vz| < 10 cm
    Vr        1.0                 # Vr < 1 cm
    nNet      "==0"               # net charge zero
  }
  .select_photon {                # photon selection
    angle_to_track    10.0        # >= 10 deg from nearest charged track
    energyThreshold_b 0.025       # > 25 MeV in the barrel
    energyThreshold_e 0.050       # > 50 MeV in the endcap
    nGam              ">=4"       # at least four photons (both pi0/eta -> gamma gamma)
  }
  .pid(method: :probability) {    # PID: probability method
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]  # K+ and K- (charge-conjugation shorthand)
    nkp ">=1"                     # at least one K+
    nkm ">=1"                     # at least one K-
  }
  # Nominal 4C kinematic fit to K+K-gammagammagammagamma with chi2 < 50
  .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 50
    invariant_mass_of(:gamma, :gamma).within(0.120, 0.150)  # |M(gammagamma)-m_pi0| < 15 MeV
    invariant_mass_of(:gamma, :gamma).within(0.518, 0.578)  # |M(gammagamma)-m_eta| < 30 MeV
    invariant_mass_of(:kp, :km).within(1.009, 1.029)        # |M(K+K-)-m_phi|     < 10 MeV
  }
  # Competing hypothesis J/psi -> K+K-pi0pi0 (store chi2 for ROOT-level veto chi2_pi0pi0 > 40)
  .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  }
  # Competing hypothesis J/psi -> K+K-etaeta (store chi2 for ROOT-level veto chi2_etaeta > 5)
  .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  }

alg_A1.with_decay_card(decay_card_A1).apply(sel_A1)
root_files_A1 = alg_A1.execute_on([jpsi_data, jpsi_incMC, exMC_A1])

# ============================================================================
# Mode A2 : J/psi -> phi eta pi0, phi -> K+K-, eta -> pi+ pi- pi0
# ============================================================================
alg_name_A2 = "JpsiToPhiEtaPi0Eta3Pi"
alg_A2 = Algorithm.new(alg_name_A2)
alg_A2.set_header(["#{alg_name_A2}Alg/#{alg_name_A2}.h"])
      .set_constant({ "ECMS" => [:double, 3.097] })

sel_A2 = Selection.new
  .select_track {                 # charged-track quality
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"               # at least two positive tracks
    nChrn     ">=2"               # at least two negative tracks
    nNet      "==0"               # net charge zero
  }
  .select_photon {                # photon selection
    angle_to_track    10.0        # >= 10 deg from nearest charged track
    energyThreshold_b 0.025       # > 25 MeV (barrel)
    energyThreshold_e 0.050       # > 50 MeV (endcap)
    nGam              ">=4"       # at least four photons
  }
  .pid(method: :probability) {    # PID: kaons and pions against each other and protons
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp   ">=1"
    nkm   ">=1"
    npip  ">=1"
    npim  ">=1"
  }
  # Reconstruct pi0 from a photon pair (needed to form M(pi+pi-pi0) for the eta window)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # Nominal 4C kinematic fit to K+K-pi+pi- + pi0 + gamma gamma with chi2 < 60
  .kinematic_fit([:kp, :km, :pip, :pim, :pi0, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 60
    invariant_mass_of(:gamma, :gamma).within(0.120, 0.150)    # pi0 window < 15 MeV
    invariant_mass_of(:pip, :pim, :pi0).within(0.528, 0.568)  # eta window < 20 MeV
    invariant_mass_of(:kp, :km).within(1.009, 1.029)          # phi window < 10 MeV
  }

alg_A2.with_decay_card(decay_card_A2).apply(sel_A2)
root_files_A2 = alg_A2.execute_on([jpsi_data, jpsi_incMC, exMC_A2])

# ============================================================================
# Mode B : psi(3686) -> gamma chi_c1, chi_c1 -> pi0 pi+ pi-
# ============================================================================
alg_name_B = "PsipToGamChiC1ToPi0PiPi"
alg_B = Algorithm.new(alg_name_B)
alg_B.set_header(["#{alg_name_B}Alg/#{alg_name_B}.h"])
     .set_constant({ "ECMS" => [:double, 3.686] })   # CM energy of the psi(3686) data

sel_B = Selection.new
  .select_track {                 # charged-track quality
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nNet      "==0"               # net charge zero
  }
  .select_photon {                # photon selection
    angle_to_track    10.0        # >= 10 deg from nearest charged track
    energyThreshold_b 0.025       # > 25 MeV (barrel)
    energyThreshold_e 0.050       # > 50 MeV (endcap)
    nGam              ">=3"       # at least three photons (radiative gamma + pi0 -> gamma gamma)
  }
  .pid(method: :probability) {    # PID: pions against kaons/protons
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=1"
    npim ">=1"
  }
  # Reconstruct pi0 from a photon pair (the pi0 mass window < 15 MeV)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # Nominal 4C fit to (gamma) pi+pi- pi0 with chi2 < 20
  .kinematic_fit([:gamma, :pip, :pim, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 20
    invariant_mass_of(:pip, :pim, :pi0).within(3.491, 3.531)  # |M(pi+pi-pi0)-m_chi_c1| < 20 MeV
  }
  # Competing-hypothesis fit pi+pi-gammagammagammagamma (extra photon); store chi2 for veto
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }
  # Competing-hypothesis fit pi+pi-gammagamma (no radiative photon); store chi2 for veto
  .kinematic_fit([:pip, :pim, :gamma, :gamma]) {
    constrain_four_momentum
  }

alg_B.with_decay_card(decay_card_B).apply(sel_B)
root_files_B = alg_B.execute_on([psip_data, psip_incMC, exMC_B])