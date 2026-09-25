# ============================================================
# J/psi -> gamma eta(1405), eta(1405) -> f0(980) pi0
#   Mode I : f0(980) -> pi+ pi-   (final state gamma pi+ pi- pi0 -> 3gamma pi+ pi-)
#   Mode II: f0(980) -> pi0 pi0   (final state gamma 3pi0 -> 7gamma)
# plus J/psi -> gamma eta' control/measurement modes (eta' -> pi+pi-pi0 / 3pi0)
# and the eta' peaking backgrounds (eta' -> gamma rho0, eta' -> gamma omega)
# at sqrt(s) = 3.097 GeV (J/psi).
# ============================================================

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # 225.2M-event J/psi real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # inclusive MC (Lund-Charm, 2.25e8)

# ---------------- Decay cards (EvtGen) ----------------

# Mode I signal: J/psi -> gamma eta(1405), eta(1405) -> f0(980) pi0, f0(980) -> pi+ pi-
decay_card_modeI = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta(1405) PHSP;
  Enddecay
  Decay eta(1405)
  1.000 f_0(980) pi0 PHSP;
  Enddecay
  Decay f_0(980)
  1.000 pi+ pi- PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode II signal: J/psi -> gamma eta(1405), eta(1405) -> f0(980) pi0, f0(980) -> pi0 pi0
decay_card_modeII = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta(1405) PHSP;
  Enddecay
  Decay eta(1405)
  1.000 f_0(980) pi0 PHSP;
  Enddecay
  Decay f_0(980)
  1.000 pi0 pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# eta' control mode (charged): J/psi -> gamma eta', eta' -> pi+ pi- pi0
decay_card_etap_charged = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta' PHSP;
  Enddecay
  Decay eta'
  1.000 pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# eta' control mode (neutral): J/psi -> gamma eta', eta' -> 3pi0
decay_card_etap_neutral = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta' PHSP;
  Enddecay
  Decay eta'
  1.000 pi0 pi0 pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# eta' peaking background 1: J/psi -> gamma eta', eta' -> gamma rho0, rho0 -> pi+ pi-
decay_card_etap_rho = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta' PHSP;
  Enddecay
  Decay eta'
  1.000 gamma rho0 PHSP;
  Enddecay
  Decay rho0
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# eta' peaking background 2: J/psi -> gamma eta', eta' -> gamma omega, omega -> pi+ pi- pi0
decay_card_etap_omega = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta' PHSP;
  Enddecay
  Decay eta'
  1.000 gamma omega PHSP;
  Enddecay
  Decay omega
  1.000 pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# ---------------- Exclusive MC ----------------
# 200k events each for the eta(1405) charged/neutral and eta' charged/neutral modes
exMC_eta1405_charged = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_eta1405_charged"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_eta1405_neutral = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_eta1405_neutral"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

exMC_etap_charged = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_etap_charged"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_etap_charged
  config.cross_section   = :default
end

exMC_etap_neutral = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_etap_neutral"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_etap_neutral
  config.cross_section   = :default
end

# 100k events each for the two eta' peaking backgrounds
exMC_etap_rho = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_etap_rho"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_etap_rho
  config.cross_section   = :default
end

exMC_etap_omega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_etap_omega"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_etap_omega
  config.cross_section   = :default
end

### ================= Event selection (BOSS) =================

# ------------------------------------------------------------
# Mode I  (gamma pi+ pi- pi0)  ->  eta(1405) charged & eta' charged
# ------------------------------------------------------------
alg_modeI = Algorithm.new("Eta1405ChargedMode")
alg_modeI.set_header(["Eta1405ChargedModeAlg/Eta1405ChargedMode.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
# competing-hypothesis chi2 values (2gamma / 4gamma) are stored for a ROOT-level veto
alg_modeI.note(:background_veto,
  "the nominal 3gamma 4C fit is compared against the 2gamma (eta'->gamma rho0) and " \
  "4gamma (eta'->gamma omega) hypotheses; the veto chi2_4c(3gamma) < chi2_4c(2gamma) " \
  "and chi2_4c(3gamma) < chi2_4c(4gamma) is applied in the ROOT analysis")

sel_modeI = Selection.new
  .select_track {                       # exactly one pi+ and one pi-, net charge zero
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        20.0                      # |Vz| < 20 cm
    Vr        2.0                       # Vr < 2 cm
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {                      # at least three photons
    tdc_emc_start     0                 # TDC 0-14
    tdc_emc_end       14
    angle_to_track    10.0              # > 10 deg from any charged track
    energyThreshold_b 0.025             # 25 MeV (barrel)
    energyThreshold_e 0.050             # 50 MeV (endcap)
    nGam              ">=3"
  }
  .pid(method: :probability) {          # pion PID against kaons and protons
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip     "==1"
    npim     "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # build pi0 from a gamma-gamma pair
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).within(0.1200, 0.1500)  # |M(gg) - m_pi0| < 15 MeV
    chi2_cut 25
    npi0     ">=1"
  }
  .invariant_mass_of(:gamma, :pi0).out_of(0.7327, 0.8327)     # veto |M(gamma pi0) - m_omega| < 50 MeV
  .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {      # nominal 4C fit: gamma gamma gamma pi+ pi-
    nominal
    constrain_four_momentum
    chi2_cut 30
  }
  .kinematic_fit([:gamma, :gamma, :pip, :pim]) {              # 2gamma hypothesis (eta' -> gamma rho0)
    constrain_four_momentum
  }
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim]) {  # 4gamma hypothesis (eta' -> gamma omega)
    constrain_four_momentum
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
root_files_modeI = alg_modeI.execute_on([
  jpsi_data, jpsi_incMC,
  exMC_eta1405_charged, exMC_etap_charged,
  exMC_etap_rho, exMC_etap_omega
])

# ------------------------------------------------------------
# Mode II  (gamma 3 pi0)  ->  eta(1405) neutral & eta' neutral
# ------------------------------------------------------------
alg_modeII = Algorithm.new("Eta1405NeutralMode")
alg_modeII.set_header(["Eta1405NeutralModeAlg/Eta1405NeutralMode.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
# |cos(theta_decay)| < 0.95 of the daughter photon in the pi0 rest frame is not expressible
alg_modeII.note(:efficiency_curve,
  "for each pi0 candidate the daughter photon is required to satisfy " \
  "|cos(theta_decay)| < 0.95 in the pi0 rest frame (combinatorial pi0 suppression); " \
  "this variable is not expressible in the DSL and must be applied in the ROOT analysis")

sel_modeII = Selection.new
  .select_track {                       # no charged tracks, net charge zero
    nChrp "==0"
    nChrn "==0"
    nNet  "==0"
  }
  .select_photon {                      # 7 or 8 photons, same energy/angle/TDC requirements
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=7"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # at least three pi0 candidates
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0     ">=3"
  }
  .kinematic_fit([:gamma, :pi0, :pi0, :pi0]) {   # nominal 4C fit: gamma pi0 pi0 pi0
    nominal
    constrain_four_momentum
    chi2_cut 60
    invariant_mass_of(:gamma, :pi0).out_of(0.7327, 0.8327)   # veto |M(gamma pi0) - m_omega| < 50 MeV
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
root_files_modeII = alg_modeII.execute_on([
  jpsi_data, jpsi_incMC,
  exMC_eta1405_neutral, exMC_etap_neutral
])