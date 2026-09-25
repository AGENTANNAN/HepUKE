# =====================================================================
#  a0(980)-f0(980) mixing study at BESIII  (BOSS / event-selection part)
#   Channel I : J/psi -> phi f0(980) -> phi a0(980) -> phi eta pi0
#                 phi -> K+K-,  eta -> gamma gamma,  pi0 -> gamma gamma
#   Channel II: psi'  -> gamma chi_c1 -> gamma pi0 a0(980) -> gamma pi0 f0(980)
#                 -> gamma pi0 pi+ pi-,  pi0 -> gamma gamma
#   The narrow mixing scalar S (m = 991.3 MeV, Gamma = 8 MeV) is the
#   intermediate that encodes the a0/f0 mixing:  S -> eta pi0 (J/psi mode)
#   and  S -> pi+ pi- (psi' mode).
# =====================================================================

### ------------------------------ Datasets ------------------------------ ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi   (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # J/psi   inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) (3.686 GeV) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # psi(2S) inclusive MC

### ------------------------------ Decay cards --------------------------- ###
# Channel I : J/psi -> phi S ; S -> eta pi0 ; phi -> K+ K- ; eta -> gg ; pi0 -> gg
decay_card_channelI = <<~DECAYCARD
  Particle S 0.9913 0.008 0.0 0.0 0 0 0 0 0 0 0 0

  Decay J/psi
  1.000 phi S PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  Decay S
  1.000 eta pi0 PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Channel II : psi(2S) -> gamma chi_c1 ; chi_c1 -> gamma pi0 S ; S -> pi+ pi- ; pi0 -> gg
decay_card_channelII = <<~DECAYCARD
  Particle S 0.9913 0.008 0.0 0.0 0 0 0 0 0 0 0 0

  Decay psi(2S)
  1.000 gamma chi_c1 P2GC1;
  Enddecay

  Decay chi_c1
  1.000 gamma pi0 S PHSP;
  Enddecay

  Decay S
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### ---------------------- Exclusive MC (200k events) -------------------- ###
exMC_channelI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_phi_a0f0_etapi0"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_channelI
  config.cross_section   = :default
end

exMC_channelII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_a0f0_gammapi0pipi"
  config.related_dataset = psip_data
  config.events          = 200_000
  config.decay_card      = decay_card_channelII
  config.cross_section   = :default
end

# Persist MC configurations for reproducibility
exMC_channelI.save_to_config(format: :yaml, file_path: 'temp_for_test')
exMC_channelII.save_to_config(format: :yaml, file_path: 'temp_for_test')

### ============ Channel I : J/psi -> phi eta pi0 (K+K- gg gg) ========== ###
alg_name_I = "A0F0MixJpsi"
alg_I = Algorithm.new(alg_name_I)
alg_I.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
     .set_constant({"ECMS" => [:double, 3.097]})          # sqrt(s) = 3.097 GeV
     .set_alias({"std::vector<double>" => "Vdouble"})

sel_I = Selection.new
sel_I.select_track {
       cos_theta 0.93      # |cos(theta)| < 0.93
       Vz        20.0      # |Vz| < 20 cm
       Vr        2.0       # Vr < 2 cm
       nChrp     "==1"     # exactly one positive track
       nChrn     "==1"     # exactly one negative track
       nNet      "==0"     # net charge = 0
     }
     .select_photon {
       tdc_emc_start     0      # TDC start = 0
       tdc_emc_end       14     # TDC end   = 14
       angle_to_track    10.0   # > 10 deg from the nearest charged track
       energyThreshold_b 0.025  # E > 25 MeV in the barrel
       energyThreshold_e 0.050  # E > 50 MeV in the endcap
       nGam              ">=4"  # at least four photons (four from eta, pi0)
     }
     .pid(method: :probability) {
       prob_cut 0.001                              # PID probability > 0.001
       identify :kaon, against: [:pion, :proton]   # one K+ and one K-
       nkp "==1"
       nkm "==1"
     }
     # eta -> gamma gamma (1C Kalman mass fit)
     .kalman_kinematic_fit([:gamma, :gamma]) {
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
       chi2_cut 25
       neta ">=1"
     }
     # pi0 -> gamma gamma (1C Kalman mass fit)
     .kalman_kinematic_fit([:gamma, :gamma]) {
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
       chi2_cut 25
       npi0 ">=1"
     }
     # nominal fit: 4C four-momentum + eta/pi0 mass constraints (via Kalman) = 6C
     .kinematic_fit([:kp, :km, :eta, :pi0]) {
       nominal
       constrain_four_momentum
       invariant_mass_of(:kp, :km).within(1.005, 1.035)  # phi window |m(K+K-)-1.02| < 0.015 GeV
       chi2_cut 200   # loose BOSS-level cut (paper applies chi2 < 60 in ROOT)
     }
     # competing hypothesis J/psi -> K+ K- eta eta ; chi2 stored for the probability veto
     .kinematic_fit([:kp, :km, :eta, :eta]) {
       constrain_four_momentum
     }
     # competing hypothesis J/psi -> K+ K- pi0 pi0 ; chi2 stored for the probability veto
     .kinematic_fit([:kp, :km, :pi0, :pi0]) {
       constrain_four_momentum
     }

alg_I.note(:background_veto,
           "ROOT-level probability veto: reject events favoured by J/psi -> K+K-eta eta " \
           "or J/psi -> K+K-pi0 pi0 using the stored competing chi2 values")

alg_I.with_decay_card(decay_card_channelI).apply(sel_I)

### ============ Channel II : psi' -> gamma pi0 pi+ pi- ================= ###
alg_name_II = "A0F0MixPsip"
alg_II = Algorithm.new(alg_name_II)
alg_II.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})         # sqrt(s) = 3.686 GeV
      .set_alias({"std::vector<double>" => "Vdouble"})

sel_II = Selection.new
sel_II.select_track {
        cos_theta 0.93      # |cos(theta)| < 0.93
        Vz        20.0      # |Vz| < 20 cm
        Vr        2.0       # Vr < 2 cm
        nChrp     "==1"     # exactly one positive track
        nChrn     "==1"     # exactly one negative track
        nNet      "==0"     # net charge = 0
      }
      .select_photon {
        tdc_emc_start     0      # TDC start = 0
        tdc_emc_end       14     # TDC end   = 14
        angle_to_track    10.0   # > 10 deg from the nearest charged track
        energyThreshold_b 0.025  # E > 25 MeV in the barrel
        energyThreshold_e 0.050  # E > 50 MeV in the endcap
        nGam              ">=3"  # at least three photons (two from pi0)
      }
      .pid(method: :probability) {
        prob_cut 0.001                              # PID probability > 0.001
        identify :pion, against: [:kaon, :proton]   # one pi+ and one pi-
        npip "==1"
        npim "==1"
      }
      # pi0 -> gamma gamma (1C Kalman mass fit)
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      }
      # nominal fit: 4C four-momentum + pi0 mass constraint (via Kalman) = 5C
      .kinematic_fit([:gamma, :pi0, :pip, :pim]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:pi0, :pip, :pim).within(3.49, 3.54)  # chi_c1 window 3.49-3.54 GeV
        chi2_cut 200   # loose BOSS-level cut (paper applies chi2 < 60 in ROOT)
      }
      # competing hypothesis psi' -> pi0 pi+ pi- ; chi2 stored for the probability veto
      .kinematic_fit([:pi0, :pip, :pim]) {
        constrain_four_momentum
      }
      # competing hypothesis psi' -> pi0 pi0 pi+ pi- ; chi2 stored for the probability veto
      .kinematic_fit([:pi0, :pi0, :pip, :pim]) {
        constrain_four_momentum
      }

alg_II.note(:background_veto,
            "off-line ROOT-level veto: reject gamma-gamma recoil J/psi " \
            "(|M_recoil - 3.097| < 0.06 GeV); probability veto using the stored " \
            "competing chi2 of psi' -> pi0 pi+ pi- and psi' -> pi0 pi0 pi+ pi-")

alg_II.with_decay_card(decay_card_channelII).apply(sel_II)

### ------------------------------- Execute ------------------------------ ###
root_files_I  = alg_I.execute_on([jpsi_data, jpsi_incMC, exMC_channelI])
root_files_II = alg_II.execute_on([psip_data, psip_incMC, exMC_channelII])