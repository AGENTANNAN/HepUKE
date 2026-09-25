# ============================================================================
# Search for Lambda_c+ -> Sigma+ eta, Sigma+ eta', Sigma+ pi0, Sigma+ omega
# at sqrt(s) = 4.600 GeV
# ============================================================================

### ---------------------------- Dataset preparation --------------------------- ###
data_4600  = DatasetManager.real_data.find("703_4600")     # 4.600 GeV real data (BOSS 703)
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")  # corresponding inclusive MC

### Decay cards (EvtGen syntax) — one per signal mode ###
# Mode (a): Lambda_c+ -> Sigma+ eta, Sigma+ -> p pi0, eta -> gamma gamma, pi0 -> gamma gamma
decay_card_a = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay
  Decay Lambda_c+
  1.0000 Sigma+ eta PHSP;
  Enddecay
  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  Decay anti-Lambda_c-
  1.0000 anti-p- pi0 PHSP;
  Enddecay
  End
DECAYCARD

# Mode (b): Lambda_c+ -> Sigma+ eta', Sigma+ -> p pi0, eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_b = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay
  Decay Lambda_c+
  1.0000 Sigma+ eta' PHSP;
  Enddecay
  Decay Sigma+
  1.0000 p+ pi0 PHSP;
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
  Decay anti-Lambda_c-
  1.0000 anti-p- pi0 PHSP;
  Enddecay
  End
DECAYCARD

# Mode (c): Lambda_c+ -> Sigma+ pi0, Sigma+ -> p pi0, pi0 -> gamma gamma
decay_card_c = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay
  Decay Lambda_c+
  1.0000 Sigma+ pi0 PHSP;
  Enddecay
  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  Decay anti-Lambda_c-
  1.0000 anti-p- pi0 PHSP;
  Enddecay
  End
DECAYCARD

# Mode (d): Lambda_c+ -> Sigma+ omega, Sigma+ -> p pi0, omega -> pi+ pi- pi0 (Dalitz)
decay_card_d = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay
  Decay Lambda_c+
  1.0000 Sigma+ omega PHSP;
  Enddecay
  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay
  Decay omega
  1.0000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  Decay anti-Lambda_c-
  1.0000 anti-p- pi0 PHSP;
  Enddecay
  End
DECAYCARD

### Exclusive MC — four 500k-event samples, one per decay mode ###
exMC_a = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4600_lc_to_sigma_eta"
  config.related_dataset = data_4600
  config.events          = 500000
  config.decay_card      = decay_card_a
  config.cross_section   = :default
end

exMC_b = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4600_lc_to_sigma_etap"
  config.related_dataset = data_4600
  config.events          = 500000
  config.decay_card      = decay_card_b
  config.cross_section   = :default
end

exMC_c = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4600_lc_to_sigma_pi0"
  config.related_dataset = data_4600
  config.events          = 500000
  config.decay_card      = decay_card_c
  config.cross_section   = :default
end

exMC_d = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4600_lc_to_sigma_omega"
  config.related_dataset = data_4600
  config.events          = 500000
  config.decay_card      = decay_card_d
  config.cross_section   = :default
end

### ------------------------------ Event selection ---------------------------- ###
# Common cuts: |cos(theta)|<0.93, |Vz|<10 cm, Vr<1 cm;
# photons: E>25 MeV (barrel) / 50 MeV (endcap), TDC 0-14 (700 ns),
# opening angle >10 deg to any charged track, >=4 photons;
# PID (probability) threshold 0.001, protons vs (K,pi), pions vs (K,p).

# ---------------------------------------------------------------------------
# Mode (a): Lambda_c+ -> Sigma+ eta  (Sigma+ -> p pi0, eta -> gamma gamma)
# ---------------------------------------------------------------------------
alg_a = Algorithm.new("LcToSigmaEta")
alg_a.set_header(["LcToSigmaEtaAlg/LcToSigmaEta.h"])
     .set_constant({ "ECMS" => [:double, 4.600] })
     .note(:delta_q_window,
           "DeltaQ = M(Sigma+ eta) - m(Lambda_c+) window [-0.032, 0.022] GeV applied at ROOT level after the 4C fit")

sel_a = Selection.new
sel_a.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     ">=1"   # at least one p (from Sigma+ -> p pi0)
        nChrn     ">=1"   # at least one p-bar recoiling against the Lambda_c+
        nNet      "==0"
      }
     .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=4"   # eta -> gamma gamma + pi0 -> gamma gamma
      }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :pion,   against: [:kaon, :proton]
        nprp ">=1"
        nprm ">=1"
      }
     .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 200
        npi0 ">=1"
      }
     .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 200
        neta ">=1"
      }
     .kinematic_fit([:prp, :pi0, :eta]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg_a.with_decay_card(decay_card_a).apply(sel_a)
alg_a.execute_on([data_4600, incMC_4600, exMC_a])

# ---------------------------------------------------------------------------
# Mode (b): Lambda_c+ -> Sigma+ eta'  (Sigma+ -> p pi0, eta' -> pi+ pi- eta)
# ---------------------------------------------------------------------------
alg_b = Algorithm.new("LcToSigmaEtap")
alg_b.set_header(["LcToSigmaEtapAlg/LcToSigmaEtap.h"])
     .set_constant({ "ECMS" => [:double, 4.600] })
     .note(:delta_q_window,
           "DeltaQ = M(Sigma+ eta') - m(Lambda_c+) window [-0.030, 0.020] GeV applied at ROOT level after the 4C fit")

sel_b = Selection.new
sel_b.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     ">=2"   # p (from Sigma+) and pi+ (from eta')
        nChrn     ">=1"   # pi- (from eta')
        nNet      "==1"
      }
     .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=4"   # eta -> gamma gamma + pi0 -> gamma gamma
      }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :pion,   against: [:kaon, :proton]
        nprp ">=1"
        npip ">=1"
        npim ">=1"
      }
     .remove([:prp <= :pip, :prm <= :pim])   # remove proton/pion candidate overlap
     .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 200
        npi0 ">=1"
      }
     .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 200
        neta ">=1"
      }
     .kinematic_fit([:prp, :pi0, :pip, :pim, :eta]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg_b.with_decay_card(decay_card_b).apply(sel_b)
alg_b.execute_on([data_4600, incMC_4600, exMC_b])

# ---------------------------------------------------------------------------
# Mode (c): Lambda_c+ -> Sigma+ pi0  (Sigma+ -> p pi0)
# ---------------------------------------------------------------------------
alg_c = Algorithm.new("LcToSigmaPi0")
alg_c.set_header(["LcToSigmaPi0Alg/LcToSigmaPi0.h"])
     .set_constant({ "ECMS" => [:double, 4.600] })
     .note(:delta_q_window,
           "DeltaQ = M(Sigma+ pi0) - m(Lambda_c+) window [-0.050, 0.030] GeV applied at ROOT level after the 4C fit")

sel_c = Selection.new
sel_c.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     ">=1"   # p (from Sigma+ -> p pi0)
        nChrn     ">=1"   # p-bar recoiling against the Lambda_c+
        nNet      "==0"
      }
     .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=4"   # two pi0 -> 4 photons
      }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :pion,   against: [:kaon, :proton]
        nprp ">=1"
        nprm ">=1"
      }
     .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 200
        npi0 ">=2"
      }
     .kinematic_fit([:prp, :pi0, :pi0]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
        invariant_mass_of(:pi0, :pi0).out_of(0.4776, 0.5176)  # veto K_S0 -> pi0 pi0 (|M - m_K_S0| > 0.02)
      }
alg_c.with_decay_card(decay_card_c).apply(sel_c)
alg_c.execute_on([data_4600, incMC_4600, exMC_c])

# ---------------------------------------------------------------------------
# Mode (d): Lambda_c+ -> Sigma+ omega  (Sigma+ -> p pi0, omega -> pi+ pi- pi0)
# ---------------------------------------------------------------------------
alg_d = Algorithm.new("LcToSigmaOmega")
alg_d.set_header(["LcToSigmaOmegaAlg/LcToSigmaOmega.h"])
     .set_constant({ "ECMS" => [:double, 4.600] })
     .note(:delta_q_window,
           "DeltaQ = M(Sigma+ omega) - m(Lambda_c+) window [-0.030, 0.020] GeV applied at ROOT level after the 4C fit")

sel_d = Selection.new
sel_d.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     ">=2"   # p (from Sigma+) and pi+ (from omega)
        nChrn     ">=1"   # pi- (from omega)
        nNet      "==1"
      }
     .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=4"   # pi0 (Sigma+) -> gamma gamma + pi0 (omega) -> gamma gamma
      }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :pion,   against: [:kaon, :proton]
        nprp ">=1"
        npip ">=1"
        npim ">=1"
      }
     .remove([:prp <= :pip, :prm <= :pim])   # remove proton/pion candidate overlap
     .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 200
        npi0 ">=2"
      }
     .kinematic_fit([:prp, :pi0, :pip, :pim, :pi0]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg_d.with_decay_card(decay_card_d).apply(sel_d)
alg_d.execute_on([data_4600, incMC_4600, exMC_d])