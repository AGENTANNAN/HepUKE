# Paper: 1811.08028v4 — Evidence for Lambda_c+ → Sigma+ eta and Sigma+ eta'
# BESIII, √s = 4.6 GeV, 567 pb⁻¹
# Ordinary analysis: 4 decay modes with Lambda_c+ reconstruction

### Dataset preparation ###
data_4600 = DatasetManager.real_data.find("703_4600")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

# ==============================================================
# Mode (a): Lambda_c+ → Sigma+ eta  (Σ+ → p π0, η → γγ)
# ==============================================================
decay_card_a = <<~DECAYCARD
  Decay psi(4260)
  1.000 anti-Lambda_c- Lambda_c+ PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 Sigma+ eta PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_a = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Lc_Sigma_eta"
  config.related_dataset = data_4600
  config.events = 500_000
  config.decay_card = decay_card_a
  config.cross_section = :default
end

alg_a = Algorithm.new("LcSigmaEta")
alg_a.set_header(["LcSigmaEtaAlg/LcSigmaEta.h"])
  .set_constant({ "ECMS" => [:double, 4.600] })
  .note(:background_veto, "anti-proton recoiling against Lambda_c+ required to suppress combinatorial background; anti-proton identified via Lp > Lpi and Lp > LK")
  .note(:tag_mode_unavailable, "DeltaQ selection window [-0.032, 0.022] GeV applied at BOSS level; M_BC fit and signal yield extraction performed at ROOT level")

sel_a = Selection.new
sel_a.select_track {
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       ">=1"
  nChrn       ">=1"
  nNet        "==0"
}
.select_photon {
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=4"
}
.pid(method: :probability) {
  prob_cut   0.001
  identify :proton, against: [:kaon, :pion]
  nprp       ">=1"
  nprm       ">=1"
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0  ">=1"
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 200
  neta  ">=1"
}
.kinematic_fit([:prp, :pi0, :eta]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_a.with_decay_card(decay_card_a).apply(sel_a)
alg_a.execute_on([data_4600, incMC_4600, exMC_a])

# ==============================================================
# Mode (b): Lambda_c+ → Sigma+ eta'  (Σ+ → p π0, η' → π+π-η, η → γγ)
# ==============================================================
decay_card_b = <<~DECAYCARD
  Decay psi(4260)
  1.000 anti-Lambda_c- Lambda_c+ PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 Sigma+ eta' PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay eta'
  1.000 pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_b = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Lc_Sigma_etap"
  config.related_dataset = data_4600
  config.events = 500_000
  config.decay_card = decay_card_b
  config.cross_section = :default
end

alg_b = Algorithm.new("LcSigmaEtap")
alg_b.set_header(["LcSigmaEtapAlg/LcSigmaEtap.h"])
  .set_constant({ "ECMS" => [:double, 4.600] })
  .note(:tag_mode_unavailable, "DeltaQ selection window [-0.030, 0.020] GeV applied at BOSS level; M_BC fit and signal yield extraction performed at ROOT level")

sel_b = Selection.new
sel_b.select_track {
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       ">=2"
  nChrn       ">=1"
  nNet        "==1"
}
.select_photon {
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=4"
}
.pid(method: :probability) {
  prob_cut   0.001
  identify :proton, against: [:kaon, :pion]
  nprp       ">=1"
  identify :pion, against: [:kaon, :proton]
  npip       ">=1"
  npim       ">=1"
}
.remove([:prp <= :chrgp])
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0  ">=1"
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 200
  neta  ">=1"
}
.kinematic_fit([:prp, :pi0, :pip, :pim, :eta]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_b.with_decay_card(decay_card_b).apply(sel_b)
alg_b.execute_on([data_4600, incMC_4600, exMC_b])

# ==============================================================
# Mode (c): Lambda_c+ → Sigma+ pi0  (Σ+ → p π0, reference mode)
# ==============================================================
decay_card_c = <<~DECAYCARD
  Decay psi(4260)
  1.000 anti-Lambda_c- Lambda_c+ PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 Sigma+ pi0 PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_c = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Lc_Sigma_pi0"
  config.related_dataset = data_4600
  config.events = 500_000
  config.decay_card = decay_card_c
  config.cross_section = :default
end

alg_c = Algorithm.new("LcSigmaPi0")
alg_c.set_header(["LcSigmaPi0Alg/LcSigmaPi0.h"])
  .set_constant({ "ECMS" => [:double, 4.600] })
  .note(:background_veto, "anti-proton recoiling against Lambda_c+ required to suppress combinatorial background; also K_S0 veto: |M(π0π0) - m_K_S0| > 0.02 GeV/c²")
  .note(:tag_mode_unavailable, "DeltaQ selection window [-0.050, 0.030] GeV applied at BOSS level; M_BC fit and signal yield extraction performed at ROOT level")

sel_c = Selection.new
sel_c.select_track {
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       ">=1"
  nChrn       ">=1"
  nNet        "==0"
}
.select_photon {
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=4"
}
.pid(method: :probability) {
  prob_cut   0.001
  identify :proton, against: [:kaon, :pion]
  nprp       ">=1"
  nprm       ">=1"
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0  ">=2"
}
.kinematic_fit([:prp, :pi0, :pi0]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_c.with_decay_card(decay_card_c).apply(sel_c)
alg_c.execute_on([data_4600, incMC_4600, exMC_c])

# ==============================================================
# Mode (d): Lambda_c+ → Sigma+ omega  (Σ+ → p π0, ω → π+π-π0)
# ==============================================================
decay_card_d = <<~DECAYCARD
  Decay psi(4260)
  1.000 anti-Lambda_c- Lambda_c+ PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 Sigma+ omega PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_d = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Lc_Sigma_omega"
  config.related_dataset = data_4600
  config.events = 500_000
  config.decay_card = decay_card_d
  config.cross_section = :default
end

alg_d = Algorithm.new("LcSigmaOmega")
alg_d.set_header(["LcSigmaOmegaAlg/LcSigmaOmega.h"])
  .set_constant({ "ECMS" => [:double, 4.600] })
  .note(:tag_mode_unavailable, "DeltaQ selection window [-0.030, 0.020] GeV applied at BOSS level; M_BC fit and signal yield extraction performed at ROOT level")

sel_d = Selection.new
sel_d.select_track {
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       ">=2"
  nChrn       ">=1"
  nNet        "==1"
}
.select_photon {
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=4"
}
.pid(method: :probability) {
  prob_cut   0.001
  identify :proton, against: [:kaon, :pion]
  nprp       ">=1"
  identify :pion, against: [:kaon, :proton]
  npip       ">=1"
  npim       ">=1"
}
.remove([:prp <= :chrgp])
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 200
  npi0  ">=2"
}
.kinematic_fit([:prp, :pi0, :pip, :pim, :pi0]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_d.with_decay_card(decay_card_d).apply(sel_d)
alg_d.execute_on([data_4600, incMC_4600, exMC_d])