# 1907.11258v2: D meson decays to φ + pseudoscalar at ψ(3770)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ====== Mode I: D+ → φπ+, φ → K+K- ======
decay_card_Dp_phipi = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D- PHSP;
  Enddecay
  Decay D+
  1.0000 K+ K- pi+ PHSP;
  Enddecay
  Decay D-
  1.0000 anti-K0 pi- PHSP;
  Enddecay
  End
DECAYCARD

exMC_Dp_phipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Dp_phipi"
  config.related_dataset = psi3770_data
  config.events = 100_000
  config.decay_card = decay_card_Dp_phipi
  config.cross_section = :default
end

alg_I = Algorithm.new("DpPhiPi")
alg_I.set_header(["DpPhiPiAlg/DpPhiPi.h"])
      .set_constant({"ECMS" => [:double, 3.773]})

sel_I = Selection.new
sel_I.select_track {
  cos_theta 0.93
  Vz 100.0
  Vr 10.0
  nChrp ">=2"
  nChrn ">=1"
  nNet "==1"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  identify :pion, against: [:kaon, :proton]
  nkp "==1"
  nkm "==1"
  npip "==1"
}
.kinematic_fit([:kp, :km, :pip]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_I.with_decay_card(decay_card_Dp_phipi).apply(sel_I)
alg_I.note(:deltaE_cut, "ΔE ∈ [-0.020, 0.019] GeV (3σ) applied in ROOT")
     .note(:multiple_candidates, "If >1 D candidate, keep the one with smallest |ΔE|")
     .note(:twoD_fit, "2D unbinned ML fit to M_BC vs M_KK; signal from MC shape convolved with Gaussian; combinatorial bg: ARGUS functions")

alg_I.execute_on([psi3770_data, psi3770_incMC, exMC_Dp_phipi])

# ====== Mode II: D+ → φK+, φ → K+K- ======
decay_card_Dp_phiK = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D- PHSP;
  Enddecay
  Decay D+
  1.0000 K+ K- K+ PHSP;
  Enddecay
  Decay D-
  1.0000 anti-K0 pi- PHSP;
  Enddecay
  End
DECAYCARD

exMC_Dp_phiK = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Dp_phiK"
  config.related_dataset = psi3770_data
  config.events = 100_000
  config.decay_card = decay_card_Dp_phiK
  config.cross_section = :default
end

alg_II = Algorithm.new("DpPhiK")
alg_II.set_header(["DpPhiKAlg/DpPhiK.h"])
       .set_constant({"ECMS" => [:double, 3.773]})

sel_II = Selection.new
sel_II.select_track {
  cos_theta 0.93
  Vz 100.0
  Vr 10.0
  nChrp ">=2"
  nChrn ">=1"
  nNet "==1"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp "==2"
  nkm "==1"
}
.kinematic_fit([:kp, :kp, :km]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_II.with_decay_card(decay_card_Dp_phiK).apply(sel_II)
alg_II.note(:deltaE_cut, "ΔE ∈ [-0.019, 0.018] GeV (3σ) applied in ROOT")
      .note(:multiple_candidates, "If >1 D candidate, keep the one with smallest |ΔE|")
      .note(:twoD_fit, "2D unbinned ML fit to M_BC vs M_KK")

alg_II.execute_on([psi3770_data, psi3770_incMC, exMC_Dp_phiK])

# ====== Mode III: D0 → φπ0, φ → K+K-, π0 → γγ ======
decay_card_D0_phipi0 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0000 K+ K- pi0 PHSP;
  Enddecay
  Decay anti-D0
  1.0000 K+ pi- PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_D0_phipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_D0_phipi0"
  config.related_dataset = psi3770_data
  config.events = 100_000
  config.decay_card = decay_card_D0_phipi0
  config.cross_section = :default
end

alg_III = Algorithm.new("D0PhiPi0")
alg_III.set_header(["D0PhiPi0Alg/D0PhiPi0.h"])
        .set_constant({"ECMS" => [:double, 3.773]})

sel_III = Selection.new
sel_III.select_track {
  cos_theta 0.93
  Vz 100.0
  Vr 10.0
  nChrp ">=1"
  nChrn ">=1"
  nNet "==0"
}
.select_photon {
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  tdc_emc_start 0
  tdc_emc_end 700
  nGam ">=2"
}
.select_isolated_photon {
  angle_to_track 10.0
  nGam ">=2"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp "==1"
  nkm "==1"
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0 ">=1"
}
.kinematic_fit([:kp, :km, :pi0]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_III.with_decay_card(decay_card_D0_phipi0).apply(sel_III)
alg_III.note(:deltaE_cut, "ΔE ∈ [-0.077, 0.035] GeV (3σ, asymmetric) applied in ROOT")
       .note(:multiple_candidates, "If >1 D candidate, keep the one with smallest |ΔE|")
       .note(:twoD_fit, "2D unbinned ML fit to M_BC vs M_KK")

alg_III.execute_on([psi3770_data, psi3770_incMC, exMC_D0_phipi0])

# ====== Mode IV: D0 → φη, φ → K+K-, η → γγ ======
decay_card_D0_phieta = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0000 K+ K- eta PHSP;
  Enddecay
  Decay anti-D0
  1.0000 K+ pi- PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_D0_phieta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_D0_phieta"
  config.related_dataset = psi3770_data
  config.events = 100_000
  config.decay_card = decay_card_D0_phieta
  config.cross_section = :default
end

alg_IV = Algorithm.new("D0PhiEta")
alg_IV.set_header(["D0PhiEtaAlg/D0PhiEta.h"])
       .set_constant({"ECMS" => [:double, 3.773]})

sel_IV = Selection.new
sel_IV.select_track {
  cos_theta 0.93
  Vz 100.0
  Vr 10.0
  nChrp ">=1"
  nChrn ">=1"
  nNet "==0"
}
.select_photon {
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  tdc_emc_start 0
  tdc_emc_end 700
  nGam ">=2"
}
.select_isolated_photon {
  angle_to_track 10.0
  nGam ">=2"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp "==1"
  nkm "==1"
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 25
  neta ">=1"
}
.kinematic_fit([:kp, :km, :eta]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_IV.with_decay_card(decay_card_D0_phieta).apply(sel_IV)
alg_IV.note(:deltaE_cut, "ΔE ∈ [-0.040, 0.038] GeV (3σ) applied in ROOT")
       .note(:multiple_candidates, "If >1 D candidate, keep the one with smallest |ΔE|")
       .note(:twoD_fit, "2D unbinned ML fit to M_BC vs M_KK")
       .note(:quantum_coherence, "QC correction (ΔN_CP = y_CP * N_CP) for D0-D0bar coherence with 1.0% uncertainty")

alg_IV.execute_on([psi3770_data, psi3770_incMC, exMC_D0_phieta])