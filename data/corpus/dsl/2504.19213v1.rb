# DSL: D^0/D+ → K 5π hadronic decays (DT method at ψ(3770))
# Paper: 2504.19213v1 — TagAnalysis surface
# Three signal decays: 1) D0→K-3π+2π-, 2) D0→K-2π+π-2π0, 3) D+→K-3π+π-π0
# Tag modes from authoritative vocabulary only (tag_modes_vocabulary.md)

# ============================================================
# Datasets
# ============================================================
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================
# Decay cards — exclusive MC with full DT chain
# ============================================================

# Signal 1: D0 → K- 3π+ 2π-  (tag anti-D0 → K+ π-)
decay_card_d0_k3pi2pi = <<~DECAYCARD
  Decay psi(3770)
  1 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1 K- pi+ pi+ pi+ pi- pi- PHSP;
  Enddecay
  Decay anti-D0
  1 K+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Signal 2: D0 → K- 2π+ π- 2π0  (tag anti-D0 → K+ π-)
decay_card_d0_k2pi1pi2pi0 = <<~DECAYCARD
  Decay psi(3770)
  1 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1 K- pi+ pi+ pi- pi0 pi0 PHSP;
  Enddecay
  Decay anti-D0
  1 K+ pi- PHSP;
  Enddecay
  Decay pi0
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Signal 3: D+ → K- 3π+ π- π0  (tag D- → K+ π- π-)
decay_card_dp_k3pi1pi_pi0 = <<~DECAYCARD
  Decay psi(3770)
  1 D+ D- PHSP;
  Enddecay
  Decay D+
  1 K- pi+ pi+ pi+ pi- pi0 PHSP;
  Enddecay
  Decay D-
  1 K+ pi- pi- PHSP;
  Enddecay
  Decay pi0
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# ============================================================
# Signal MC
# ============================================================
sig_d0_k3pi2pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_d0_k3pi2pi"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_d0_k3pi2pi
  config.cross_section   = :default
end

sig_d0_k2pi1pi2pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_d0_k2pi1pi2pi0"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_d0_k2pi1pi2pi0
  config.cross_section   = :default
end

sig_dp_k3pi1pi_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_dp_k3pi1pi_pi0"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_dp_k3pi1pi_pi0
  config.cross_section   = :default
end

# ============================================================
# Common tag modes (used by all three algorithms)
# D0 ST: D0→K+π-, D0→K+π-π0, D0→K+π-π-π+
# D+ ST: D-→K+π-π-, D-→K_S0π-, D-→K+π-π-π0, D-→K_S0π-π0, D-→K_S0π+π-π-, D-→K+K-π-
# ============================================================

# ============================================================
# Algorithm 1: D0 → K- 3π+ 2π-  (all charged, 6 tracks)
# ============================================================
alg1 = TagAnalysis.new("D0ToK3Pi2Pi_DT")
alg1.set_header(["D0ToK3Pi2PiAlg/D0ToK3Pi2Pi.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })

alg1.tag_side(:D0) do
  modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
end

alg1.tag_side(:Dplus) do
  modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
end

alg1.signal_side do
  charged km: 1, pip: 3, pim: 2
  photons 0
end

alg1.fit do
  constrain_four_momentum
  chi2_cut 200
end

alg1.with_decay_card(decay_card_d0_k3pi2pi).apply
alg1.note(:signal_mode, "D0 → K- 3π+ 2π- (all charged). χ_c²(DT) < 200 in BOSS; optimal cut tuned in ROOT.")
alg1.note(:dt_method, "Double tag at ψ(3770): tag D (D0 or D+) via ST hadronic modes, signal D reconstructed from remaining tracks. ΔE and M_BC stored for both tag and signal sides; window cuts applied in ROOT.")
alg1.note(:ks_veto, "K_S0 veto on signal side: events with K_S0 candidates (π+π- vertex, |M_ππ-m_K_S0|<12 MeV) in the signal side are rejected in ROOT.")
alg1.execute_on([data_3773, incMC_3773, sig_d0_k3pi2pi])

# ============================================================
# Algorithm 2: D0 → K- 2π+ π- 2π0  (4 charged + 4γ)
# ============================================================
alg2 = TagAnalysis.new("D0ToK2Pi1Pi2Pi0_DT")
alg2.set_header(["D0ToK2Pi1Pi2Pi0Alg/D0ToK2Pi1Pi2Pi0.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })

alg2.tag_side(:D0) do
  modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
end

alg2.tag_side(:Dplus) do
  modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
end

alg2.signal_side do
  charged km: 1, pip: 2, pim: 1
  photons 4
end

alg2.fit do
  constrain_four_momentum
  chi2_cut 200
end

alg2.with_decay_card(decay_card_d0_k2pi1pi2pi0).apply
alg2.note(:signal_mode, "D0 → K- 2π+ π- 2π0. π0 → γγ. χ_c²(DT) < 200 in BOSS; optimal cut tuned in ROOT.")
alg2.note(:dt_method, "Double tag at ψ(3770): same tag sides as Algorithm 1. ΔE and M_BC stored for both tag and signal sides; window cuts applied in ROOT.")
alg2.note(:ks_veto, "K_S0 veto on signal side applied in ROOT (as in Algorithm 1).")
alg2.execute_on([data_3773, incMC_3773, sig_d0_k2pi1pi2pi0])

# ============================================================
# Algorithm 3: D+ → K- 3π+ π- π0  (5 charged + 2γ)
# ============================================================
alg3 = TagAnalysis.new("DpToK3Pi1PiPi0_DT")
alg3.set_header(["DpToK3Pi1PiPi0Alg/DpToK3Pi1PiPi0.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })

alg3.tag_side(:D0) do
  modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
end

alg3.tag_side(:Dplus) do
  modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
end

alg3.signal_side do
  charged km: 1, pip: 3, pim: 1
  photons 2
end

alg3.fit do
  constrain_four_momentum
  chi2_cut 200
end

alg3.with_decay_card(decay_card_dp_k3pi1pi_pi0).apply
alg3.note(:signal_mode, "D+ → K- 3π+ π- π0. π0 → γγ. χ_c²(DT) < 200 in BOSS; optimal cut tuned in ROOT.")
alg3.note(:dt_method, "Double tag at ψ(3770): same tag sides as Algorithm 1. ΔE and M_BC stored for both tag and signal sides; window cuts applied in ROOT.")
alg3.note(:ks_veto, "K_S0 veto on signal side applied in ROOT (as in Algorithm 1).")
alg3.execute_on([data_3773, incMC_3773, sig_dp_k3pi1pi_pi0])