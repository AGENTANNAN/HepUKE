# DSL: 2503.05382v1 — Measurement of BF of D+ → K+K-π+π+π-, φπ+π+π-, K_S0K+π+π-π0, K_S0K+η, K_S0K+ω
# BESIII at ψ(3770), 20.3 fb-1, double-tag method
# TagAnalysis — DT pattern

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ═══════════════════════════════════════════════════════════════════════
# Mode 1: D+ → K+K-π+π+π-  (expressible as :DptoKKPiPiPi)
# ═══════════════════════════════════════════════════════════════════════
decay_card_m1 = <<~DECAYCARD
  Decay D+
  1.000 K+ K- pi+ pi+ pi- PHSP;
  Enddecay
DECAYCARD

exMC_m1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Dp_KKpipipi"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_m1
  config.cross_section   = :default
end

alg_m1 = TagAnalysis.new("DpToKKPiPiPi")
alg_m1.set_header(["DpToKKPiPiPiAlg/DpToKKPiPiPi.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })

# Tag side 1: D- tag (6 modes, charm -1)
alg_m1.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :deltaE, abs: 0.06    # ≈ ±3.5σ, paper Table 2
end

# Tag side 2: D+ signal (charm +1)
alg_m1.tag_side(:Dplus) do |t|
  t.modes :DptoKKPiPiPi
  t.charm 1
end

# No extra signal-side content (pure DT — all tracks consumed by the two tags)
alg_m1.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_m1.note(:tag_deltaE_windows, "Mode-dependent ΔE windows from Table 2: Kππ (-25,24), K_S0π (-25,26), Kπππ0 (-57,46), K_S0ππ0 (-62,49), K_S0πππ (-28,27), KKπ (-24,23) MeV")
alg_m1.note(:mbc_signal_region, "M_BC ∈ (1.863, 1.877) GeV/c² signal region applied in ROOT")
alg_m1.note(:st_yield_fit, "ST yields from ML fit to M_BC distributions; ARGUS background")
alg_m1.note(:bf_extraction, "BF = N_DT / Σ(N_ST^i × ε_DT^i / ε_ST^i); DT fit in ROOT")

alg_m1.apply
alg_m1.execute_on([data_3773, incMC_3773, exMC_m1])

# ═══════════════════════════════════════════════════════════════════════
# Mode 2: D+ → φπ+π+π-, φ→K+K-
# Same final state as mode 1; DTagAlg mode :DptoKKPiPiPi (φ indistinguishable)
# ═══════════════════════════════════════════════════════════════════════
decay_card_m2 = <<~DECAYCARD
  Decay D+
  1.000 phi pi+ pi+ pi- PHSP;
  Enddecay
  Decay phi
  1.000 K+ K- VSS;
  Enddecay
DECAYCARD

exMC_m2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Dp_phi_pipipi"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_m2
  config.cross_section   = :default
end

alg_m2 = TagAnalysis.new("DpToPhiPiPiPi")
alg_m2.set_header(["DpToPhiPiPiPiAlg/DpToPhiPiPiPi.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })

alg_m2.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :deltaE, abs: 0.06
end

alg_m2.tag_side(:Dplus) do |t|
  t.modes :DptoKKPiPiPi       # φ→K+K- collapses to KK in DTagAlg naming
  t.charm 1
end

alg_m2.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_m2.note(:phi_intermediate, "φ→K+K- intermediate resonance; DTagAlg :DptoKKPiPiPi does not distinguish φ from non-φ KK; φ mass window analysis in ROOT")
alg_m2.note(:tag_deltaE_windows, "Mode-dependent ΔE windows from Table 2")
alg_m2.note(:mbc_signal_region, "M_BC ∈ (1.863, 1.877) GeV/c²")
alg_m2.note(:bf_extraction, "BF = N_DT / Σ(N_ST^i × ε_DT^i / ε_ST^i)")

alg_m2.apply
alg_m2.execute_on([data_3773, incMC_3773, exMC_m2])

# ═══════════════════════════════════════════════════════════════════════
# Modes 3-5: D+ → K_S0K+π+π-π0 / K_S0K+η / K_S0K+ω
# Tag modes NOT in authoritative DTagAlg vocabulary — unavailable.
# ═══════════════════════════════════════════════════════════════════════

# Mode 3: D+ → K_S0K+π+π-π0 (K_S0→π+π-, π0→γγ)
# Required tag mode: would be D+→K_S0 K+ π+ π- π0, no matching vocabulary entry.
# Mode 4: D+ → K_S0K+η (η→γγ)
# Required tag mode: would be D+→K_S0 K+ η, no matching vocabulary entry.
# Mode 5: D+ → K_S0K+ω (ω→π+π-π0, π0→γγ)
# Required tag mode: would be D+→K_S0 K+ ω, no matching vocabulary entry.

# Placeholder algorithm with notes for modes 3-5
alg_m3 = Algorithm.new("DpToKsKPipiPi0")
alg_m3.set_header(["DpToKsKPipiPi0Alg/DpToKsKPipiPi0.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })
alg_m3.note(:tag_mode_unavailable, "D+→K_S0K+π+π-π0 tag mode not in DTagAlg vocabulary; cannot express as TagAnalysis")
alg_m3.note(:tag_mode_unavailable, "D+→K_S0K+η tag mode not in DTagAlg vocabulary; cannot express as TagAnalysis")
alg_m3.note(:tag_mode_unavailable, "D+→K_S0K+ω tag mode not in DTagAlg vocabulary; cannot express as TagAnalysis")
alg_m3.note(:paper_signal_modes, "Three signal modes (K_S0K+π+π-π0, K_S0K+η, K_S0K+ω) use DT method at ψ(3770) with same 6 D- tag modes; BF extraction via N_DT/(N_ST×ε_DT/ε_ST)")
# No execute_on — placeholder only