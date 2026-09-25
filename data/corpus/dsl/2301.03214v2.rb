# DSL for 2301.03214v2: D^0/D^+ → π^+π^+π^-X inclusive BFs at ψ(3770)
# Tag-based: D^0/D^+ DT with ST modes K^+π^- (flavor) and K^+π^-π^-
# Signal side: π^+π^+π^-X inclusive (K_S^0 veto)
# QC correction factor for D^0 tag, efficiency matrix unfolding
# 2.93 fb^{-1} at ψ(3770), sample 712_3773

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================
# D^0 analysis: tag anti-D^0 → K^+π^- (flavor tag)
# Signal: D^0 → π^+π^+π^-X inclusive
# ============================================================

decay_card_d0 = <<~DECAY
Decay psi(3770)
  1.0  D0  anti-D0  VSS;
Enddecay
Decay D0
  1.0  pi+  pi+  pi-  X  PHSP;
Enddecay
Decay anti-D0
  1.0  K+  pi-  PHSP;
Enddecay
DECAY

sigMC_d0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_D0_3piX_flavorTag"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_d0
  config.cross_section   = :default
end

alg_d0 = TagAnalysis.new("D0Tag3piX")
alg_d0.set_header(["D0Tag3piXAlg/D0Tag3piX.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })

alg_d0.tag_side(:D0) do |t|
  t.modes :D0toKPi
  t.charm -1
end

alg_d0.signal_side do |s|
  s.charged(pip: 2, pim: 1, at_least: true)
  s.require_charge 1
end

alg_d0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0.note(:inclusive, "Inclusive signal: D^0 signal side must have at least π^+π^+π^- X. Efficiency matrix unfolding in 9 bins of M_3π from 0.50 to 2.05 GeV/c^2")
alg_d0.note(:ks_veto, "K_S^0 veto on signal side: reject events where any π^+π^- pair has invariant mass in [0.485, 0.510] GeV/c^2")
alg_d0.note(:qc_correction, "QC (quantum coherence) correction factor for D^0 tag at ψ(3770): accounts for C=-1 initial state correlation. Applied as multiplicative correction to observed yields")
alg_d0.note(:tag_mode, "Flavor tag: anti-D^0 → K^+π^-. DCS contamination negligible for inclusive measurement")
alg_d0.note(:tag_selection, "Tag side M_BC and ΔE windows applied. Best tag candidate by minimum |ΔE|")
alg_d0.note(:efficiency_matrix, "9×9 efficiency matrix for bin migration in M_3π. Unfolding via SVD method with L-curve regularization")
alg_d0.note(:beam_energy, "ψ(3770) at 3.773 GeV, BOSS 712. Luminosity 2.93 fb^{-1}")

alg_d0.apply
alg_d0.execute_on([data_3773, incMC_3773, sigMC_d0])

# ============================================================
# D^+ analysis: tag D^- → K^+π^-π^- (hadronic tag)
# Signal: D^+ → π^+π^+π^-X inclusive
# ============================================================

decay_card_dp = <<~DECAY
Decay psi(3770)
  1.0  D+  D-  VSS;
Enddecay
Decay D+
  1.0  pi+  pi+  pi-  X  PHSP;
Enddecay
Decay D-
  1.0  K+  pi-  pi-  PHSP;
Enddecay
DECAY

sigMC_dp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Dp_3piX_tagKKpipi"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_dp
  config.cross_section   = :default
end

alg_dp = TagAnalysis.new("DpTag3piX")
alg_dp.set_header(["DpTag3piXAlg/DpTag3piX.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })

alg_dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi
  t.charm -1
end

alg_dp.signal_side do |s|
  s.charged(pip: 2, pim: 1, at_least: true)
  s.require_charge 1
end

alg_dp.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_dp.note(:inclusive, "Inclusive signal: D^+ signal side must have at least π^+π^+π^- X. Efficiency matrix unfolding in 10 bins of M_3π from 0.55 to 2.15 GeV/c^2")
alg_dp.note(:ks_veto, "K_S^0 veto on signal side: reject events where any π^+π^- pair has invariant mass in [0.485, 0.510] GeV/c^2")
alg_dp.note(:tag_mode, "Hadronic tag: D^- → K^+π^-π^-. No quantum correlation effects for charged D mesons")
alg_dp.note(:tag_selection, "Tag side M_BC and ΔE windows applied. Best tag candidate by minimum |ΔE|")
alg_dp.note(:efficiency_matrix, "10×10 efficiency matrix for bin migration in M_3π. Unfolding via SVD method with L-curve regularization")
alg_dp.note(:beam_energy, "ψ(3770) at 3.773 GeV, BOSS 712. Luminosity 2.93 fb^{-1}")
alg_dp.note(:combined_fit, "D^0 and D^+ results combined for isospin-averaged partial BFs")

alg_dp.apply
alg_dp.execute_on([data_3773, incMC_3773, sigMC_dp])