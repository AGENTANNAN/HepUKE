# =============================================================================
# BOSS / DSL spec — tag-based (double-tag style) search for rare D -> h(h') e+e-
# at psi(3770), sqrt(s) = 3.773 GeV.  Real data + inclusive MC (712_3773) plus
# 200k-event exclusive MC for each of the eleven signal modes.
# =============================================================================

### --------------------------- Datasets --------------------------- ###
data_3773  = DatasetManager.real_data.find("712_3773")      # psi(3770) real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC

### ------------------------- Decay cards ------------------------- ###
# ---- D+ signal modes (tag side: D-) ----
card_dp_pipi0ee = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D- PHSP;
  Enddecay
  Decay D+
  1.0 pi+ pi0 e+ e- PHSP;
  Enddecay
  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay
  Decay D-
  1.0 K+ pi- pi- PHSP;   # representative hadronic tag decay
  Enddecay
  End
DECAYCARD

card_dp_kpi0ee = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D- PHSP;
  Enddecay
  Decay D+
  1.0 K+ pi0 e+ e- PHSP;
  Enddecay
  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay
  Decay D-
  1.0 K+ pi- pi- PHSP;
  Enddecay
  End
DECAYCARD

card_dp_kspiee = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D- PHSP;
  Enddecay
  Decay D+
  1.0 K_S0 pi+ e+ e- PHSP;
  Enddecay
  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay
  Decay D-
  1.0 K+ pi- pi- PHSP;
  Enddecay
  End
DECAYCARD

card_dp_kskeee = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D- PHSP;
  Enddecay
  Decay D+
  1.0 K_S0 K+ e+ e- PHSP;
  Enddecay
  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay
  Decay D-
  1.0 K+ pi- pi- PHSP;
  Enddecay
  End
DECAYCARD

# ---- D0 signal modes (tag side: anti-D0) ----
card_d0_kkee = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0 K- K+ e+ e- PHSP;
  Enddecay
  Decay anti-D0
  1.0 K+ pi- pi0 PHSP;   # representative hadronic tag decay
  Enddecay
  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

card_d0_pipimee = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0 pi+ pi- e+ e- PHSP;
  Enddecay
  Decay anti-D0
  1.0 K+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

card_d0_kpiee = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0 K- pi+ e+ e- PHSP;
  Enddecay
  Decay anti-D0
  1.0 K+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

card_d0_pi0ee = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0 pi0 e+ e- PHSP;
  Enddecay
  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay
  Decay anti-D0
  1.0 K+ pi- pi0 PHSP;
  Enddecay
  End
DECAYCARD

card_d0_etaee = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0 eta e+ e- PHSP;
  Enddecay
  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay
  Decay anti-D0
  1.0 K+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

card_d0_omegaee = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0 omega e+ e- PHSP;
  Enddecay
  Decay omega
  1.0 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay
  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay
  Decay anti-D0
  1.0 K+ pi- pi0 PHSP;
  Enddecay
  End
DECAYCARD

card_d0_ksee = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0 K_S0 e+ e- PHSP;
  Enddecay
  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay
  Decay anti-D0
  1.0 K+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

### ------------------- Exclusive MC (200k events) ------------------- ###
exMC_dp_pipi0ee = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_Dp_pipi0ee"; c.related_dataset = data_3773
  c.events = 200_000; c.decay_card = card_dp_pipi0ee; c.cross_section = :default
end
exMC_dp_kpi0ee = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_Dp_kpi0ee"; c.related_dataset = data_3773
  c.events = 200_000; c.decay_card = card_dp_kpi0ee; c.cross_section = :default
end
exMC_dp_kspiee = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_Dp_kspiee"; c.related_dataset = data_3773
  c.events = 200_000; c.decay_card = card_dp_kspiee; c.cross_section = :default
end
exMC_dp_kskeee = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_Dp_kskeee"; c.related_dataset = data_3773
  c.events = 200_000; c.decay_card = card_dp_kskeee; c.cross_section = :default
end
exMC_d0_kkee = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_D0_kkee"; c.related_dataset = data_3773
  c.events = 200_000; c.decay_card = card_d0_kkee; c.cross_section = :default
end
exMC_d0_pipimee = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_D0_pipimee"; c.related_dataset = data_3773
  c.events = 200_000; c.decay_card = card_d0_pipimee; c.cross_section = :default
end
exMC_d0_kpiee = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_D0_kpiee"; c.related_dataset = data_3773
  c.events = 200_000; c.decay_card = card_d0_kpiee; c.cross_section = :default
end
exMC_d0_pi0ee = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_D0_pi0ee"; c.related_dataset = data_3773
  c.events = 200_000; c.decay_card = card_d0_pi0ee; c.cross_section = :default
end
exMC_d0_etaee = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_D0_etaee"; c.related_dataset = data_3773
  c.events = 200_000; c.decay_card = card_d0_etaee; c.cross_section = :default
end
exMC_d0_omegaee = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_D0_omegaee"; c.related_dataset = data_3773
  c.events = 200_000; c.decay_card = card_d0_omegaee; c.cross_section = :default
end
exMC_d0_ksee = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "exmc_3773_D0_ksee"; c.related_dataset = data_3773
  c.events = 200_000; c.decay_card = card_d0_ksee; c.cross_section = :default
end

### ==================== Tag-based event selection ==================== ###
# The tag D is reconstructed from pre-stored DTag candidates (single-tag), the
# opposite-charm D is fully reconstructed as the signal side.  Store-not-cut:
# the tag DeltaE is windowed only because the description asks for it.

# Common signal-side / fit notes reused by every mode.
COMMON_NOTES = {
  background_veto: "phi veto: reject events with M(e+e-) in (0.935, 1.053) GeV/c^2 to suppress gamma-conversion background",
  pid_correction_method: "electron PID uses the fixed v1 lepton thresholds of the tag framework; the additional E/pc > 0.8 and combined e-PID criteria are not implemented in the DSL"
}
VERTEX_NOTE = "e+e- pair required to have vertex R_xy outside (2.0, 8.0) cm (gamma-conversion suppression)"

def tag_side_dplus(t)
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.window :deltaE, min: -0.060, max: 0.034
end

def tag_side_d0(t)
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1                      # require the tag D0 and signal D0 to have opposite charm
  t.window :deltaE, min: -0.064, max: 0.035
end

# ---------------------------- D+ signal modes ---------------------------- #
alg = TagAnalysis.new("DpRarePipi0EE")
alg.set_header(["DpRarePipi0EEAlg/DpRarePipi0EE.h"]).set_constant({"ECMS" => [:double, 3.773]}).with_decay_card(card_dp_pipi0ee)
alg.tag_side(:Dplus) { |t| tag_side_dplus(t) }
alg.signal_side { |s| s.charged(pip: 1, ep: 1, em: 1); s.photons 2 }
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end
COMMON_NOTES.each { |k, v| alg.note(k, v) }
alg.note(:background_veto, VERTEX_NOTE)
alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_dp_pipi0ee])

alg = TagAnalysis.new("DpRareKpi0EE")
alg.set_header(["DpRareKpi0EEAlg/DpRareKpi0EE.h"]).set_constant({"ECMS" => [:double, 3.773]}).with_decay_card(card_dp_kpi0ee)
alg.tag_side(:Dplus) { |t| tag_side_dplus(t) }
alg.signal_side { |s| s.charged(kp: 1, ep: 1, em: 1); s.photons 2 }
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end
COMMON_NOTES.each { |k, v| alg.note(k, v) }
alg.note(:background_veto, VERTEX_NOTE)
alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_dp_kpi0ee])

alg = TagAnalysis.new("DpRareKsPiEE")
alg.set_header(["DpRareKsPiEEAlg/DpRareKsPiEE.h"]).set_constant({"ECMS" => [:double, 3.773]}).with_decay_card(card_dp_kspiee)
alg.tag_side(:Dplus) { |t| tag_side_dplus(t) }
# K_S0 -> pi+pi- contributes 2 tracks; signal content = K_S0 pi+ e+ e-
alg.signal_side { |s| s.charged(pip: 2, pim: 1, ep: 1, em: 1) }
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end
COMMON_NOTES.each { |k, v| alg.note(k, v) }
alg.note(:background_veto, VERTEX_NOTE)
alg.note(:ks0_reconstruction, "signal-side K_S0 reconstructed via displaced vertex, L/sigma_L > 2; not expressible in signal_side")
alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_dp_kspiee])

alg = TagAnalysis.new("DpRareKsKEE")
alg.set_header(["DpRareKsKEEAlg/DpRareKsKEE.h"]).set_constant({"ECMS" => [:double, 3.773]}).with_decay_card(card_dp_kskeee)
alg.tag_side(:Dplus) { |t| tag_side_dplus(t) }
alg.signal_side { |s| s.charged(kp: 1, pip: 1, pim: 1, ep: 1, em: 1) }
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end
COMMON_NOTES.each { |k, v| alg.note(k, v) }
alg.note(:background_veto, VERTEX_NOTE)
alg.note(:ks0_reconstruction, "signal-side K_S0 reconstructed via displaced vertex, L/sigma_L > 2; not expressible in signal_side")
alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_dp_kskeee])

# ---------------------------- D0 signal modes ---------------------------- #
alg = TagAnalysis.new("D0RareKKEE")
alg.set_header(["D0RareKKEEAlg/D0RareKKEE.h"]).set_constant({"ECMS" => [:double, 3.773]}).with_decay_card(card_d0_kkee)
alg.tag_side(:D0) { |t| tag_side_d0(t) }
alg.signal_side { |s| s.charged(km: 1, kp: 1, ep: 1, em: 1) }
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end
COMMON_NOTES.each { |k, v| alg.note(k, v) }
alg.note(:background_veto, VERTEX_NOTE)
alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_d0_kkee])

alg = TagAnalysis.new("D0RarePiPiEE")
alg.set_header(["D0RarePiPiEEAlg/D0RarePiPiEE.h"]).set_constant({"ECMS" => [:double, 3.773]}).with_decay_card(card_d0_pipimee)
alg.tag_side(:D0) { |t| tag_side_d0(t) }
alg.signal_side { |s| s.charged(pip: 1, pim: 1, ep: 1, em: 1) }
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end
COMMON_NOTES.each { |k, v| alg.note(k, v) }
alg.note(:background_veto, VERTEX_NOTE)
alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_d0_pipimee])

alg = TagAnalysis.new("D0RareKPiEE")
alg.set_header(["D0RareKPiEEAlg/D0RareKPiEE.h"]).set_constant({"ECMS" => [:double, 3.773]}).with_decay_card(card_d0_kpiee)
alg.tag_side(:D0) { |t| tag_side_d0(t) }
alg.signal_side { |s| s.charged(km: 1, pip: 1, ep: 1, em: 1) }
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end
COMMON_NOTES.each { |k, v| alg.note(k, v) }
alg.note(:background_veto, VERTEX_NOTE)
alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_d0_kpiee])

alg = TagAnalysis.new("D0RarePi0EE")
alg.set_header(["D0RarePi0EEAlg/D0RarePi0EE.h"]).set_constant({"ECMS" => [:double, 3.773]}).with_decay_card(card_d0_pi0ee)
alg.tag_side(:D0) { |t| tag_side_d0(t) }
alg.signal_side { |s| s.charged(ep: 1, em: 1); s.photons 2 }
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end
COMMON_NOTES.each { |k, v| alg.note(k, v) }
alg.note(:background_veto, VERTEX_NOTE)
alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_d0_pi0ee])

alg = TagAnalysis.new("D0RareEtaEE")
alg.set_header(["D0RareEtaEEAlg/D0RareEtaEE.h"]).set_constant({"ECMS" => [:double, 3.773]}).with_decay_card(card_d0_etaee)
alg.tag_side(:D0) { |t| tag_side_d0(t) }
alg.signal_side { |s| s.charged(ep: 1, em: 1); s.photons 2 }
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)   # eta -> gamma gamma
  f.chi2_cut 200
end
COMMON_NOTES.each { |k, v| alg.note(k, v) }
alg.note(:background_veto, VERTEX_NOTE)
alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_d0_etaee])

alg = TagAnalysis.new("D0RareOmegaEE")
alg.set_header(["D0RareOmegaEEAlg/D0RareOmegaEE.h"]).set_constant({"ECMS" => [:double, 3.773]}).with_decay_card(card_d0_omegaee)
alg.tag_side(:D0) { |t| tag_side_d0(t) }
alg.signal_side { |s| s.charged(pip: 1, pim: 1, ep: 1, em: 1); s.photons 2 }   # omega -> pi+pi-pi0, pi0 -> gamma gamma
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)       # constrain the pi0 from omega
  f.chi2_cut 200
end
COMMON_NOTES.each { |k, v| alg.note(k, v) }
alg.note(:background_veto, VERTEX_NOTE)
alg.note(:omega_mass_window, "require M(pi+pi-pi0) in (0.720, 0.840) GeV/c^2 before the kinematic fit (signal-side omega selection; not expressible in signal_side)")
alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_d0_omegaee])

alg = TagAnalysis.new("D0RareKsEE")
alg.set_header(["D0RareKsEEAlg/D0RareKsEE.h"]).set_constant({"ECMS" => [:double, 3.773]}).with_decay_card(card_d0_ksee)
alg.tag_side(:D0) { |t| tag_side_d0(t) }
# K_S0 -> pi+pi- contributes 2 tracks; signal content = K_S0 e+ e-
alg.signal_side { |s| s.charged(pip: 1, pim: 1, ep: 1, em: 1) }
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end
COMMON_NOTES.each { |k, v| alg.note(k, v) }
alg.note(:background_veto, VERTEX_NOTE)
alg.note(:ks0_reconstruction, "signal-side K_S0 reconstructed via displaced vertex, L/sigma_L > 2; not expressible in signal_side")
alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_d0_ksee])