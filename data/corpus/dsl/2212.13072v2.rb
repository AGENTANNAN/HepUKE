# DSL for 2212.13072v2: D_s^+ → π^+π^+π^-X inclusive BF at 4.178 GeV
# Tag-based: D_s DT with ST modes K_S^0K^- and K^-K^+π^-
# Signal side: π^+π^+π^-X inclusive (K_S^0 excluded, electron veto)
# Efficiency matrix unfolding in 11 M_3π intervals

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_4180 = DatasetManager.real_data.find("703_4180")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")

# ============================================================
# Tag Mode 1: D_s^- → K_S^0 K^-  (ST + inclusive 3π signal)
# ============================================================

decay_card_st1 = <<~DECAY
Decay e+ e-
  1.0  D_s+  D_s-  VSS;
Enddecay
Decay D_s+
  1.0  pi+  pi+  pi-  X  PHSP;
Enddecay
Decay D_s-
  1.0  K_S0  K-  PHSP;
Enddecay
Decay K_S0
  1.0  pi+  pi-  PHSP;
Enddecay
DECAY

sigMC_st1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Ds_3piX_ST_KsK"
  config.related_dataset = data_4180
  config.events          = 500_000
  config.decay_card      = decay_card_st1
  config.cross_section   = :default
end

alg_st1 = TagAnalysis.new("DsTag3piX_KsK")
alg_st1.set_header(["DsTag3piXAlg/DsTag3piX.h"])
        .set_constant({ "ECMS" => [:double, 4.178] })

alg_st1.tag_side(:Ds) do |t|
  t.modes :DstoKsK
  t.charm -1
end

alg_st1.signal_side do |s|
  s.charged(pip: 2, pim: 1)
  s.require_charge 1
end

alg_st1.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_st1.note(:inclusive, "Inclusive signal: any event where the D_s^+ signal side has at least π^+π^+π^-X. Efficiency matrix unfolding in 11 bins of M_3π from 0.55 to 2.15 GeV/c^2")
alg_st1.note(:ks_veto, "K_S^0 veto on signal side: reject events where any π^+π^- pair falls in [0.485, 0.510] GeV/c^2")
alg_st1.note(:electron_veto, "Electron veto on signal side pions: E/p < 0.8 to reject e^+e^--contaminated tracks")
alg_st1.note(:tag_selection, "ST D_s^- tagged via K_S^0K^- and K^-K^+π^-. M_BC and ΔE windows: M_BC > 2.05 GeV/c^2 for K_S^0K^- mode. Best tag by minimum |ΔE|")
alg_st1.note(:efficiency_matrix, "11×11 efficiency matrix accounts for bin migration in M_3π. Unfolding via SVD method with L-curve optimization")
alg_st1.note(:bias_study, "Efficiency and unfolding validated using inclusive MC and mixed-signal generic MC")
alg_st1.note(:beam_energy, "4.178 GeV at BOSS 703. Luminosity 3.19 fb^{-1}")

alg_st1.apply
alg_st1.execute_on([data_4180, incMC_4180, sigMC_st1])

# ============================================================
# Tag Mode 2: D_s^- → K^- K^+ π^-  (ST + inclusive 3π signal)
# ============================================================

decay_card_st2 = <<~DECAY
Decay e+ e-
  1.0  D_s+  D_s-  VSS;
Enddecay
Decay D_s+
  1.0  pi+  pi+  pi-  X  PHSP;
Enddecay
Decay D_s-
  1.0  K-  K+  pi-  PHSP;
Enddecay
DECAY

sigMC_st2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Ds_3piX_ST_KKpi"
  config.related_dataset = data_4180
  config.events          = 500_000
  config.decay_card      = decay_card_st2
  config.cross_section   = :default
end

alg_st2 = TagAnalysis.new("DsTag3piX_KKpi")
alg_st2.set_header(["DsTag3piXAlg/DsTag3piX.h"])
        .set_constant({ "ECMS" => [:double, 4.178] })

alg_st2.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

alg_st2.signal_side do |s|
  s.charged(pip: 2, pim: 1)
  s.require_charge 1
end

alg_st2.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_st2.note(:inclusive, "Inclusive signal: any event where the D_s^+ signal side has at least π^+π^+π^-X. Efficiency matrix unfolding in 11 bins of M_3π from 0.55 to 2.15 GeV/c^2")
alg_st2.note(:ks_veto, "K_S^0 veto on signal side: reject events where any π^+π^- pair falls in [0.485, 0.510] GeV/c^2")
alg_st2.note(:electron_veto, "Electron veto on signal side pions: E/p < 0.8 to reject e^+e^--contaminated tracks")
alg_st2.note(:tag_selection, "ST D_s^- tagged via K_S^0K^- and K^-K^+π^-. M_BC and ΔE windows: M_BC > 2.05 GeV/c^2 for K_S^0K^- mode. Best tag by minimum |ΔE|")
alg_st2.note(:efficiency_matrix, "11×11 efficiency matrix accounts for bin migration in M_3π. Unfolding via SVD method with L-curve optimization")
alg_st2.note(:bias_study, "Efficiency and unfolding validated using inclusive MC and mixed-signal generic MC")
alg_st2.note(:beam_energy, "4.178 GeV at BOSS 703. Luminosity 3.19 fb^{-1}")

alg_st2.apply
alg_st2.execute_on([data_4180, incMC_4180, sigMC_st2])