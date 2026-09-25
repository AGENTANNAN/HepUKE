# DSL for arXiv:2110.10999v2
# Observation of D+ -> K+ pi0 pi0 and D+ -> K+ pi0 eta DCS decays at psi(3770)
# TagAnalysis (DT): D- tagged via 3 hadronic modes, signal D+ reconstructed from remaining particles
# Dataset: psi(3770) at 3.773 GeV

psipp_data = DatasetManager.real_data.find("712_3773")
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")

# === Signal mode 1: D+ -> K+ pi0 pi0 ===
alg_kpi0pi0 = TagAnalysis.new("DpToKPi0Pi0")
alg_kpi0pi0.set_header(["DpToKPi0Pi0Alg/DpToKPi0Pi0.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            .note(:dcs_decay, "Doubly Cabibbo-suppressed decay; signal mode not in standard DTagAlg channel list")
            .note(:ks0_veto, "KS0 veto: pi0 pi0 invariant mass outside [0.388, 0.588] GeV/c^2 to reject D+ -> K+ KS0(->pi0pi0); applied in ROOT")
            .note(:opening_angle, "D+ D- opening angle > 167 deg to suppress non-DDbar background; applied in ROOT")
            .note(:sub_resonance, "K*+ -> K+ pi0 sub-resonance studied via simultaneous fit in K*+ signal/sideband regions; ROOT analysis")

# Tag side: D- reconstructed via 3 hadronic DTagAlg modes
alg_kpi0pi0.tag_side(:D) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0
  t.charm -1
end

# Signal side: D+ -> K+ pi0 pi0
alg_kpi0pi0.signal_side do |s|
  s.photons 4          # 2 pi0 -> 4 photons
  s.charged(kp: 1)    # K+
end

# 4C kinematic fit: tag D- + signal particles constrained to ECMS
alg_kpi0pi0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_kpi0pi0.dtag_reconstruction do |d|
  d.beam_energy :db
  d.local true
end

alg_kpi0pi0.apply

# === Signal mode 2: D+ -> K+ pi0 eta ===
alg_kpi0eta = TagAnalysis.new("DpToKPi0Eta")
alg_kpi0eta.set_header(["DpToKPi0EtaAlg/DpToKPi0Eta.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            .note(:dcs_decay, "Doubly Cabibbo-suppressed decay; signal mode not in standard DTagAlg channel list")
            .note(:opening_angle, "D+ D- opening angle > 167 deg; applied in ROOT")
            .note(:sub_resonance, "K*+ -> K+ pi0 sub-resonance via simultaneous fit; ROOT analysis")

alg_kpi0eta.tag_side(:D) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0
  t.charm -1
end

alg_kpi0eta.signal_side do |s|
  s.photons 2          # pi0 -> 2 photons (eta also decays to 2 photons)
  s.charged(kp: 1)    # K+
end

alg_kpi0eta.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_kpi0eta.dtag_reconstruction do |d|
  d.beam_energy :db
  d.local true
end

alg_kpi0eta.apply

# Execute both algorithms on the same datasets
alg_kpi0pi0.execute_on([psipp_data, psipp_incMC])
alg_kpi0eta.execute_on([psipp_data, psipp_incMC])