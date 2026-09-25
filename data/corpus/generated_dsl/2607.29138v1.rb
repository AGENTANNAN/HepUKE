# =====================================================================
# e+e- -> D0 D- pi+  at sqrt(s) = 4.682 GeV  (sample 706_4680)
# Single-tag (ST) and double-tag (DT) analyses with DTag (evtRecDTag).
# ST samples: ST0_Kpi, ST0_K3pi, STm_K2pi ; DT = D0 tag x D- tag x direct pi+.
# No explicit track / photon / PID cut values are given, so the tag
# algorithms' standard object reconstruction is used throughout.
# =====================================================================

### Dataset preparation ###
data_4680  = DatasetManager.real_data.find("706_4680")          # real data, sqrt(s) = 4.682 GeV
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")       # corresponding inclusive MC

# Signal decay table: e+e- -> D0 D- pi+ with all tag modes of D0 and D-
# (relative weights used for the listed modes only).
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D0 D- pi+ PHSP;
    Enddecay

    Decay D0
    0.15 K- pi+ PHSP;
    0.31 K- pi+ pi+ pi- PHSP;
    0.54 K- pi+ pi0 PHSP;
    Enddecay

    Decay D-
    0.20 K+ pi- pi- PHSP;
    0.20 K_S0 pi- PHSP;
    0.20 K+ pi- pi- pi0 PHSP;
    0.20 K_S0 pi- pi0 PHSP;
    0.20 K_S0 pi- pi- pi+ PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal process (500k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_4680_D0Dmpip"
  config.related_dataset = data_4680
  config.events         = 500000
  config.decay_card     = decay_card_signal
  config.cross_section  = :default
end

# =====================================================================
# DT analysis: reconstruct D0 (tag1) and D- (tag2) on the tag sides and
# exactly one direct pi+ on the signal side.
# =====================================================================
alg_dt = TagAnalysis.new("D0DmPiDT")
alg_dt.set_header(["D0DmPiDTAlg/D0DmPiDT.h"])
      .set_constant({"ECMS" => [:double, 4.682]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      .with_decay_card(decay_card_signal)

# Tag side 1: D0, charm = +1
alg_dt.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPiPi, :D0toKPiPi0     # Kpi, K3pi, Kpipi0
  t.charm 1
end

# Tag side 2: D-, charm = -1, ranked by invariant mass
alg_dt.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm(-1)
  t.rank_by :inv
end

# Signal side: exactly one direct pi+, net charge +1, no extra charged tracks
alg_dt.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge 1
end

# 4C kinematic fit, loose BOSS-level chi2 < 200 (paper applies chi2 < 40 in ROOT)
alg_dt.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)    # constrain D0 mass to PDG
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:"D+")  # constrain D- mass to PDG
  f.chi2_cut 200
end

# D* background veto (recoil / invariant-mass observable - applied in ROOT)
alg_dt.note(:background_veto,
  "D* background veto: RM(D-) > 2.04 GeV/c^2 OR M(D0 pi_d+) > 2.04 GeV/c^2; " \
  "the recoil-against-D- and D0 pi_d+ invariant masses are stored and the window " \
  "is applied post-fit in the ROOT analysis")

alg_dt.apply
alg_dt.execute_on([data_4680, incMC_4680, exMC_signal])

# =====================================================================
# ST0_Kpi: D0 -> K- pi+ tag, one direct pi+, recoiling D- missing
# =====================================================================
alg_st0_kpi = TagAnalysis.new("ST0KPi")
alg_st0_kpi.set_header(["ST0KPiAlg/ST0KPi.h"])
           .set_constant({"ECMS" => [:double, 4.682]})
           .with_decay_card(decay_card_signal)

alg_st0_kpi.tag_side(:D0) do |t|
  t.modes :D0toKPi
end

alg_st0_kpi.signal_side do |s|
  s.charged(pip: 1)      # one reconstructed direct pi+
  s.missing :Dm          # recoiling D- is not reconstructed
end

alg_st0_kpi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_st0_kpi.apply
alg_st0_kpi.execute_on([data_4680, incMC_4680, exMC_signal])

# =====================================================================
# ST0_K3pi: D0 -> K- pi+ pi+ pi- tag, one direct pi+, recoiling D- missing
# =====================================================================
alg_st0_k3pi = TagAnalysis.new("ST0K3Pi")
alg_st0_k3pi.set_header(["ST0K3PiAlg/ST0K3Pi.h"])
            .set_constant({"ECMS" => [:double, 4.682]})
            .with_decay_card(decay_card_signal)

alg_st0_k3pi.tag_side(:D0) do |t|
  t.modes :D0toKPiPiPi          # D0 -> K- pi+ pi+ pi-
end

alg_st0_k3pi.signal_side do |s|
  s.charged(pip: 1)
  s.missing :Dm
end

alg_st0_k3pi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# Extra background-reflection vetoes for this ST sample (applied in ROOT)
alg_st0_k3pi.note(:background_veto,
  "K3pi ST reflection vetoes: |M(K- pi_d+ pi+ pi-) - m_D0| > 10 MeV/c^2 AND " \
  "|M(K- pi_d+ pi+) - m_D+| > 16 MeV/c^2, where pi_d is the direct pion; " \
  "these masses use the direct pion substituted into the tag decay products and " \
  "are windowed in the ROOT analysis")

alg_st0_k3pi.apply
alg_st0_k3pi.execute_on([data_4680, incMC_4680, exMC_signal])

# =====================================================================
# STm_K2pi: D- -> K+ pi- pi- tag, one direct pi+, recoiling D0 missing
# =====================================================================
alg_stm_k2pi = TagAnalysis.new("STmK2Pi")
alg_stm_k2pi.set_header(["STmK2PiAlg/STmK2Pi.h"])
            .set_constant({"ECMS" => [:double, 4.682]})
            .with_decay_card(decay_card_signal)

alg_stm_k2pi.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi             # D- -> K+ pi- pi-
  t.charm(-1)
end

alg_stm_k2pi.signal_side do |s|
  s.charged(pip: 1)              # one reconstructed direct pi+
  s.missing :D0                  # recoiling D0 is not reconstructed
end

alg_stm_k2pi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"D+")  # D- mass to PDG
  f.chi2_cut 200
end

alg_stm_k2pi.apply
alg_stm_k2pi.execute_on([data_4680, incMC_4680, exMC_signal])