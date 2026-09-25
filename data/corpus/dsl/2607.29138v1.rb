# Dataset: BESIII 1667 pb^-1 e+e- data at sqrt(s) = 4.682 GeV (BOSS 706)
data_4680  = DatasetManager.real_data.find("706_4680")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")

# Signal decay card for e+e- -> D0 D- pi+ (charge conjugate implied)
# Analysed via single-tag (ST) and double-tag (DT) reconstruction:
#   D0 tag modes:   K- pi+ ; K- pi+ pi+ pi- ; K- pi+ pi0
#   D-  tag modes:  K+ pi- pi- ; Ks0 pi- ; K+ pi- pi- pi0 ; Ks0 pi- pi0 ; Ks0 pi- pi- pi+ ; K+ pi- pi- (K+ K- pi-)
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D0 D- pi+                     PHSP;
  Enddecay

  Decay D0
  0.3333 K- pi+                        PHSP;
  0.3333 K- pi+ pi+ pi-                PHSP;
  0.3334 K- pi+ pi0                    D_DALITZ;
  Enddecay

  Decay D-
  0.2000 K+ pi- pi-                    D_DALITZ;
  0.2000 K_S0 pi-                      PHSP;
  0.2000 K+ pi- pi- pi0                PHSP;
  0.2000 K_S0 pi- pi0                  PHSP;
  0.2000 K_S0 pi- pi- pi+              PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-                       PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                   PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "ee_D0DmPi_4682"
  config.related_dataset = data_4680
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# -----------------------------------------------------------------------
# Double-tag (DT) analysis: fully reconstruct D0, D-, and the direct pi+
# using the tag-based DSL layer (D0 tag + D- tag; no missing particle).
# -----------------------------------------------------------------------
alg_dt = TagAnalysis.new("D0DmPiDT")
alg_dt.set_header(["D0DmPiDTAlg/D0DmPiDT.h"])
      .set_constant({ "ECMS" => [:double, 4.682] })
      .with_decay_card(decay_card_signal)

alg_dt.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPiPi, :D0toKPiPi0
  t.charm 1                # D0 (not anti-D0)
end

alg_dt.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm -1               # D- (D+ with charm = -1)
  t.rank_by :inv
end

alg_dt.signal_side do |s|
  s.charged(pip: 1)        # the direct pi+ from e+e- annihilation
  s.require_charge 1
end

alg_dt.fit do |f|
  f.constrain_four_momentum    # 4C kinematic fit; paper requires chi2 < 40
  f.chi2_cut 200               # loose BOSS-level cut; tight cut applied in ROOT
end

alg_dt.note(:pi_d_soft_pion,
            "The signal-side pi+ (pi_d) is the direct pion from the e+e- -> D0 D- pi+ vertex; " \
            "no additional charged tracks are allowed in the DT event.")
alg_dt.note(:dstar_veto,
            "Background from e+e- -> D*+ D- suppressed by requiring RM(D-) > 2.04 GeV/c^2 " \
            "(equivalently M(D0 pi_d+) > 2.04 GeV/c^2).")
alg_dt.note(:d_mass_constraint,
            "In the DT sample, invariant masses of the D0 and D- are constrained to their " \
            "nominal PDG values to improve the mass resolution used in the PWA.")

alg_dt.apply
alg_dt.execute_on([data_4680, incMC_4680, exMC_signal])

# -----------------------------------------------------------------------
# Single-tag (ST) samples: tag one D meson + the direct pi+; the other D
# meson is inferred from the recoil (missing).
# Three ST categories: ST^0_Kpi (D0->Kpi), ST^0_K3pi (D0->K3pi), ST^-_K2pi (D- ->K+pi-pi-).
# -----------------------------------------------------------------------

# ST^0_Kpi : tag D0 -> K- pi+, direct pi+ reconstructed, missing D-
alg_st_D0Kpi = TagAnalysis.new("ST_D0Kpi_MissDm")
alg_st_D0Kpi.set_header(["ST_D0Kpi_MissDmAlg/ST_D0Kpi_MissDm.h"])
            .set_constant({ "ECMS" => [:double, 4.682] })
            .with_decay_card(decay_card_signal)

alg_st_D0Kpi.tag_side(:D0) do |t|
  t.modes :D0toKPi
  t.charm 1
end

alg_st_D0Kpi.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge 0            # +1 (pi+) + (-1) missing D- = 0
  s.missing :"D-"
end

alg_st_D0Kpi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_st_D0Kpi.apply
alg_st_D0Kpi.execute_on([data_4680, incMC_4680, exMC_signal])

# ST^0_K3pi : tag D0 -> K- pi+ pi+ pi-, direct pi+ reconstructed, missing D-
alg_st_D0K3pi = TagAnalysis.new("ST_D0K3pi_MissDm")
alg_st_D0K3pi.set_header(["ST_D0K3pi_MissDmAlg/ST_D0K3pi_MissDm.h"])
             .set_constant({ "ECMS" => [:double, 4.682] })
             .with_decay_card(decay_card_signal)

alg_st_D0K3pi.tag_side(:D0) do |t|
  t.modes :D0toKPiPiPi
  t.charm 1
end

alg_st_D0K3pi.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge 0
  s.missing :"D-"
end

alg_st_D0K3pi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_st_D0K3pi.note(:comb_bkg_veto,
                   "For ST_K3pi sample, combinatorial backgrounds from e+e- -> D0bar D0 pi+ and " \
                   "e+e- -> pi+ pi- D+ D- suppressed by requiring " \
                   "|M(K- pi_d+ pi_{1/2}+ pi-) - m_D0| > 10 MeV/c^2 and " \
                   "|M(K- pi_d+ pi_{1/2}+) - m_D+| > 16 MeV/c^2.")

alg_st_D0K3pi.apply
alg_st_D0K3pi.execute_on([data_4680, incMC_4680, exMC_signal])

# ST^-_K2pi : tag D- -> K+ pi- pi-, direct pi+ reconstructed, missing D0
alg_st_DmK2pi = TagAnalysis.new("ST_DmK2pi_MissD0")
alg_st_DmK2pi.set_header(["ST_DmK2pi_MissD0Alg/ST_DmK2pi_MissD0.h"])
             .set_constant({ "ECMS" => [:double, 4.682] })
             .with_decay_card(decay_card_signal)

alg_st_DmK2pi.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi
  t.charm -1
end

alg_st_DmK2pi.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge 0
  s.missing :D0
end

alg_st_DmK2pi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_st_DmK2pi.apply
alg_st_DmK2pi.execute_on([data_4680, incMC_4680, exMC_signal])
