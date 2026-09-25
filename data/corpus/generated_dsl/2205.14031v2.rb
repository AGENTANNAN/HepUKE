# ============================================================
# Datasets — psi(3770), sqrt(s) = 3.773 GeV  (2.93 fb^-1 sample)
# No exclusive-MC generation is configured for this analysis.
# ============================================================
data_3773  = DatasetManager.real_data.find("712_3773")    # psi(3770) real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773") # matching inclusive MC

# ============================================================
# Chain 1 : e+e- -> D0 D0bar  (double tag)
# Tag D0bar in K+pi-, K+pi-pi0, K+pi-pi-pi+
# Signal D0 in K_S0 pi0 pi0 pi0, K- pi+ pi0 pi0 pi0, K_S0 pi+ pi- pi0 pi0
# ============================================================
alg_d0 = TagAnalysis.new("DTagD0")
alg_d0.set_header(["DTagD0Alg/DTagD0.h"])
      .set_constant({"ECMS" => [:double, 3.773]})   # sqrt(s) = 3.773 GeV (1.8865 GeV per beam)

# --- tag side : D0bar ---
alg_d0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1                                        # pin the tagged side to the D0bar
end

# --- signal side of the double tag : D0 ---
alg_d0.tag_side(:D0) do |t|
  t.modes :D0toKsPi0Pi0Pi0, :D0toKPiPi0Pi0Pi0, :D0toKsPiPiPi0Pi0
  t.charm 1                                         # pin the signal side to the D0
end

# BOSS-side reconstruction steps with no dedicated DSL construct
alg_d0.note(:ks0_reconstruction,
            "K_S0 candidates are formed from pi+ pi- pairs by a secondary-vertex fit; the pair invariant mass is required to be in [0.486, 0.510] GeV/c^2")
     .note(:pi0_reconstruction,
            "pi0 candidates are formed from gamma gamma pairs with invariant mass in [0.115, 0.150] GeV/c^2")
     .note(:background_veto,
            "pi+ pi- pairs with invariant mass in [0.468, 0.528] GeV/c^2 are vetoed to suppress K_S0 contamination")

# --- kinematic fit : 4-momentum conservation + both tag and signal D masses
#     constrained to the nominal D0 mass, chi2 < 200 ---
alg_d0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)   # tag D0bar -> M(D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)   # signal D0  -> M(D0)
  f.chi2_cut 200
end

alg_d0.apply                                        # no Selection argument for tag analyses
root_files_d0 = alg_d0.execute_on([data_3773, incMC_3773])

# ============================================================
# Chain 2 : e+e- -> D+ D-  (double tag)
# Tag D- in K+pi-pi-, K_S0 pi-, K+pi-pi-pi0, K_S0 pi- pi0,
#            K_S0 pi+ pi- pi-, K+ K- pi-
# Signal D+ in K_S0 pi+ pi0 pi0, K_S0 pi+ pi+ pi- pi0,
#              K_S0 pi+ pi0 pi0 pi0, K- pi+ pi+ pi0 pi0
# ============================================================
alg_dp = TagAnalysis.new("DTagDp")
alg_dp.set_header(["DTagDpAlg/DTagDp.h"])
      .set_constant({"ECMS" => [:double, 3.773]})   # sqrt(s) = 3.773 GeV (1.8865 GeV per beam)

# --- tag side : D- ---
alg_dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1                                        # pin the tagged side to the D-
end

# --- signal side of the double tag : D+ ---
alg_dp.tag_side(:Dplus) do |t|
  t.modes :DptoKsPiPi0Pi0, :DptoKsPiPiPiPi0, :DptoKsPiPi0Pi0Pi0, :DptoKPiPiPi0Pi0
  t.charm 1                                         # pin the signal side to the D+
end

# BOSS-side reconstruction steps with no dedicated DSL construct
alg_dp.note(:ks0_reconstruction,
            "K_S0 candidates are formed from pi+ pi- pairs by a secondary-vertex fit; the pair invariant mass is required to be in [0.486, 0.510] GeV/c^2")
     .note(:pi0_reconstruction,
            "pi0 candidates are formed from gamma gamma pairs with invariant mass in [0.115, 0.150] GeV/c^2")
     .note(:background_veto,
            "pi+ pi- pairs with invariant mass in [0.468, 0.528] GeV/c^2 are vetoed to suppress K_S0 contamination")

# --- kinematic fit : 4-momentum conservation + both tag and signal D masses
#     constrained to the nominal D+ mass, chi2 < 200 ---
alg_dp.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"D+")   # tag D-    -> M(D+)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:"D+")   # signal D+ -> M(D+)
  f.chi2_cut 200
end

alg_dp.apply                                        # no Selection argument for tag analyses
root_files_dp = alg_dp.execute_on([data_3773, incMC_3773])