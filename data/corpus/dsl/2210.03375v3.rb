# 2210.03375v3: Observations of Cabibbo-Suppressed decays Λc+ → nπ+π0, nπ+π-π+, nK-π+π+
# Double-tag (DT) method with ST Λc- (12 hadronic tag modes)
# Signal: Λc+ decay in system recoiling against ST Λc-
# Neutron reconstructed via missing-mass technique
# 7 energy points: 4599.53-4698.82 MeV

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# Data at 7 c.m. energies
data_4600 = DatasetManager.real_data.find("703_4600")
data_4612 = DatasetManager.real_data.find("706_4610")
data_4628 = DatasetManager.real_data.find("706_4620")
data_4641 = DatasetManager.real_data.find("706_4640")
data_4661 = DatasetManager.real_data.find("706_4660")
data_4682 = DatasetManager.real_data.find("706_4680")
data_4699 = DatasetManager.real_data.find("706_4700")

all_data = [data_4600, data_4612, data_4628, data_4641, data_4661, data_4682, data_4699]

# ======================================================================
# SIGNAL MODE 1: Λc+ → nπ+π0
# Signal side: 1 π+ + 1 π0, missing neutron
# ST: Λc- with 10 tag modes (excluding pK+π-π0 and pπ-π+)
# ======================================================================
alg_mode1 = TagAnalysis.new("LcToNPiPi0")
alg_mode1.set_header(["LcToNPiPi0Alg/LcToNPiPi0.h"])
          .set_constant({ "ECMS" => [:double, 4.640] })

alg_mode1.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoPKs,
          :LambdacPtoKPiP,
          :LambdacPtoPKsPi0,
          :LambdacPtoPKsPiPi,
          :LambdacPtoLambdaPi,
          :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi,
          :LambdacPtoSigmaPi0,
          :LambdacPtoSigmaPi0Pi0,
          :LambdacPtoSigmaPiPi
  t.charm -1
end

alg_mode1.signal_side do |s|
  s.photons 2
  s.charged(pip: 1)
  s.require_charge 1
  s.min_photon_angle 10.0
  s.missing :n0
end

alg_mode1.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mode1.note(:tag_mode_unavailable,
  "ST modes pK+π-π0 and pπ-π+ excluded; see Table 2 of paper")
alg_mode1.note(:peaking_background_veto,
  "Λc+→Λπ+ with Λ→nπ0, Λc+→Σ+π0 with Σ+→nπ+, Λc+→Σ0π+ with Σ0→γΛ(→nπ0) vetoed. " \
  "Mmiss(π+)>1300 MeV/c² and Mmiss(π0)>1370 MeV/c²")

alg_mode1.apply

# ======================================================================
# SIGNAL MODE 2: Λc+ → nπ+π-π+
# Signal side: 2 π+ + 1 π-, missing neutron
# ST: Λc- with 11 tag modes (excluding pπ-π+)
# ======================================================================
alg_mode2 = TagAnalysis.new("LcToNPiPiPi")
alg_mode2.set_header(["LcToNPiPiPiAlg/LcToNPiPiPi.h"])
          .set_constant({ "ECMS" => [:double, 4.640] })

alg_mode2.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoPKs,
          :LambdacPtoKPiP,
          :LambdacPtoPKsPi0,
          :LambdacPtoPKsPiPi,
          :LambdacPtoKPiPPi0,
          :LambdacPtoLambdaPi,
          :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi,
          :LambdacPtoSigmaPi0,
          :LambdacPtoSigmaPi0Pi0,
          :LambdacPtoSigmaPiPi
  t.charm -1
end

alg_mode2.signal_side do |s|
  s.charged(pip: 2, pim: 1)
  s.require_charge 1
  s.missing :n0
end

alg_mode2.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mode2.note(:tag_mode_unavailable,
  "ST mode pπ-π+ excluded; see Table 2 of paper")
alg_mode2.note(:peaking_background_veto,
  "KS0→π+π- veto (Mπ+π- ∉ [487,511] MeV/c²). " \
  "Σ+→nπ+ veto (Mmiss(π+π-) ∉ [1150,1250] MeV/c²). " \
  "Σ-→nπ- veto (Mmiss(π+π+) ∉ [1150,1250] MeV/c²)")

alg_mode2.apply

# ======================================================================
# SIGNAL MODE 3: Λc+ → nK-π+π+
# Signal side: 1 K- + 2 π+, missing neutron
# ST: Λc- with 12 tag modes (all)
# ======================================================================
alg_mode3 = TagAnalysis.new("LcToNKPiPi")
alg_mode3.set_header(["LcToNKPiPiAlg/LcToNKPiPi.h"])
          .set_constant({ "ECMS" => [:double, 4.640] })

alg_mode3.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoPKs,
          :LambdacPtoKPiP,
          :LambdacPtoPKsPi0,
          :LambdacPtoPKsPiPi,
          :LambdacPtoKPiPPi0,
          :LambdacPtoLambdaPi,
          :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi,
          :LambdacPtoSigmaPi0,
          :LambdacPtoSigmaPi0Pi0,
          :LambdacPtoSigmaPiPi,
          :LambdacPtoPPiPi
  t.charm -1
end

alg_mode3.signal_side do |s|
  s.charged(km: 1, pip: 2)
  s.require_charge 1
  s.missing :n0
end

alg_mode3.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mode3.apply

# Execute all three modes
[alg_mode1, alg_mode2, alg_mode3].each do |alg|
  alg.execute_on(all_data)
end

# Note: MBC signal region 2.275-2.300 GeV/c² at 4.600 GeV,
# 2.275-2.306 GeV/c² at 4.612-4.641 GeV,
# 2.275-2.310 GeV/c² at 4.661-4.699 GeV
# Tag ΔE window: mode-dependent (~3σ), see Table 2