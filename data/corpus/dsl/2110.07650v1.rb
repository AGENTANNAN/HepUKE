# DSL for arXiv:2110.07650v1
# Amplitude analysis of Ds+ -> KS0 KS0 pi+ at sqrt(s)=4.178-4.226 GeV
# TagAnalysis (DT): D-tag technique with 8 hadronic tag modes for Ds-
# Multi-energy: 4.178, 4.189, 4.199, 4.209, 4.219, 4.226 GeV

data_4180 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")

incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

all_data = [data_4180, data_4190, data_4200, data_4210, data_4220, data_4230]
all_incMC = [incMC_4180, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4230]

alg = TagAnalysis.new("DsToKsKsPi")
alg.set_header(["DsToKsKsPiAlg/DsToKsKsPi.h"])
   .set_constant({"ECMS" => [:double, 4.178]})
   .note(:multi_energy, "Analysis uses 6 energy points: 4.178, 4.189, 4.199, 4.209, 4.219, 4.226 GeV")
   .note(:amplitude_analysis, "Amplitude analysis (isobar model) performed in ROOT; DSL covers event selection only")
   .note(:bf_measurement, "BF measured via DT method without kinematic fit (separate selection from amplitude analysis)")
   .note(:signal_mode_invented, "Signal mode Ds+ -> KS0 KS0 pi+ not in authoritative DTagAlg mode list; requires custom DTagAlg channel")

# Tag side: Ds- reconstructed via 8 hadronic DTagAlg modes
# Paper modes: KS0K-, K+K-pi-, K+K-pi-pi0, KS0K-pi-pi+, KS0K+pi-pi-, pi-pi-pi+, pi-eta', K-pi+pi-
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKKPiPi0, :DstoKsKPiPi,
          :DstoPiPiPi, :DstoPiEta, :DstoKPiPi
  t.charm -1
end

# Note: :DstoKsKPiPi covers two modes (KS0K-pi-pi+ and KS0K+pi-pi-)
# :DstoPiEta used as proxy for pi-eta' mode
alg.note(:mode_coverage, "8 tag modes mapped to 7 DTagAlg symbols; DstoKsKPiPi covers 2 modes, DstoPiEta proxies for pi-eta'")

# Signal side: Ds+ -> KS0 KS0 pi+ (reconstructed via DTagAlg with custom mode)
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsKPiPi   # proxy for Ds+ -> KS0 KS0 pi+ (custom mode, not in authoritative list)
  t.charm 1
end

# 8C kinematic fit: 4C + 2x KS0 mass + tag Ds mass + Ds* mass constraints
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:Ds)
  f.chi2_cut 200
end

alg.note(:fit_8c, "Full 8C fit: 4C + 2x KS0 mass + tag Ds + signal Ds + Ds* mass constraints. DSL covers 4C + Ds mass constraints. KS0 mass and Ds* mass constraints applied in BOSS jobOptions.")
alg.note(:signal_mass_window, "Signal Ds+ invariant mass window [1.950, 1.990] GeV/c^2 applied in ROOT")
alg.note(:bf_no_kinfit, "BF measurement uses DT selection without kinematic fit, requiring soft pion momentum > 0.1 GeV/c to suppress D*+ decays")

alg.dtag_reconstruction do |d|
  d.beam_energy :db
  d.local true
end

alg.apply
alg.execute_on(all_data + all_incMC)