# Search for FCNC J/ψ → D0 μ+ μ- with three D0 tag modes
# J/ψ, 1.0087×10^10 events, arXiv:2501.08080v2
# Tag-based analysis: D0 reconstructed from pre-stored DTagAlg candidates

jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card = <<~DECAYCARD
  Decay J/psi
  1.000 D0 mu+ mu- PHSP;
  Enddecay

  Decay D0
  1.000 K- pi+ PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_D0_mumu_fcnc"
  config.related_dataset = jpsi_data
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

# TagAnalysis: tag D0 from pre-stored candidates, signal side μ+μ-
alg = TagAnalysis.new("JpsiD0MuMuFCNC")
alg.set_header(["JpsiD0MuMuFCNCAlg/JpsiD0MuMuFCNC.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })

# Three D0 tag modes: Kπ (Mode I), Kππ0 (Mode II), Kπππ (Mode III)
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
end

# Signal side: μ+ μ- from FCNC decay
alg.signal_side do |s|
  s.charged(mup: 1, mum: 1)
  s.require_charge 0
end

# 5C kinematic fit: 4C + D0 mass constraint
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg.apply

# Additional cuts applied in ROOT:
# - D0 mass windows per mode: Mode I [1.84,1.89], Mode II [1.80,1.91], Mode III [1.84,1.89] GeV/c²
# - χ²_D0 cuts per mode: 4.5, 3.1, 3.7
# - Muon PID: E_μ in [0.11, 0.25] GeV (EMC energy)
# - |p_miss| < 0.05 GeV/c
# - K_S0 veto, M_4π veto, M_recoil vetoes
# - Signal region: M(D0μ+μ-) in [3.05, 3.15] GeV/c²
alg.note(:root_cuts, "ROOT-level: D0 mass windows per mode (I: [1.84,1.89], II: [1.80,1.91], III: [1.84,1.89] GeV/c²), χ²_D0 cuts (4.5, 3.1, 3.7), E_μ in [0.11,0.25] GeV, |p_miss|<0.05 GeV/c, K_S0/M_4π/M_recoil vetoes, signal region M(D0μ+μ-) in [3.05,3.15] GeV/c²")
alg.note(:signal_extraction, "Signal yield from fit to M(D0μ+μ-) distribution; upper limits set via Bayesian method")

alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])