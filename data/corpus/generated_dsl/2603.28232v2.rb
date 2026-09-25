# frozen_string_literal: true

### Dataset preparation ###
# Real data at 4.600-4.843 GeV (11 energy points)
data_points = [
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("706_4610"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4740"),
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
  DatasetManager.real_data.find("707_4840")
]

# Matching inclusive MC samples (available for ten of the energy points)
incMC_points = [
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
  DatasetManager.inclusive_mc.find("707_4740"),
  DatasetManager.inclusive_mc.find("707_4750"),
  DatasetManager.inclusive_mc.find("707_4780"),
  DatasetManager.inclusive_mc.find("707_4840")
]

# Decay card (EvtGen):
#   e+e- -> Lambda_c+ anti-Lambda_c-
#   Lambda_c+ -> X pi+ eta with X = n, Lambda0, Sigma0 in equal parts
#   anti-Lambda_c- -> anti-p K+ pi-
#   eta -> gamma gamma / pi+ pi- pi0 / pi0 pi0 pi0
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.00000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  0.33333 n0       pi+ eta PHSP;
  0.33333 Lambda0  pi+ eta PHSP;
  0.33334 Sigma0   pi+ eta PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.00000 anti-p- K+ pi- PHSP;
  Enddecay

  Decay eta
  0.39400 gamma gamma PHSP;
  0.22900 pi+ pi- pi0 PHSP;
  0.37700 pi0 pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.00000 gamma gamma PHSP;
  Enddecay

  Decay Lambda0
  1.00000 p+ pi- HypWK;
  Enddecay

  Decay Sigma0
  1.00000 Lambda0 gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive signal MC: 500k events per energy point (same decay card for all points)
exMC_signals = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lambdac_npiEta"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection - tag-based analysis ###
# Single tag: a hadronically reconstructed anti-Lambda_c- tags the event; the signal side is
# Lambda_c+ -> n pi+ eta. Two independent signal-side chains are used:
#   chain 1: eta -> gamma gamma
#   chain 2: eta -> pi+ pi- pi0

# ---- Chain 1: eta -> gamma gamma ----
alg_gg = TagAnalysis.new("LambdaCtagEtaToGG")
alg_gg.set_header(["LambdaCtagEtaToGGAlg/LambdaCtagEtaToGG.h"])
      .set_constant({"ECMS" => [:double, 4.720]})   # per-run measured beam energy used via MeasuredEcmsSvc
      .with_decay_card(decay_card_signal)

alg_gg.tag_side(:Lambdac) do |t|
  t.mode_group :hadronic   # 11 hadronic anti-Lambda_c- tag modes
  t.charm -1               # pin the tagged side to anti-Lambda_c-
end

alg_gg.signal_side do |s|
  s.photons 2              # two photons from eta -> gamma gamma
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.charged(pip: 1)        # one positive signal track
  s.require_charge 1       # net signal-side charge +1
  s.missing :n             # missing neutron
end

alg_gg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)   # gamma gamma -> eta
  f.chi2_cut 200
end

alg_gg.note(:efficiency_curve, "M(gamma gamma) required in [0.505, 0.575] GeV/c^2 after a standalone 1C eta-mass fit with chi2 < 20")
alg_gg.note(:background_veto, "M_recoil(eta) in (1.13, 1.25) GeV/c^2 vetoed to suppress the Lambda_c+ -> Sigma+ eta background")
alg_gg.note(:missing_mass_observable, "M_miss = sqrt(E_miss^2 - |p_miss|^2) of the missing neutron is the extraction observable (auto-stored)")

alg_gg.apply
alg_gg.execute_on(data_points + incMC_points + exMC_signals)

# ---- Chain 2: eta -> pi+ pi- pi0 ----
alg_pipimpi0 = TagAnalysis.new("LambdaCtagEtaToPiPiPi0")
alg_pipimpi0.set_header(["LambdaCtagEtaToPiPiPi0Alg/LambdaCtagEtaToPiPiPi0.h"])
            .set_constant({"ECMS" => [:double, 4.720]})
            .with_decay_card(decay_card_signal)

alg_pipimpi0.tag_side(:Lambdac) do |t|
  t.mode_group :hadronic   # 11 hadronic anti-Lambda_c- tag modes
  t.charm -1               # pin the tagged side to anti-Lambda_c-
end

alg_pipimpi0.signal_side do |s|
  s.photons 2                    # two photons from pi0 -> gamma gamma
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.charged(pip: 2, pim: 1)      # two pi+ and one pi-
  s.require_charge 1             # net signal-side charge +1
  s.missing :n                   # missing neutron
end

alg_pipimpi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # gamma gamma -> pi0
  f.chi2_cut 200
end

alg_pipimpi0.note(:efficiency_curve, "M(gamma gamma) in [0.115, 0.150] GeV/c^2 (1C pi0-mass fit chi2 < 200) and M(pi+ pi- pi0) in [0.505, 0.575] GeV/c^2 (1C eta-mass fit chi2 < 20) required; the pi+ pi- pi0 -> eta mass constraint is not expressible in the tag fit")
alg_pipimpi0.note(:background_veto, "M_recoil(eta) in (1.16, 1.22) GeV/c^2 and M_recoil(pi+ pi- pi+) in (1.10, 1.13) GeV/c^2 vetoed")
alg_pipimpi0.note(:missing_mass_observable, "M_miss of the missing neutron is the extraction observable (auto-stored)")

alg_pipimpi0.apply
alg_pipimpi0.execute_on(data_points + incMC_points + exMC_signals)