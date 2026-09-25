# BESIII DSL: Observation of D -> K1(1270) mu+ nu_mu and test of LFU
# Paper: 2502.03828v2
# CMS energy: 3.773 GeV (psi(3770)), luminosity 7.93 fb^-1
# Method: Single-tag + missing neutrino (semileptonic)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================
# Decay cards for exclusive signal MC
# ============================================================

decay_card_dp_signal = <<~DECAYCARD
  Decay psi(3770)
  1.000 D+ D- PHSP;
  Enddecay
  Decay D+
  1.000 K- pi+ pi0 mu+ nu_mu PHSP;
  Enddecay
  Decay D-
  1.000 K+ pi- pi- PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

decay_card_d0_signal = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.000 K- pi+ pi- mu+ nu_mu PHSP;
  Enddecay
  Decay anti-D0
  1.000 K+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# ============================================================
# Exclusive MC samples
# ============================================================

sig_mc_dp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Dp_K1muNu"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_dp_signal
  config.cross_section   = :default
end

sig_mc_d0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_D0_K1muNu"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_d0_signal
  config.cross_section   = :default
end

# ============================================================
# Algorithm 1: D+ -> K1(1270)0 mu+ nu_mu
# ST + missing: tag D- (anti-D+), signal D+ -> K- pi+ pi0 mu+ nu_mu
# K1(1270)0 -> K- pi+ pi0
# ============================================================

alg_dp = TagAnalysis.new("DpK1muNu")
alg_dp.set_header(["DpK1muNuAlg/DpK1muNu.h"])
     .set_constant({ "ECMS" => [:double, 3.773] })
     .note(:signal_mc_generator,
       "Signal MC uses a dedicated generator from amplitude analysis of D+ -> K- pi+ pi0 e+ nu")

alg_dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :mBC, min: 1.863, max: 1.877
end

alg_dp.signal_side do |s|
  s.photons 2
  s.charged(km: 1, pip: 1, mup: 1)
  s.require_charge 1
  s.missing :nu_mu
end

alg_dp.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_dp.note(:per_mode_deltaE,
  "Per-mode deltaE windows: DtoKPiPi (-0.025,0.024), DtoKsPi (-0.025,0.026), " \
  "DtoKPiPiPi0 (-0.057,0.046), DtoKsPiPi0 (-0.062,0.049), " \
  "DtoKsPiPiPi (-0.028,0.027), DtoKKPi (-0.024,0.023)")
alg_dp.note(:extra_shower_veto,
  "N_extra_pi0 = 0; E_extra_gamma_max < 0.20 GeV; cos_theta_pmiss_gamma < 0.69")
alg_dp.note(:signal_mass_window,
  "M(K pi pi0 mu) < 1.72 GeV/c^2; K1 mass window (1.163, 1.343) GeV/c^2")
alg_dp.note(:muon_pid_cuts,
  "CL_mu > 0.001, CL_mu > CL_K, CL_mu > CL_pi; EMC deposit < 0.28 GeV")
alg_dp.note(:background_veto,
  "Peaking background from D+ -> K- pi+ pi+ pi0 suppressed via mass and PID cuts")

alg_dp.apply
alg_dp.execute_on([data_3773, incMC_3773, sig_mc_dp])

# ============================================================
# Algorithm 2: D0 -> K1(1270)- mu+ nu_mu
# ST + missing: tag anti-D0, signal D0 -> K- pi+ pi- mu+ nu_mu
# K1(1270)- -> K- pi+ pi-
# ============================================================

alg_d0 = TagAnalysis.new("D0K1muNu")
alg_d0.set_header(["D0K1muNuAlg/D0K1muNu.h"])
     .set_constant({ "ECMS" => [:double, 3.773] })
     .note(:signal_mc_generator,
       "Signal MC uses a dedicated generator from amplitude analysis of D0 -> K- pi+ pi- e+ nu")

alg_d0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :mBC, min: 1.859, max: 1.873
end

alg_d0.signal_side do |s|
  s.charged(km: 1, pip: 1, pim: 1, mup: 1)
  s.require_charge 0
  s.missing :nu_mu
end

alg_d0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0.note(:per_mode_deltaE,
  "Per-mode deltaE windows: D0toKPi (-0.027,0.027), D0toKPiPi0 (-0.062,0.049), " \
  "D0toKPiPiPi (-0.026,0.024)")
alg_d0.note(:extra_shower_veto,
  "N_extra_pi0 = 0; E_extra_gamma_max < 0.20 GeV; cos_theta_pmiss_gamma < 0.54")
alg_d0.note(:signal_mass_window,
  "M(K pi pi mu) < 1.65 GeV/c^2; K1 mass window (1.163, 1.343) GeV/c^2")
alg_d0.note(:muon_pid_cuts,
  "CL_mu > 0.001, CL_mu > CL_K, CL_mu > CL_pi; EMC deposit < 0.28 GeV")
alg_d0.note(:background_veto,
  "Peaking backgrounds from D0 -> K- pi- eta'(pi+ pi- gamma), " \
  "D0 -> K- pi+ K+ pi-, D0 -> K- pi+ pi- pi+ pi0 suppressed via mass and PID cuts")

alg_d0.apply
alg_d0.execute_on([data_3773, incMC_3773, sig_mc_d0])