# Paper: 1811.00392v1 — Observation of D_s+ → ω π+ and evidence for D_s+ → ω K+
# BESIII, √s = 4.178 GeV, 3.19 fb⁻¹
# Tag-based analysis: D_s tag ST + fully reconstructed signal side

### Dataset preparation ###
data_4180 = DatasetManager.real_data.find("703_4180")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")

# Common decay card base for D_s*+ D_s- production
decay_card_omega_pi = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s- PHSP;
  Enddecay

  Decay D_s*+
  1.000 D_s+ gamma VSP_PWAVE;
  Enddecay

  Decay D_s+
  1.000 omega pi+ PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

decay_card_omega_K = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s- PHSP;
  Enddecay

  Decay D_s*+
  1.000 D_s+ gamma VSP_PWAVE;
  Enddecay

  Decay D_s+
  1.000 omega K+ PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC samples
exMC_omega_pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Ds_omega_pi"
  config.related_dataset = data_4180
  config.events = 2_000_000
  config.decay_card = decay_card_omega_pi
  config.cross_section = :default
end

exMC_omega_K = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Ds_omega_K"
  config.related_dataset = data_4180
  config.events = 2_000_000
  config.decay_card = decay_card_omega_K
  config.cross_section = :default
end

# ==============================================================
# Mode I: D_s+ → ω π+ (signal side: ω→π+π-π0, π+)
# ==============================================================
alg_omega_pi = TagAnalysis.new("DsOmegaPi")
alg_omega_pi.set_header(["DsOmegaPiAlg/DsOmegaPi.h"])
  .set_constant({ "ECMS" => [:double, 4.178] })
  .with_decay_card(decay_card_omega_pi)

alg_omega_pi.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi
end

alg_omega_pi.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0
  s.charged(pip: 2, pim: 1)
  s.require_charge 1
end

alg_omega_pi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_omega_pi.apply
alg_omega_pi.execute_on([data_4180, incMC_4180, exMC_omega_pi])

# ==============================================================
# Mode II: D_s+ → ω K+ (signal side: ω→π+π-π0, K+)
# ==============================================================
alg_omega_K = TagAnalysis.new("DsOmegaK")
alg_omega_K.set_header(["DsOmegaKAlg/DsOmegaK.h"])
  .set_constant({ "ECMS" => [:double, 4.178] })
  .with_decay_card(decay_card_omega_K)
  .note(:background_veto, "K_S0 veto: events with |M(π+π-) - m_K_S0| < 0.03 GeV/c² and L/σ_L > 2 are vetoed to suppress D_s+ → K_S0 K+ π0 background")
  .note(:tag_mode_unavailable, "omega mass constraint and pi0 mass constraint applied at ROOT level via 2D fit to M(π+π-π0) and M(D_s)")

alg_omega_K.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi
end

alg_omega_K.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0
  s.charged(kp: 1, pip: 1, pim: 1)
  s.require_charge 1
end

alg_omega_K.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_omega_K.apply
alg_omega_K.execute_on([data_4180, incMC_4180, exMC_omega_K])