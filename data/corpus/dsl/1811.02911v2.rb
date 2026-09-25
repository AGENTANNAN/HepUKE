# Paper: 1811.02911v2 — Semileptonic D_s+ decays: D_s+ → K0 e+ ν_e and D_s+ → K*0 e+ ν_e
# BESIII, √s = 4.178 GeV, 3.19 fb⁻¹
# Tag-based analysis: D_s ST + missing (ν_e)

### Dataset preparation ###
data_4180 = DatasetManager.real_data.find("703_4180")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")

# ==============================================================
# Mode I: D_s+ → K0 e+ ν_e (K0 → K_S0 → π+π-)
# ==============================================================
decay_card_K0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s- PHSP;
  Enddecay

  Decay D_s*+
  1.000 D_s+ gamma VSP_PWAVE;
  Enddecay

  Decay D_s+
  1.000 K0 e+ nu_e ISGW2;
  Enddecay

  Decay K0
  1.000 K_S0 PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

exMC_K0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Ds_K0_enu"
  config.related_dataset = data_4180
  config.events = 2_000_000
  config.decay_card = decay_card_K0
  config.cross_section = :default
end

alg_K0 = TagAnalysis.new("DsK0Enu")
alg_K0.set_header(["DsK0EnuAlg/DsK0Enu.h"])
  .set_constant({ "ECMS" => [:double, 4.178] })
  .with_decay_card(decay_card_K0)

# 13 ST tag modes
alg_K0.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKsKPi0, :DstoKKPiPi0,
          :DstoKsKminusPiPi, :DstoPiPiPi, :DstoPiEta,
          :DstoPiPi0Eta, :DstoPiEPPiPiEta, :DstoPiEPRhoGam
end

alg_K0.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

alg_K0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_K0.apply
alg_K0.execute_on([data_4180, incMC_4180, exMC_K0])

# ==============================================================
# Mode II: D_s+ → K*0 e+ ν_e (K*0 → K+ π-)
# ==============================================================
decay_card_Kstar = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s- PHSP;
  Enddecay

  Decay D_s*+
  1.000 D_s+ gamma VSP_PWAVE;
  Enddecay

  Decay D_s+
  1.000 K*0 e+ nu_e ISGW2;
  Enddecay

  Decay K*0
  1.000 K+ pi- VSS;
  Enddecay

  End
DECAYCARD

exMC_Kstar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Ds_Kstar_enu"
  config.related_dataset = data_4180
  config.events = 2_000_000
  config.decay_card = decay_card_Kstar
  config.cross_section = :default
end

alg_Kstar = TagAnalysis.new("DsKstarEnu")
alg_Kstar.set_header(["DsKstarEnuAlg/DsKstarEnu.h"])
  .set_constant({ "ECMS" => [:double, 4.178] })
  .with_decay_card(decay_card_Kstar)

alg_Kstar.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKsKPi0, :DstoKKPiPi0,
          :DstoKsKminusPiPi, :DstoPiPiPi, :DstoPiEta,
          :DstoPiPi0Eta, :DstoPiEPPiPiEta, :DstoPiEPRhoGam
end

alg_Kstar.signal_side do |s|
  s.charged(kp: 1, pim: 1, ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

alg_Kstar.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_Kstar.apply
alg_Kstar.execute_on([data_4180, incMC_4180, exMC_Kstar])