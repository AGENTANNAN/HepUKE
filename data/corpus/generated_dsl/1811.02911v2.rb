# ===== Dataset preparation =====
# Real 4.178 GeV data (3.19 fb^-1) and its inclusive MC (BOSS 703)
data_4180  = DatasetManager.real_data.find("703_4180")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")

# Decay card: signal mode I — D_s+ -> K0 e+ nu_e, K0 -> K_S0 -> pi+ pi-
# Production: e+e- -> D_s*+ D_s- with D_s*+ -> gamma D_s+  (psi(4260) = KKMC top mother)
decay_card_K0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s*+ D_s-  PHSP;
  Enddecay

  Decay D_s*+
  1.0000 gamma D_s+  PHSP;
  Enddecay

  Decay D_s+
  1.0000 K_S0 e+ nu_e  PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-  PHSP;
  Enddecay

  End
DECAYCARD

# Decay card: signal mode II — D_s+ -> K*0 e+ nu_e, K*0 -> K+ pi-
decay_card_Kstar0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s*+ D_s-  PHSP;
  Enddecay

  Decay D_s*+
  1.0000 gamma D_s+  PHSP;
  Enddecay

  Decay D_s+
  1.0000 K*0 e+ nu_e  PHSP;
  Enddecay

  Decay K*0
  1.0000 K+ pi-  VSS;
  Enddecay

  End
DECAYCARD

# 2M-event exclusive MC samples, one per signal mode, associated with the 4.178 GeV data
exMC_K0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "ds_semilep_K0_enu"
  config.related_dataset = data_4180
  config.events          = 2_000_000
  config.decay_card      = decay_card_K0
  config.cross_section   = :default
end

exMC_Kstar0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "ds_semilep_Kstar0_enu"
  config.related_dataset = data_4180
  config.events          = 2_000_000
  config.decay_card      = decay_card_Kstar0
  config.cross_section   = :default
end

# ===== Mode I: D_s+ -> K0 (-> K_S0 -> pi+ pi-) e+ nu_e =====
alg_K0 = TagAnalysis.new("DsTagK0Enu")
alg_K0.set_header(["DsTagK0EnuAlg/DsTagK0Enu.h"])
      .set_constant({ "ECMS" => [:double, 4.178] })
      .with_decay_card(decay_card_K0)

# Tag side: single D_s- tag in 13 hadronic ST modes (tag pinned to the anti-particle)
alg_K0.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKsKPi0, :DstoKKPiPi0, :DstoKsKPiPi,
          :DstoKsKPiPiPi0, :DstoPiPiPi, :DstoPiPiPiPi0, :DstoPiEta,
          :DstoPiPi0Eta, :DstoPiPiPiEta, :DstoKPiPi, :DstoKPiPiPi0
  t.charm -1
end

# Signal side: one pi+ (from K_S0), one pi- (from K_S0), one e+; net charge +1; neutrino missing
alg_K0.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

# 4C kinematic constraint (tag + signal + missing nu = CMS four-momentum), chi2 < 200
alg_K0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_K0.apply
alg_K0.execute_on([data_4180, incMC_4180, exMC_K0])

# ===== Mode II: D_s+ -> K*0 (-> K+ pi-) e+ nu_e =====
alg_Kstar0 = TagAnalysis.new("DsTagKstar0Enu")
alg_Kstar0.set_header(["DsTagKstar0EnuAlg/DsTagKstar0Enu.h"])
         .set_constant({ "ECMS" => [:double, 4.178] })
         .with_decay_card(decay_card_Kstar0)

# Same tag side and tag modes as mode I
alg_Kstar0.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKsKPi0, :DstoKKPiPi0, :DstoKsKPiPi,
          :DstoKsKPiPiPi0, :DstoPiPiPi, :DstoPiPiPiPi0, :DstoPiEta,
          :DstoPiPi0Eta, :DstoPiPiPiEta, :DstoKPiPi, :DstoKPiPiPi0
  t.charm -1
end

# Signal side: one K+ (from K*0), one pi- (from K*0), one e+; net charge +1; neutrino missing
alg_Kstar0.signal_side do |s|
  s.charged(kp: 1, pim: 1, ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

# Same 4C kinematic constraint, chi2 < 200
alg_Kstar0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_Kstar0.apply
alg_Kstar0.execute_on([data_4180, incMC_4180, exMC_Kstar0])