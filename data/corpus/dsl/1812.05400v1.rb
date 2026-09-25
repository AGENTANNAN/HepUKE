# DSL for arxiv:1812.05400v1
# D+ -> K_{S,L}^0 K+ (pi0) at sqrt(s) = 3.773 GeV (BOSS 712, 3773)
# Tag-based: ST D- in 6 hadronic modes, DT signal: 4 independent channels
# Signal modes: K_S0 K+, K_S0 K+ pi0, K_L0 K+, K_L0 K+ pi0

data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# =====================
# Shared: decay cards
# =====================

decay_card_ks_k  = <<~DECAYCARD
  Decay psi(3770)
  1 D+ D- PHSP;
  Enddecay
  Decay D+
  1 K_S0 K+ PHSP;
  Enddecay
  End
DECAYCARD

decay_card_ks_kpi0 = <<~DECAYCARD
  Decay psi(3770)
  1 D+ D- PHSP;
  Enddecay
  Decay D+
  1 K_S0 K+ pi0 PHSP;
  Enddecay
  End
DECAYCARD

decay_card_kl_k = <<~DECAYCARD
  Decay psi(3770)
  1 D+ D- PHSP;
  Enddecay
  Decay D+
  1 K_L0 K+ PHSP;
  Enddecay
  End
DECAYCARD

decay_card_kl_kpi0 = <<~DECAYCARD
  Decay psi(3770)
  1 D+ D- PHSP;
  Enddecay
  Decay D+
  1 K_L0 K+ pi0 PHSP;
  Enddecay
  End
DECAYCARD

# Shared 6 ST D- tag modes
# D± → K∓π±π± → DptoKPiPi, D± → K∓π±π±π0 → DptoKPiPiPi0
# D± → K_S^0π∓ → DptoKsPi, D± → K_S^0π∓π0 → DptoKsPiPi0
# D± → K_S^0π±π∓π∓ → DptoKsPiPiPi, D± → K∓K±π∓ → DptoKKPi
ST_TAG_MODES = [:DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0,
                :DptoKsPiPiPi, :DptoKKPi]

# =====================
# Signal Mode 1: D+ -> K_S0 K+
# =====================

sig_ks_k = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_dp_ks0kp"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_ks_k
  config.cross_section   = :default
end

alg_ks_k = TagAnalysis.new("DpTagKs0Kp")
alg_ks_k.set_header(["DpTagKs0KpAlg/DpTagKs0Kp.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })
  .with_decay_card(decay_card_ks_k)

alg_ks_k.tag_side(:Dplus) do |t|
  t.modes(*ST_TAG_MODES)
  t.charm -1     # D- tag
end

alg_ks_k.signal_side do |s|
  s.charged(pip: 1, pim: 1, kp: 1)   # K_S0 -> pi+ pi-, plus signal K+
  s.require_charge 1                   # K_S0(0) + K+(+1) = +1
end

alg_ks_k.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_ks_k.note(:ks_reco, "K_S0 -> pi+ pi- reconstructed via secondary vertex fit. M(pi+pi-) within ±12 MeV of nominal K_S0 mass. Decay length significance > 2 sigma.")
alg_ks_k.note(:dt_2d_fit, "DT yield from unbinned 2D fit on M_BC_tag vs M_BC_sig — ROOT-level analysis")
alg_ks_k.note(:st_deltae, "Tag-mode-dependent DeltaE requirements on ST side (see Table I in paper)")
alg_ks_k.note(:ks_efficiency_correction, "K_S0 reconstruction efficiency corrected for data-MC differences (~2% correction)")

alg_ks_k.apply
alg_ks_k.execute_on([data_3773, incMC_3773, sig_ks_k])

# =====================
# Signal Mode 2: D+ -> K_S0 K+ pi0
# =====================

sig_ks_kpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_dp_ks0kppi0"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_ks_kpi0
  config.cross_section   = :default
end

alg_ks_kpi0 = TagAnalysis.new("DpTagKs0KpPi0")
alg_ks_kpi0.set_header(["DpTagKs0KpPi0Alg/DpTagKs0KpPi0.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })
  .with_decay_card(decay_card_ks_kpi0)

alg_ks_kpi0.tag_side(:Dplus) do |t|
  t.modes(*ST_TAG_MODES)
  t.charm -1
end

alg_ks_kpi0.signal_side do |s|
  s.charged(pip: 1, pim: 1, kp: 1)   # K_S0 -> pi+ pi-, plus signal K+
  s.photons 2                          # pi0 -> gamma gamma
  s.min_photon_energy 0.025
  s.min_photon_angle 10.0
  s.require_charge 1
end

alg_ks_kpi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_ks_kpi0.note(:pi0_sel, "pi0 -> gamma gamma: M(gamma gamma) in (0.110, 0.155) GeV/c^2, chi2 from kinematic fit to nominal pi0 mass < 20")
alg_ks_kpi0.note(:ks_reco, "K_S0 reconstruction as above")
alg_ks_kpi0.note(:dt_2d_fit, "DT yield from unbinned 2D fit on M_BC_tag vs M_BC_sig")
alg_ks_kpi0.note(:signal_deltae, "DeltaE_sig in (-0.057, 0.040) GeV for K_S0 K+ pi0 mode")

alg_ks_kpi0.apply
alg_ks_kpi0.execute_on([data_3773, incMC_3773, sig_ks_kpi0])

# =====================
# Signal Mode 3: D+ -> K_L0 K+
# =====================

sig_kl_k = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_dp_kl0kp"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_kl_k
  config.cross_section   = :default
end

alg_kl_k = TagAnalysis.new("DpTagKl0Kp")
alg_kl_k.set_header(["DpTagKl0KpAlg/DpTagKl0Kp.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })
  .with_decay_card(decay_card_kl_k)

alg_kl_k.tag_side(:Dplus) do |t|
  t.modes(*ST_TAG_MODES)
  t.charm -1
end

alg_kl_k.signal_side do |s|
  s.charged(kp: 1)                     # signal K+
  s.require_charge 1
  s.missing :K_L0                      # massive missing K_L0 (known mass, direction from EMC shower)
end

alg_kl_k.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_kl_k.note(:kl_shower_sel, "K_L0 direction from EMC shower not used in tag; E_shower > 0.1 GeV. Best candidate selected by minimal chi2 from kinematic fit (mode-dependent chi2 cuts, see Table II). Momentum inferred from DeltaE_sig = 0 constraint.")
alg_kl_k.note(:kl_efficiency_correction, "K_L0 reconstruction efficiency corrected for data-MC differences (~10% correction)")
alg_kl_k.note(:ks_peaking_bkg, "Peaking background D+ -> K_S0 K+ with K_S0 -> pi0 pi0 (~3% of signal) fixed in 2D M_BC fit")
alg_kl_k.note(:dt_2d_fit, "DT yield from unbinned 2D fit on M_BC_tag vs M_BC_sig")

alg_kl_k.apply
alg_kl_k.execute_on([data_3773, incMC_3773, sig_kl_k])

# =====================
# Signal Mode 4: D+ -> K_L0 K+ pi0
# =====================

sig_kl_kpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_dp_kl0kppi0"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_kl_kpi0
  config.cross_section   = :default
end

alg_kl_kpi0 = TagAnalysis.new("DpTagKl0KpPi0")
alg_kl_kpi0.set_header(["DpTagKl0KpPi0Alg/DpTagKl0KpPi0.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })
  .with_decay_card(decay_card_kl_kpi0)

alg_kl_kpi0.tag_side(:Dplus) do |t|
  t.modes(*ST_TAG_MODES)
  t.charm -1
end

alg_kl_kpi0.signal_side do |s|
  s.charged(kp: 1)                     # signal K+
  s.photons 2                          # pi0 -> gamma gamma
  s.min_photon_energy 0.025
  s.min_photon_angle 10.0
  s.require_charge 1
  s.missing :K_L0                      # massive missing K_L0
end

alg_kl_kpi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_kl_kpi0.note(:pi0_sel, "pi0 -> gamma gamma as above")
alg_kl_kpi0.note(:kl_shower_sel, "K_L0 selection as above; mode-dependent chi2 cuts per Table II")
alg_kl_kpi0.note(:ks_peaking_bkg, "Peaking background D+ -> K_S0 K+ pi0 with K_S0 -> pi0 pi0 (~5%) fixed in 2D M_BC fit")
alg_kl_kpi0.note(:dt_2d_fit, "DT yield from unbinned 2D fit on M_BC_tag vs M_BC_sig")

alg_kl_kpi0.apply
alg_kl_kpi0.execute_on([data_3773, incMC_3773, sig_kl_kpi0])