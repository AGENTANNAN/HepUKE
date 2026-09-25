# DSL for arxiv:1811.10890v2
# Ds+ -> mu+ nu_mu at sqrt(s) = 4.178 GeV (BOSS 703, 4180)
# Tag-based: ST Ds- with 13 unique hadronic tag modes, DT signal: gamma/pi0 + mu+ + nu_mu

decay_card_ds_munu = <<~DECAYCARD
  Decay psi(4260)
  1 D_s+ D_s*- PHSP;
  Enddecay
  Decay D_s+
  1 mu+ nu_mu PHSP;
  Enddecay
  Decay D_s*-
  1 gamma D_s- PHSP;
  Enddecay
  End
DECAYCARD

data_4180  = DatasetManager.real_data.find("703_4180")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")

sig_ds_munu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_ds_munu"
  config.related_dataset = data_4180
  config.events          = 100_000
  config.decay_card      = decay_card_ds_munu
  config.cross_section   = :default
end

alg = TagAnalysis.new("DsTagMuNu")
alg.set_header(["DsTagMuNuAlg/DsTagMuNu.h"])
  .set_constant({ "ECMS" => [:double, 4.178] })
  .with_decay_card(decay_card_ds_munu)

# ST tag side: Ds- reconstructed in 13 hadronic modes
# Note: eta_{gamma gamma} pi- and eta_{pi0 pi+ pi-} pi- are the same DTagAlg mode
# (Ds -> pi eta), as DTagAlg handles the eta decay internally
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKKPiPi0, :DstoKsK, :DstoKsKPi0,
          :DstoKsKsPi, :DstoKsKplusPiPi, :DstoKsKminusPiPi,
          :DstoKPiPi, :DstoPiPiPi, :DstoPiEta,
          :DstoPiEPPiPiEta, :DstoPiEPRhoGam, :DstoPiPi0Eta
  t.charm -1
end

# Signal side: remaining tracks -> mu+ + missing nu_mu
# The soft gamma/pi0 from Ds* -> gamma/pi0 Ds is handled by DTagAlg internally
alg.signal_side do |s|
  s.charged(mup: 1)
  s.require_charge 1
  s.missing :nu_mu              # massless neutrino
end

# Kinematic fit: 4C + Ds mass constraints
# The tag-side gamma/pi0 from Ds* and Ds mass constraints are internal to DTagAlg
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.note(:ds_star_gamma_pi0, "Soft gamma/pi0 from Ds* -> gamma/pi0 Ds selected by minimum |DeltaE|; Ds* mass constrained in kinematic fit internally")
alg.note(:mm2_fit, "DT yield extracted from unbinned constrained fit to MM^2 distribution — ROOT-level analysis")
alg.note(:muon_pid, "Muon PID: depositted energy in EMC within (0.0, 0.3) GeV, 2D hit depth requirement in muon counter")
alg.note(:extra_photon_veto, "Maximum energy of unused showers < 0.4 GeV; no additional good charged tracks allowed")
alg.note(:tag_deltae_window, "Paper uses M_BC window (2.010, 2.073) GeV/c^2 for ST selection — applied in ROOT, not here (store-not-cut)")
alg.note(:eta_sub_decays, "eta_{gamma gamma}pi and eta_{pi0 pi+ pi-}pi are both DstoPiEta in DTagAlg; the sub-decay selection is internal to DTagTool")
alg.note(:signal_deltae, "DeltaE in (-0.05, 0.10) GeV for DT selection")

alg.apply
alg.execute_on([data_4180, incMC_4180, sig_ds_munu])