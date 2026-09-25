### Dataset description ###
# psi(3770) at 3.773 GeV: 20.3 fb^-1 of real data + inclusive MC
d3770_data  = DatasetManager.real_data.find("712_3773")
d3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the electron mode: D0 -> anti-K0 pi- e+ nu_e (signal),
# anti-D0 -> hadronic tag modes K+pi-, K+pi-pi+pi-, K+pi-pi0 (pi0 -> gamma gamma)
decay_card_e = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 anti-K0 pi- e+ nu_e PHSP;
    Enddecay

    Decay anti-K0
    1.0000 K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    0.3333 K+ pi- PHSP;
    0.3334 K+ pi- pi+ pi- PHSP;
    0.3333 K+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the muon mode: D0 -> anti-K0 pi- mu+ nu_mu (signal)
decay_card_mu = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 anti-K0 pi- mu+ nu_mu PHSP;
    Enddecay

    Decay anti-K0
    1.0000 K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    0.3333 K+ pi- PHSP;
    0.3334 K+ pi- pi+ pi- PHSP;
    0.3333 K+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 500k events per lepton mode
exMC_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_d0_kspimenu_e"
  config.related_dataset = d3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_e
  config.cross_section   = :default
end

exMC_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_d0_kspimenu_mu"
  config.related_dataset = d3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_mu
  config.cross_section   = :default
end

### Event selection (BOSS) — D0 tag (ST + missing neutrino) analysis ###
# -------------------- Electron mode --------------------
alg_e = TagAnalysis.new("D0KsPimENu")
alg_e.set_header(["D0KsPimENuAlg/D0KsPimENu.h"])
     .set_constant({"ECMS" => [:double, 3.773]})
     .set_alias({"std::vector<double>" => "Vdouble"})
     .note(:deltaE_windows, "mode-dependent tag-side deltaE windows applied in the tag scan (not expressible as a single per-side window declaration): D0->K+pi- [-0.027, 0.027] GeV, D0->K+pi-pi+pi- [-0.026, 0.024] GeV, D0->K+pi-pi0 [-0.062, 0.049] GeV")
     .note(:ks0_reconstruction, "K_S0 candidates are formed from pi+pi- with a secondary-vertex fit and accepted within a mass window around the nominal M(K_S0)")
     .note(:background_veto, "extra-photon E_gamma_max veto and an M(K_S0 pi- e+) window applied to suppress the D0 -> K_S0 pi- pi+ pi0 peaking background")
     .with_decay_card(decay_card_e)

# Tag side: anti-D0 (charm +1) reconstructed from three hadronic tag modes;
# explicit M_BC window (opt-in override of the default store-not-cut)
alg_e.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPiPi, :D0toKPiPi0
  t.charm 1
  t.window :mBC, min: 1.859, max: 1.873
end

# Signal side: one pi+, two pi-, one e+ (net charge 0) + missing massless nu_e
alg_e.signal_side do |s|
  s.charged(pip: 1, pim: 2, ep: 1)
  s.require_charge 0
  s.missing :nu_e
end

# Four-momentum constraint with chi2 < 200
alg_e.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_e.apply
alg_e.execute_on([d3770_data, d3770_incMC, exMC_e])

# -------------------- Muon mode --------------------
alg_mu = TagAnalysis.new("D0KsPimMuNu")
alg_mu.set_header(["D0KsPimMuNuAlg/D0KsPimMuNu.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      .note(:deltaE_windows, "mode-dependent tag-side deltaE windows applied in the tag scan (not expressible as a single per-side window declaration): D0->K+pi- [-0.027, 0.027] GeV, D0->K+pi-pi+pi- [-0.026, 0.024] GeV, D0->K+pi-pi0 [-0.062, 0.049] GeV")
      .note(:ks0_reconstruction, "K_S0 candidates are formed from pi+pi- with a secondary-vertex fit and accepted within a mass window around the nominal M(K_S0)")
      .note(:background_veto, "extra-photon E_gamma_max veto and an M(K_S0 pi- mu+ (pi0)) window applied to suppress the D0 -> K_S0 pi- pi+ pi0 peaking background")
      .with_decay_card(decay_card_mu)

alg_mu.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPiPi, :D0toKPiPi0
  t.charm 1
  t.window :mBC, min: 1.859, max: 1.873
end

# Signal side: one pi+, two pi-, one mu+ (net charge 0) + missing massless nu_mu
alg_mu.signal_side do |s|
  s.charged(pip: 1, pim: 2, mup: 1)
  s.require_charge 0
  s.missing :nu_mu
end

alg_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mu.apply
alg_mu.execute_on([d3770_data, d3770_incMC, exMC_mu])