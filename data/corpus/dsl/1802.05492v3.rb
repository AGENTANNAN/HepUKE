# 1802.05492v3: D^{0(+)} → π^{-(0)} μ^{+} ν_{μ} at ψ(3770)
# Tag-based semileptonic analysis via ST D mesons

psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ---- Decay card for exclusive MC (D0 signal) ----
decay_card_d0_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi- mu+ nu_mu PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_d0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0_to_pi_mu_nu_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_d0_signal
  config.cross_section = :default
end

# ---- Decay card for exclusive MC (D+ signal) ----
decay_card_dp_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 pi0 mu+ nu_mu PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_dp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_to_pi0_mu_nu_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_dp_signal
  config.cross_section = :default
end

# ============================================================
# Algorithm 1: D0 → π⁻ μ⁺ ν_μ (tag anti-D0 with hadronic modes)
# ============================================================
alg_d0_mu = TagAnalysis.new("D0ToPiMuNu")
alg_d0_mu.set_header(["D0PiMuNuAlg/D0PiMuNu.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .note(:tag_mode_unavailable, "MUC hit-depth/polar-angle dependent muon ID criteria not expressible in tag DSL; muon PID uses fixed v1 thresholds")
          .note(:background_veto, "M(π⁻μ⁺) < 1.7 GeV/c² and E_max_extra_gamma < 0.07 GeV applied on signal side to suppress hadronic backgrounds; K_S⁰ veto via M(π⁻μ⁺) outside (0.46,0.50) GeV/c²")
          .with_decay_card(decay_card_d0_signal)

# Tag side: anti-D0 reconstructed via three hadronic modes
alg_d0_mu.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :deltaE, min: -0.055, max: 0.040   # for modes with π⁰
  # Note: for Kπ mode without π⁰, deltaE window is (-0.025, 0.025)
  # but DSL allows only one window per observable per side; use the looser
  # ST mBC signal region: (1.859, 1.873) GeV/c² applied in ROOT
end

# Signal side: D0 → π⁻ μ⁺ ν_μ
alg_d0_mu.signal_side do |s|
  s.charged(pim: 1, mup: 1)
  s.require_charge 0    # -1 + 1 = 0
  s.missing :nu_mu      # massless neutrino
end

alg_d0_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_mu.apply
alg_d0_mu.execute_on([psi3770_data, psi3770_incMC, exMC_d0])

# ============================================================
# Algorithm 2: D+ → π⁰ μ⁺ ν_μ (tag D- with hadronic modes)
# ============================================================
alg_dp_mu = TagAnalysis.new("DpToPi0MuNu")
alg_dp_mu.set_header(["DpPi0MuNuAlg/DpPi0MuNu.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .note(:tag_mode_unavailable, "MUC hit-depth/polar-angle dependent muon ID criteria not expressible in tag DSL; muon PID uses fixed v1 thresholds")
          .note(:background_veto, "M(π⁰μ⁺) < 1.7 GeV/c² and E_max_extra_gamma < 0.07 GeV applied; K_S⁰ veto via M_recoil(D⁻μ⁺) outside (0.45,0.55) GeV/c²")
          .with_decay_card(decay_card_dp_signal)

# Tag side: D- reconstructed via six hadronic modes
alg_dp_mu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :deltaE, min: -0.055, max: 0.040   # looser window covering modes with π⁰
  # ST mBC signal region: (1.863, 1.877) GeV/c² applied in ROOT
end

# Signal side: D+ → π⁰ μ⁺ ν_μ
# π⁰ → γγ reconstructed from photon pair
alg_dp_mu.signal_side do |s|
  s.photons 2
  s.charged(mup: 1)
  s.require_charge 1    # μ⁺ = +1
  s.missing :nu_mu      # massless neutrino
end

alg_dp_mu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_dp_mu.apply
alg_dp_mu.execute_on([psi3770_data, psi3770_incMC, exMC_dp])