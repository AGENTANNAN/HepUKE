# 2411.07730v2: D⁰ → a₀(980)⁻ e⁺ ν_e at ψ(3770)
# Tag-based ST+missing semileptonic analysis
# Similar structure to 1803.02166v1.rb D⁰ channel

psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ---- Decay card for exclusive MC (D0 signal) ----
# D0 → a₀(980)⁻ e⁺ ν_e, a₀(980)⁻ → ηπ⁻, η → γγ
decay_card_d0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 a_0(980)- e+ nu_e PHSP;
    Enddecay

    Decay a_0(980)-
    1.000 eta pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_d0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0_to_a0_etapi_e_nu_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_d0
  config.cross_section = :default
end

# ============================================================
# D0 → a₀(980)⁻ e⁺ ν_e (tag anti-D0)
# ============================================================
alg_d0_a0 = TagAnalysis.new("D0ToA0980minuseNu")
alg_d0_a0.set_header(["D0A0ENuAlg/D0A0ENu.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .note(:tag_mode_unavailable, "E/pc > 0.8 and combined-e-PID (L_e'/(L_e'+L_pi'+L_K') > 0.8) criteria for e± not expressible in tag DSL; electron PID uses fixed v1 thresholds")
          .note(:background_veto, "K_L⁰ veto via lateral moment cut on higher-energy photon from η; extra π⁰ veto; extra charged track veto; all applied in ROOT")
          .note(:signal_extraction, "2-D unbinned ML fit to M(ηπ) vs U in ROOT; Flatté formula for a₀(980) line shape; U = E_miss - c|p_miss|")
          .note(:paper_info, "7.93 fb⁻¹ at ψ(3770); improved statistics over 1803.02166v1")
          .with_decay_card(decay_card_d0)

alg_d0_a0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :deltaE, min: -0.064, max: 0.035
end

alg_d0_a0.signal_side do |s|
  s.charged(pim: 1, ep: 1)
  s.photons 2
  s.require_charge 0
  s.missing :nu_e
end

alg_d0_a0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_d0_a0.apply
alg_d0_a0.execute_on([psi3770_data, psi3770_incMC, exMC_d0])