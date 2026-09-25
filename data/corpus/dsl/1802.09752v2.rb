# 1802.09752v2: D → h(h') e⁺ e⁻ rare decays at ψ(3770)
# Tag-based double-tag analysis

psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ---- Common decay card for exclusive MC (D+ signal, D0 signal) ----
# D+ → π⁺π⁰ e⁺e⁻
decay_card_dp_pi_pi0_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 pi+ pi0 e+ e- PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_dp_pi_pi0_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_to_pi_pi0_ee_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_dp_pi_pi0_ee
  config.cross_section = :default
end

# D+ → K⁺π⁰ e⁺e⁻
decay_card_dp_k_pi0_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K+ pi0 e+ e- PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_dp_k_pi0_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_to_K_pi0_ee_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_dp_k_pi0_ee
  config.cross_section = :default
end

# D+ → K_S⁰ π⁺ e⁺e⁻
decay_card_dp_ks_pi_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K_S0 pi+ e+ e- PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_dp_ks_pi_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_to_Ks_pi_ee_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_dp_ks_pi_ee
  config.cross_section = :default
end

# D+ → K_S⁰ K⁺ e⁺e⁻
decay_card_dp_ks_k_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K_S0 K+ e+ e- PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_dp_ks_k_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_to_Ks_K_ee_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_dp_ks_k_ee
  config.cross_section = :default
end

# D0 → K⁻K⁺ e⁺e⁻
decay_card_d0_kk_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K- K+ e+ e- PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_d0_kk_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0_to_KK_ee_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_d0_kk_ee
  config.cross_section = :default
end

# D0 → π⁺π⁻ e⁺e⁻
decay_card_d0_pipi_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi+ pi- e+ e- PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_d0_pipi_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0_to_pipi_ee_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_d0_pipi_ee
  config.cross_section = :default
end

# D0 → K⁻π⁺ e⁺e⁻
decay_card_d0_kpi_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ e+ e- PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_d0_kpi_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0_to_Kpi_ee_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_d0_kpi_ee
  config.cross_section = :default
end

# D0 → π⁰ e⁺e⁻
decay_card_d0_pi0_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi0 e+ e- PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_d0_pi0_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0_to_pi0_ee_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_d0_pi0_ee
  config.cross_section = :default
end

# D0 → η e⁺e⁻
decay_card_d0_eta_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 eta e+ e- PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_d0_eta_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0_to_eta_ee_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_d0_eta_ee
  config.cross_section = :default
end

# D0 → ω e⁺e⁻
decay_card_d0_omega_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 omega e+ e- PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_d0_omega_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0_to_omega_ee_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_d0_omega_ee
  config.cross_section = :default
end

# D0 → K_S⁰ e⁺e⁻
decay_card_d0_ks_ee = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K_S0 e+ e- PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_d0_ks_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0_to_Ks_ee_exclusive_mc"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_d0_ks_ee
  config.cross_section = :default
end

# ============================================================
# D+ signal modes (tag D-)
# ============================================================

# D+ → π⁺π⁰ e⁺e⁻
alg_dp_pi_pi0_ee = TagAnalysis.new("DpToPiPi0EE")
alg_dp_pi_pi0_ee.set_header(["DpPiPi0EEAlg/DpPiPi0EE.h"])
                 .set_constant({"ECMS" => [:double, 3.773]})
                 .note(:background_veto, "φ veto: M(e⁺e⁻) outside (0.935, 1.053) GeV/c²; e⁺e⁻ vertex R_xy veto outside (2.0, 8.0) cm to suppress γ-conversion")
                 .note(:tag_mode_unavailable, "E/pc > 0.8 and combined-e-PID (L_e/(L_e+L_pi+L_K) > 0.8) criteria for e± not expressible in tag DSL; electron PID uses fixed v1 thresholds")
                 .with_decay_card(decay_card_dp_pi_pi0_ee)

alg_dp_pi_pi0_ee.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :deltaE, min: -0.060, max: 0.034   # approximate 3σ for modes with and without π⁰
end

alg_dp_pi_pi0_ee.signal_side do |s|
  s.charged(pip: 1, ep: 1, em: 1)
  s.photons 2
end

alg_dp_pi_pi0_ee.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_dp_pi_pi0_ee.apply
alg_dp_pi_pi0_ee.execute_on([psi3770_data, psi3770_incMC, exMC_dp_pi_pi0_ee])

# D+ → K⁺π⁰ e⁺e⁻
alg_dp_k_pi0_ee = TagAnalysis.new("DpToKPi0EE")
alg_dp_k_pi0_ee.set_header(["DpKPi0EEAlg/DpKPi0EE.h"])
                .set_constant({"ECMS" => [:double, 3.773]})
                .note(:background_veto, "φ veto: M(e⁺e⁻) outside (0.935, 1.053) GeV/c²; e⁺e⁻ vertex R_xy veto outside (2.0, 8.0) cm to suppress γ-conversion")
                .note(:tag_mode_unavailable, "E/pc > 0.8 and combined-e-PID criteria for e± not expressible in tag DSL")
                .with_decay_card(decay_card_dp_k_pi0_ee)

alg_dp_k_pi0_ee.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :deltaE, min: -0.060, max: 0.034
end

alg_dp_k_pi0_ee.signal_side do |s|
  s.charged(kp: 1, ep: 1, em: 1)
  s.photons 2
end

alg_dp_k_pi0_ee.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_dp_k_pi0_ee.apply
alg_dp_k_pi0_ee.execute_on([psi3770_data, psi3770_incMC, exMC_dp_k_pi0_ee])

# D+ → K_S⁰ π⁺ e⁺e⁻
alg_dp_ks_pi_ee = TagAnalysis.new("DpToKsPiEE")
alg_dp_ks_pi_ee.set_header(["DpKsPiEEAlg/DpKsPiEE.h"])
                .set_constant({"ECMS" => [:double, 3.773]})
                .note(:background_veto, "φ veto: M(e⁺e⁻) outside (0.935, 1.053) GeV/c²; e⁺e⁻ vertex R_xy veto outside (2.0, 8.0) cm to suppress γ-conversion")
                .note(:tag_mode_unavailable, "E/pc > 0.8 and combined-e-PID criteria for e± not expressible in tag DSL; K_S⁰ uses L/σ_L > 2 requirement")
                .with_decay_card(decay_card_dp_ks_pi_ee)

alg_dp_ks_pi_ee.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :deltaE, min: -0.060, max: 0.034
end

alg_dp_ks_pi_ee.signal_side do |s|
  s.charged(pip: 1, ep: 1, em: 1)
end

alg_dp_ks_pi_ee.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_dp_ks_pi_ee.apply
alg_dp_ks_pi_ee.execute_on([psi3770_data, psi3770_incMC, exMC_dp_ks_pi_ee])

# D+ → K_S⁰ K⁺ e⁺e⁻
alg_dp_ks_k_ee = TagAnalysis.new("DpToKsKEE")
alg_dp_ks_k_ee.set_header(["DpKsKEEAlg/DpKsKEE.h"])
               .set_constant({"ECMS" => [:double, 3.773]})
               .note(:background_veto, "φ veto: M(e⁺e⁻) outside (0.935, 1.053) GeV/c²; e⁺e⁻ vertex R_xy veto outside (2.0, 8.0) cm to suppress γ-conversion")
               .note(:tag_mode_unavailable, "E/pc > 0.8 and combined-e-PID criteria for e± not expressible in tag DSL; K_S⁰ uses L/σ_L > 2 requirement")
               .with_decay_card(decay_card_dp_ks_k_ee)

alg_dp_ks_k_ee.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :deltaE, min: -0.060, max: 0.034
end

alg_dp_ks_k_ee.signal_side do |s|
  s.charged(kp: 1, ep: 1, em: 1)
end

alg_dp_ks_k_ee.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_dp_ks_k_ee.apply
alg_dp_ks_k_ee.execute_on([psi3770_data, psi3770_incMC, exMC_dp_ks_k_ee])

# ============================================================
# D0 signal modes (tag anti-D0)
# ============================================================

# D0 → K⁻K⁺ e⁺e⁻
alg_d0_kk_ee = TagAnalysis.new("D0ToKKEE")
alg_d0_kk_ee.set_header(["D0KKEEAlg/D0KKEE.h"])
             .set_constant({"ECMS" => [:double, 3.773]})
             .note(:background_veto, "φ veto: M(e⁺e⁻) outside (0.935, 1.053) GeV/c²; e⁺e⁻ vertex R_xy veto outside (2.0, 8.0) cm to suppress γ-conversion")
             .note(:tag_mode_unavailable, "E/pc > 0.8 and combined-e-PID criteria for e± not expressible in tag DSL")
             .with_decay_card(decay_card_d0_kk_ee)

alg_d0_kk_ee.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :deltaE, min: -0.064, max: 0.035   # looser covering modes with π⁰
end

alg_d0_kk_ee.signal_side do |s|
  s.charged(km: 1, kp: 1, ep: 1, em: 1)
end

alg_d0_kk_ee.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_kk_ee.apply
alg_d0_kk_ee.execute_on([psi3770_data, psi3770_incMC, exMC_d0_kk_ee])

# D0 → π⁺π⁻ e⁺e⁻
alg_d0_pipi_ee = TagAnalysis.new("D0ToPiPiEE")
alg_d0_pipi_ee.set_header(["D0PiPiEEAlg/D0PiPiEE.h"])
               .set_constant({"ECMS" => [:double, 3.773]})
               .note(:background_veto, "φ veto: M(e⁺e⁻) outside (0.935, 1.053) GeV/c²; e⁺e⁻ vertex R_xy veto outside (2.0, 8.0) cm to suppress γ-conversion")
               .note(:tag_mode_unavailable, "E/pc > 0.8 and combined-e-PID criteria for e± not expressible in tag DSL")
               .with_decay_card(decay_card_d0_pipi_ee)

alg_d0_pipi_ee.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :deltaE, min: -0.064, max: 0.035
end

alg_d0_pipi_ee.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1, em: 1)
end

alg_d0_pipi_ee.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_pipi_ee.apply
alg_d0_pipi_ee.execute_on([psi3770_data, psi3770_incMC, exMC_d0_pipi_ee])

# D0 → K⁻π⁺ e⁺e⁻
alg_d0_kpi_ee = TagAnalysis.new("D0ToKPiEE")
alg_d0_kpi_ee.set_header(["D0KPiEEAlg/D0KPiEE.h"])
              .set_constant({"ECMS" => [:double, 3.773]})
              .note(:background_veto, "φ veto: M(e⁺e⁻) outside (0.935, 1.053) GeV/c²; e⁺e⁻ vertex R_xy veto outside (2.0, 8.0) cm to suppress γ-conversion")
              .note(:tag_mode_unavailable, "E/pc > 0.8 and combined-e-PID criteria for e± not expressible in tag DSL")
              .with_decay_card(decay_card_d0_kpi_ee)

alg_d0_kpi_ee.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :deltaE, min: -0.064, max: 0.035
end

alg_d0_kpi_ee.signal_side do |s|
  s.charged(km: 1, pip: 1, ep: 1, em: 1)
end

alg_d0_kpi_ee.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_kpi_ee.apply
alg_d0_kpi_ee.execute_on([psi3770_data, psi3770_incMC, exMC_d0_kpi_ee])

# D0 → π⁰ e⁺e⁻
alg_d0_pi0_ee = TagAnalysis.new("D0ToPi0EE")
alg_d0_pi0_ee.set_header(["D0Pi0EEAlg/D0Pi0EE.h"])
              .set_constant({"ECMS" => [:double, 3.773]})
              .note(:background_veto, "φ veto: M(e⁺e⁻) outside (0.935, 1.053) GeV/c²; e⁺e⁻ vertex R_xy veto outside (2.0, 8.0) cm to suppress γ-conversion")
              .note(:tag_mode_unavailable, "E/pc > 0.8 and combined-e-PID criteria for e± not expressible in tag DSL")
              .with_decay_card(decay_card_d0_pi0_ee)

alg_d0_pi0_ee.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :deltaE, min: -0.064, max: 0.035
end

alg_d0_pi0_ee.signal_side do |s|
  s.photons 2
  s.charged(ep: 1, em: 1)
end

alg_d0_pi0_ee.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_d0_pi0_ee.apply
alg_d0_pi0_ee.execute_on([psi3770_data, psi3770_incMC, exMC_d0_pi0_ee])

# D0 → η e⁺e⁻
alg_d0_eta_ee = TagAnalysis.new("D0ToEtaEE")
alg_d0_eta_ee.set_header(["D0EtaEEAlg/D0EtaEE.h"])
              .set_constant({"ECMS" => [:double, 3.773]})
              .note(:background_veto, "φ veto: M(e⁺e⁻) outside (0.935, 1.053) GeV/c²; e⁺e⁻ vertex R_xy veto outside (2.0, 8.0) cm to suppress γ-conversion")
              .note(:tag_mode_unavailable, "E/pc > 0.8 and combined-e-PID criteria for e± not expressible in tag DSL")
              .with_decay_card(decay_card_d0_eta_ee)

alg_d0_eta_ee.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :deltaE, min: -0.064, max: 0.035
end

alg_d0_eta_ee.signal_side do |s|
  s.photons 2
  s.charged(ep: 1, em: 1)
end

alg_d0_eta_ee.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_d0_eta_ee.apply
alg_d0_eta_ee.execute_on([psi3770_data, psi3770_incMC, exMC_d0_eta_ee])

# D0 → ω e⁺e⁻
alg_d0_omega_ee = TagAnalysis.new("D0ToOmegaEE")
alg_d0_omega_ee.set_header(["D0OmegaEEAlg/D0OmegaEE.h"])
                .set_constant({"ECMS" => [:double, 3.773]})
                .note(:background_veto, "φ veto: M(e⁺e⁻) outside (0.935, 1.053) GeV/c²; e⁺e⁻ vertex R_xy veto outside (2.0, 8.0) cm to suppress γ-conversion")
                .note(:tag_mode_unavailable, "E/pc > 0.8 and combined-e-PID criteria for e± not expressible in tag DSL; ω → π⁺π⁻π⁰ with M_3pi window (0.720, 0.840) GeV/c² applied in ROOT")
                .with_decay_card(decay_card_d0_omega_ee)

alg_d0_omega_ee.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :deltaE, min: -0.064, max: 0.035
end

alg_d0_omega_ee.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1, em: 1)
  s.photons 2
end

alg_d0_omega_ee.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_d0_omega_ee.apply
alg_d0_omega_ee.execute_on([psi3770_data, psi3770_incMC, exMC_d0_omega_ee])

# D0 → K_S⁰ e⁺e⁻
alg_d0_ks_ee = TagAnalysis.new("D0ToKsEE")
alg_d0_ks_ee.set_header(["D0KsEEAlg/D0KsEE.h"])
             .set_constant({"ECMS" => [:double, 3.773]})
             .note(:background_veto, "φ veto: M(e⁺e⁻) outside (0.935, 1.053) GeV/c²; e⁺e⁻ vertex R_xy veto outside (2.0, 8.0) cm to suppress γ-conversion")
             .note(:tag_mode_unavailable, "E/pc > 0.8 and combined-e-PID criteria for e± not expressible in tag DSL; K_S⁰ uses L/σ_L > 2 requirement")
             .with_decay_card(decay_card_d0_ks_ee)

alg_d0_ks_ee.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :deltaE, min: -0.064, max: 0.035
end

alg_d0_ks_ee.signal_side do |s|
  s.charged(ep: 1, em: 1)
end

alg_d0_ks_ee.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_ks_ee.apply
alg_d0_ks_ee.execute_on([psi3770_data, psi3770_incMC, exMC_d0_ks_ee])