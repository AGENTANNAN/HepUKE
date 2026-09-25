# =============================================================================
# ψ(3770) → D⁺D⁻  —  double-tag (DT) search for
#   D⁺ → γρ⁺  with ρ⁺ → π⁺π⁰   (π⁰ → γγ)   [signal mode I]
#   D⁺ → γK*⁺ with K*⁺ → K⁺π⁰  (π⁰ → γγ)   [signal mode II]
# Tag side: six hadronic D⁻ modes, charm −1.
# 20.3 fb⁻¹ real data @ 3.773 GeV + inclusive MC + 1M exclusive MC per signal mode.
# =============================================================================

### ---------------------------- Datasets ---------------------------- ###
data_3773  = DatasetManager.real_data.find("712_3773")     # 20.3 fb⁻¹ ψ(3770) data @ 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # matching inclusive MC sample

### ------------------------ Decay cards (signal MC) ------------------------ ###
# Tag D⁻ fixed to K⁺π⁻π⁻ in the exclusive MC; the signal D⁺ decay differs per mode.

# Mode I: D⁺ → γρ⁺, ρ⁺ → π⁺π⁰, π⁰ → γγ
decay_card_rho = <<~DECAYCARD
    Decay psi(3770)
    1.000  D+  D-  PHSP;
    Enddecay

    Decay D-
    1.000  K+  pi-  pi-  PHSP;
    Enddecay

    Decay D+
    1.000  gamma  rho+  PHSP;
    Enddecay

    Decay rho+
    1.000  pi+  pi0  VSS;
    Enddecay

    Decay pi0
    1.000  gamma  gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: D⁺ → γK*⁺, K*⁺ → K⁺π⁰, π⁰ → γγ
decay_card_Kstar = <<~DECAYCARD
    Decay psi(3770)
    1.000  D+  D-  PHSP;
    Enddecay

    Decay D-
    1.000  K+  pi-  pi-  PHSP;
    Enddecay

    Decay D+
    1.000  gamma  K*+  PHSP;
    Enddecay

    Decay K*+
    1.000  K+  pi0  VSS;
    Enddecay

    Decay pi0
    1.000  gamma  gamma  PHSP;
    Enddecay

    End
DECAYCARD

### -------------------- Exclusive MC (1M events / mode) -------------------- ###
exMC_rho = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_Dp_gamma_rho"
  config.related_dataset = data_3773           # associate with the real data sample
  config.events          = 1_000_000           # 1 million events
  config.decay_card      = decay_card_rho
  config.cross_section   = :default
end

exMC_Kstar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_Dp_gamma_Kstar"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_Kstar
  config.cross_section   = :default
end

### ====================================================================== ###
###  Signal mode I : D⁺ → γρ⁺ ,  ρ⁺ → π⁺π⁰                               ###
### ====================================================================== ###
alg_rho = TagAnalysis.new("DpToGammaRho")
alg_rho.set_header(["DpToGammaRhoAlg/DpToGammaRho.h"])
       .set_constant({"ECMS" => [:double, 3.773]})
       .with_decay_card(decay_card_rho)
       # Standard tag-side (DTagAlg) quality configuration — not expressible as DSL cuts
       .note(:tag_track_quality,
             "tag side uses the standard charged-track quality cuts |cosθ| < 0.93, |Vz| < 10 cm, Vr < 1 cm (DTagAlg defaults)")
       .note(:tag_photon_thresholds,
             "tag reconstruction uses the standard EMC photon thresholds 25 MeV (barrel) / 50 MeV (endcap) (DTagAlg defaults)")
       .note(:tag_pid_method,
             "tag side uses probability-based PID (DTagAlg standard recipe)")

# Tag side: six hadronic D⁻ modes, pin the D⁻ (charm −1)
alg_rho.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKKPi, :DptoKsPi, :DptoKsPiPi, :DptoKsPiPiPi
  t.charm -1
end

# Signal side: 3 photons (1 radiative + 2 from π⁰), isolated by ≥10°, plus exactly one π⁺
alg_rho.signal_side do |s|
  s.photons 3
  s.min_photon_angle 10.0
  s.charged(pip: 1)
end

# Kinematic fit: 4C constraint + π⁰ mass constraint on the two photons, χ² < 200
alg_rho.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_rho.apply                       # tag spec validation/render (takes no Selection)
root_files_rho = alg_rho.execute_on([data_3773, incMC_3773, exMC_rho])

### ====================================================================== ###
###  Signal mode II : D⁺ → γK*⁺ ,  K*⁺ → K⁺π⁰                              ###
### ====================================================================== ###
alg_Kstar = TagAnalysis.new("DpToGammaKstar")
alg_Kstar.set_header(["DpToGammaKstarAlg/DpToGammaKstar.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_Kstar)
         .note(:tag_track_quality,
               "tag side uses the standard charged-track quality cuts |cosθ| < 0.93, |Vz| < 10 cm, Vr < 1 cm (DTagAlg defaults)")
         .note(:tag_photon_thresholds,
               "tag reconstruction uses the standard EMC photon thresholds 25 MeV (barrel) / 50 MeV (endcap) (DTagAlg defaults)")
         .note(:tag_pid_method,
               "tag side uses probability-based PID (DTagAlg standard recipe)")

# Same tag side as mode I
alg_Kstar.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKKPi, :DptoKsPi, :DptoKsPiPi, :DptoKsPiPiPi
  t.charm -1
end

# Signal side: 3 photons (1 radiative + 2 from π⁰), isolated by ≥10°, plus exactly one K⁺
alg_Kstar.signal_side do |s|
  s.photons 3
  s.min_photon_angle 10.0
  s.charged(kp: 1)
end

# Same fit procedure as mode I
alg_Kstar.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_Kstar.apply
root_files_Kstar = alg_Kstar.execute_on([data_3773, incMC_3773, exMC_Kstar])