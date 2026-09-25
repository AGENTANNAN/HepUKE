# =============================================================================
#  ψ(3770) semileptonic D decays with hadronic tag-side reconstruction
#    channel 1 :  D0  → π⁻ μ⁺ ν_μ     (tag the anti-D0)
#    channel 2 :  D+  → π⁰ μ⁺ ν_μ     (tag the D⁻)
# =============================================================================

### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")      # real ψ(3770) data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # matching inclusive MC

# decay card — ψ(3770) → D0 anti-D0, D0 → π⁻ μ⁺ ν_μ (signal), anti-D0 → K⁺ π⁻ (tag)
decay_card_d0 = <<~DECAYCARD
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

# decay card — ψ(3770) → D+ D-, D+ → π⁰ μ⁺ ν_μ (signal), D- → K⁺ π⁻ π⁻ (tag)
decay_card_dp = <<~DECAYCARD
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

# 200k-event exclusive MC samples, one per channel
exMC_d0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0_pimupinu"
  config.related_dataset = data_3773
  config.events          = 200_000
  config.decay_card      = decay_card_d0
  config.cross_section   = :default
end

exMC_dp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dp_pi0mupinu"
  config.related_dataset = data_3773
  config.events          = 200_000
  config.decay_card      = decay_card_dp
  config.cross_section   = :default
end

### Channel 1 — D0 → π⁻ μ⁺ ν_μ (tag the anti-D0) ###
alg_name_d0 = "D0SemilepTag"
alg_d0 = TagAnalysis.new(alg_name_d0)
alg_d0.set_header(["#{alg_name_d0}Alg/#{alg_name_d0}.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_d0)

# tag side — hadronic anti-D0 taken from the pre-stored DTag candidates
alg_d0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi      # Kπ, Kππ⁰, Kπππ
  t.charm -1                                       # tag the anti-D0
  t.window :deltaE, min: -0.055, max: 0.040        # ΔE ∈ (−0.055, 0.040) GeV
  t.window :mBC,    min: 1.859,  max: 1.873        # mBC ∈ (1.859, 1.873) GeV/c²
end

# signal side — exactly one π⁻ and one μ⁺, net charge 0, one missing neutrino
alg_d0.signal_side do |s|
  s.charged(pim: 1, mup: 1)
  s.require_charge 0
  s.missing :nu_mu
end

# kinematic fit — 4C constraint, χ² < 200; M(π⁻μ⁺) window < 1.7 GeV/c²
alg_d0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pim, :mup).between(0.0, 1.7)
  f.chi2_cut 200
end

alg_d0
  .note(:background_veto,
        "K_S0 veto: events whose M(pi- mu+) falls inside (0.46, 0.50) GeV/c^2 " \
        "are rejected to suppress K_S0 -> pi+ pi- contamination")
  .note(:background_veto,
        "extra-photon veto: the energy of any photon not associated with the " \
        "signal must be < 0.07 GeV (suppresses backgrounds with extra pi0 / gamma)")
  .note(:pid_correction_method,
        "muon identification uses the fixed v1 ParticleID thresholds; the detailed " \
        "MUC hit-depth / polar-angle criteria are not expressible in the DSL")
alg_d0.apply

### Channel 2 — D+ → π⁰ μ⁺ ν_μ (tag the D⁻) ###
alg_name_dp = "DpSemilepTag"
alg_dp = TagAnalysis.new(alg_name_dp)
alg_dp.set_header(["#{alg_name_dp}Alg/#{alg_name_dp}.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_dp)

# tag side — hadronic D⁻ taken from the pre-stored DTag candidates
alg_dp.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi    # Kππ, K_Sπ, Kπππ⁰, K_Sππ⁰, K_Sπππ, KKπ
  t.window :deltaE, min: -0.055, max: 0.040         # ΔE ∈ (−0.055, 0.040) GeV
  t.window :mBC,    min: 1.863,  max: 1.877         # mBC ∈ (1.863, 1.877) GeV/c²
end

# signal side — π⁰ → γγ (two photons), one μ⁺, net charge +1, one missing neutrino
alg_dp.signal_side do |s|
  s.photons 2
  s.charged(mup: 1)
  s.require_charge 1
  s.missing :nu_mu
end

# kinematic fit — 4C constraint, π⁰ mass constraint on the γγ pair, χ² < 200
alg_dp.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_dp
  .note(:background_veto,
        "M(pi0 mu+) < 1.7 GeV/c^2 requirement, suppressing hadronic D backgrounds")
  .note(:background_veto,
        "extra-photon veto: the energy of any photon beyond the pi0 -> gamma gamma " \
        "pair must be < 0.07 GeV")
  .note(:background_veto,
        "K_S0 veto: events whose recoil mass M(recoil D- mu+) falls inside " \
        "(0.45, 0.55) GeV/c^2 are rejected")
  .note(:pid_correction_method,
        "muon identification uses the fixed v1 ParticleID thresholds; the detailed " \
        "MUC hit-depth / polar-angle criteria are not expressible in the DSL")
alg_dp.apply

### Execute ###
root_files_d0 = alg_d0.execute_on([data_3773, incMC_3773, exMC_d0])
root_files_dp = alg_dp.execute_on([data_3773, incMC_3773, exMC_dp])