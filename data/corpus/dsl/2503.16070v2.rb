# DSL: 2503.16070v2 — Search for D+ → γ e+ ν_e at ψ(3770), 20.3 fb-1
# Double-tag (DT) method: tag D- via 6 hadronic modes, signal side γ e+ + missing ν_e
# TagAnalysis — semileptonic ST+missing pattern

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ── Exclusive MC: D+ → γ e+ ν_e ──────────────────────────────────────
decay_card_sig = <<~DECAYCARD
  Decay D+
  1.000 gamma e+ nu_e PHSP;
  Enddecay
DECAYCARD

exMC_sig = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Dp_gamma_e_nue"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_sig
  config.cross_section   = :default
end

# ── Exclusive MC: validation channel D+ → π0 e+ ν_e ──────────────────
decay_card_val = <<~DECAYCARD
  Decay D+
  1.000 pi0 e+ nu_e PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
DECAYCARD

exMC_val = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Dp_pi0_e_nue"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_val
  config.cross_section   = :default
end

# ── TagAnalysis: D+ → γ e+ ν_e (signal) ──────────────────────────────
alg_sig = TagAnalysis.new("DpToGammaENuE")
alg_sig.set_header(["DpToGammaENuEAlg/DpToGammaENuE.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })

# Tag side: D- via 6 hadronic modes
alg_sig.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  # Paper Table 2: mode-dependent ΔE windows ~±25-62 MeV (≈±3.5σ)
  # Use a loose common window; precise windows applied in ROOT
  t.window :deltaE, abs: 0.06
end

# Signal side: one positron, one radiative photon, missing neutrino
alg_sig.signal_side do |s|
  s.photons 1
  s.charged(ep: 1)
  s.missing :nu_e              # massless — semileptonic ST+missing pattern
  s.require_charge 1           # e+ has charge +1
end

alg_sig.fit do |f|
  f.constrain_four_momentum    # 4C: tag + γ + e+ + ν_e = ecms_lab
  f.chi2_cut 200
end

# Inexpressible BOSS-side procedures
alg_sig.note(:dnn_signal_id, "DNN (ParT architecture) for signal/background discrimination; 3-class output: signal/pi0enu/other; DNN score cuts applied in ROOT")
alg_sig.note(:fsr_recovery, "FSR recovery: photons within 5° of positron track added back to e+ momentum")
alg_sig.note(:umiss_window, "U_miss ∈ [-0.2, 0.2] GeV window applied in ROOT analysis")
alg_sig.note(:photon_energy, "Most energetic photon selected if multiple candidates pass selection")
alg_sig.note(:signal_yield_extraction, "DT yield from U_miss fit; upper limit set via Bayesian method")

alg_sig.apply
alg_sig.execute_on([data_3773, incMC_3773, exMC_sig])

# ── TagAnalysis: D+ → π0 e+ ν_e (validation channel) ─────────────────
# Note: this validation channel requires π0 reconstruction which is not
# fully expressible in the TagAnalysis signal side; captured as note.
alg_val = TagAnalysis.new("DpToPi0ENuE")
alg_val.set_header(["DpToPi0ENuEAlg/DpToPi0ENuE.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })

alg_val.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :deltaE, abs: 0.06
end

alg_val.signal_side do |s|
  s.photons 2                 # π0 → γγ
  s.charged(ep: 1)
  s.missing :nu_e
  s.require_charge 1
end

alg_val.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_val.note(:pi0_reconstruction, "π0 → γγ: best candidate chosen by minimum χ² from kinematic fit; M(γγ) ∈ [115, 150] MeV/c² before fit")
alg_val.note(:dnn_validation, "DNN vetoes partially reversed for validation: require DNN pi0enu score > 0.15, other-bkg score < 0.05")
alg_val.note(:umiss_prime, "U'_miss observable for validation channel")

alg_val.apply
alg_val.execute_on([data_3773, incMC_3773, exMC_val])