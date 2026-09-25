# frozen_string_literal: true
### Dataset description ###
data_3773 = DatasetManager.real_data.find("712_3773")       # ψ(3770) real data (7.93 fb⁻¹ at √s = 3.773 GeV)
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")    # Matching inclusive MC

# Decay card for the signal process: D+ → K_S0 π0 e+ ν_e (K_S0 → π+π-, π0 → γγ)
decay_card_signal = <<~DECAYCARD
  Decay D+
  1.000 K_S0 pi0 e+ nu_e PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for the signal channel (500k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_Dp_KsPi0Enu"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
alg_name = "DpKsPi0EnuDTag"
tag_alg = TagAnalysis.new(alg_name)
tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })               # √s = 3.773 GeV
       .set_alias({ "std::vector<double>" => "Vdouble" })
       .with_decay_card(decay_card_signal)

# Tag side: D- reconstructed in six hadronic modes (charm pinned to -1)
tag_alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  # Explicit tag-side mBC window requested by the analysis
  t.window :mBC, min: 1.863, max: 1.877
end

# Signal side: D+ → K_S0 π0 e+ ν_e, with K_S0 → π+π- and π0 → γγ
# (2-48 photons, one π+, one π-, one e+, net charge +1, one missing massless ν_e)
tag_alg.signal_side do |s|
  s.photons 2..48
  s.charged(pip: 1, pim: 1, ep: 1)
  s.missing :nu_e
  s.require_charge 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# Kinematic fit: 4-momentum conservation + K_S0 and π0 mass constraints
tag_alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # π0 → γγ
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)      # K_S0 → π+π-
  f.chi2_cut 200
end

# The tag-side ΔE requirement is mode dependent: it cannot be expressed as a
# single BOSS-level window and is therefore stored and applied per-tag-mode in ROOT.
tag_alg.note(:deltaE_window,
             "tag-side ΔE requirement is reconstructed-mode dependent; cannot be " \
             "expressed as a single BOSS-level window, stored and windowed per mode in ROOT")

tag_alg.apply
tag_alg.execute_on([data_3773, incMC_3773, exMC_signal])