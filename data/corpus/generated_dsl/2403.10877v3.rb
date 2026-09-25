# Core DSL classes and dependencies will be loaded automatically at execution
### Datasets ###
# psi(3770) real data and inclusive MC (BOSS 712, sqrt(s) = 3.773 GeV)
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the signal process (EvtGen format, EvtGen particle names)
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi0 mu+ nu_mu PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 500k events from psi(3770) -> D0 anti-D0 (phase space)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0D0bar_signal"
  config.related_dataset = data_3773                      # associated real dataset
  config.events          = 500000                         # 500k-event exclusive MC
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
alg_name = "D0ToKPi0MuNu"
tag_alg  = TagAnalysis.new(alg_name)
tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })      # sqrt(s) = 3.773 GeV
       .with_decay_card(decay_card_signal)

# Tag side: opposite-side anti-D0 tag (charm -1) in the six hadronic tag modes.
# One tag_side declaration -> single-tag pattern; mBC/deltaE are stored, not cut.
tag_alg.tag_side(:D0) do |t|
  t.modes :D0toKPi,         # K+ pi-
          :D0toKPiPi0,       # K+ pi- pi0
          :D0toKPiPiPi,      # K+ pi+ pi- pi-
          :D0toKsPiPi,       # K_S0 pi+ pi-
          :D0toKPiPiPi0,     # K+ pi- pi- pi0
          :D0toKPiPiPiPi0    # K+ pi+ pi- pi- pi0
  t.charm -1                 # pin the tagged side to the anti-D0
end

# Signal side: the D0 -> K- pi0 mu+ nu_mu content not used by the tag.
tag_alg.signal_side do |s|
  s.photons 2                # two photons from pi0 -> gamma gamma
  s.charged(km: 1, mup: 1)   # one K- and one mu+
  s.require_charge 0         # net charge zero on the signal side
  s.missing :nu_mu           # unconstrained (massless) missing nu_mu
end

# Kinematic fit: 4-momentum conservation over tag + signal + neutrino at the CMS
# energy, pi0 nominal-mass constraint on the two photons, chi2 < 200.
tag_alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

# Render the tag specification (no Selection argument for a TagAnalysis)
tag_alg.apply

# Apply the same selection to real data, inclusive MC and the exclusive signal MC
root_files = tag_alg.execute_on([data_3773, incMC_3773, exMC_signal])