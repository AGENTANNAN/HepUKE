# ============================================================================
# ψ(3770) dineutrino search:  e+e- → D0 D0bar,  D0 → π0 ν νbar  (double-tag)
# Tag side  : anti-D0 (charm = -1) in three hadronic modes
# Signal side: all-neutral D0 → π0 ν νbar (π0 → γγ) + missing neutrino
# ============================================================================

### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")      # ψ(3770) real data at √s = 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC

# Decay card for the process under study: ψ(3770) → D0 D0bar (EvtGen format)
decay_card = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC: 300k events of ψ(3770) → D0 D0bar
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3770_D0D0bar"
  config.related_dataset = data_3773          # tie the MC to the 3.773 GeV real data
  config.events          = 300_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
alg_name = "D0ToPi0NuNuBar"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(decay_card)

# Tag side: reconstruct the anti-D0 (charm = -1) in the three hadronic modes.
# D0 → K−π+,  D0 → K−π+π0,  D0 → K−π+π+π−
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: the all-neutral D0 → π0 ν νbar opposite the tag.
# π0 → γγ gives exactly two tag-unused photons; zero additional charged tracks is
# the default (no charged declaration); one massless neutrino via the recoil method.
alg.signal_side do |s|
  s.photons 2                 # the γγ pair from π0 → γγ not used by the tag
  s.missing :nu               # one massless missing neutrino (recoil reconstruction)
end

# Kinematic fit: 4-momentum conservation + γγ invariant-mass constraint to m(π0)
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

# Physics-motivated veto with no dedicated DSL primitive: the missing-mass-squared
# cut that rejects D0 → K_L^0 π0 and D0 → Kbar*(892)^0 π0 backgrounds. The tag layer
# stores M²_miss unconditionally (store-not-cut), so the window is applied at ROOT level.
alg.note(:background_veto,
  "Require M^2_miss = (p_e+e- - p_tag - p_pi0)^2 in (1.1, 1.9) GeV^2/c^4 to reject " \
  "D0 -> K_L^0 pi0 and D0 -> Kbar*(892)^0 pi0 backgrounds; M^2_miss is stored by the " \
  "tag layer and windowed in the ROOT analysis.")

# Render the tag specification (no Selection argument for a tag analysis)
alg.apply

### Execute on real data, inclusive MC and signal exclusive MC ###
root_files = alg.execute_on([data_3773, incMC_3773, exMC_signal])