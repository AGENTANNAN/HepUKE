# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")      # ψ(3770) real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC sample

# Decay card for the DCS signal mode D0 → K+π−π0
# (ψ(3770) → D0 anti-D0; signal D0 → K+π−π0; recoiling anti-D0 → K+π−; π0 → γγ; all phase space)
decay_card_signal_pi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K+ pi- pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the DCS signal mode D0 → K+π−π0π0
# (ψ(3770) → D0 anti-D0; signal D0 → K+π−π0π0; recoiling anti-D0 → K+π−; π0 → γγ; all phase space)
decay_card_signal_pi0pi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K+ pi- pi0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# One exclusive MC sample per signal mode
exMC_signal_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKPiPi0"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_signal_pi0
  config.cross_section   = :default
end

exMC_signal_pi0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKPiPi0Pi0"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_signal_pi0pi0
  config.cross_section   = :default
end

### Event selection (tag-based) ###
# The two signal modes are independent decay modes (different photon multiplicity and
# different fit constraints) → two separate TagAnalysis objects (Rule T1).
# Both share the SAME hadronic tag chain: one D0 tagged via D0 → K+π−, K+π−π0, K+π−π−π+.

# ---------------- Mode I : D0 → K+π−π0 ----------------
alg_pi0 = TagAnalysis.new("D0DCSKPiPi0")
alg_pi0.set_header(["D0DCSKPiPi0Alg/D0DCSKPiPi0.h"])
       .set_constant({"ECMS" => [:double, 3.773]})   # ψ(3770) centre-of-mass energy (GeV)
       .with_decay_card(decay_card_signal_pi0)

# Tag side: one tagged D0 through its hadronic modes (single tag; tag mBC/ΔE are
# stored unconditionally and windowed later in the ROOT analysis).
alg_pi0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
end

# Signal side: D0 → K+π−π0 → K+ π− γγ
alg_pi0.signal_side do |s|
  s.charged(kp: 1, pim: 1)   # exactly one K+ and one π−
  s.photons 2                # two photons forming the π0
end

# Kinematic fit: 4-momentum constraint + γγ invariant mass constrained to the
# nominal π0 mass (this constraint is applied to the K+π−π0 channel only)
alg_pi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_pi0.apply   # tag spec renders without a Selection argument

# ---------------- Mode II : D0 → K+π−π0π0 ----------------
alg_pi0pi0 = TagAnalysis.new("D0DCSKPiPi0Pi0")
alg_pi0pi0.set_header(["D0DCSKPiPi0Pi0Alg/D0DCSKPiPi0Pi0.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .with_decay_card(decay_card_signal_pi0pi0)

# Same hadronic tag selection chain as Mode I
alg_pi0pi0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
end

# Signal side: D0 → K+π−π0π0 → K+ π− γγγγ
alg_pi0pi0.signal_side do |s|
  s.charged(kp: 1, pim: 1)   # exactly one K+ and one π−
  s.photons 4                # four photons forming the two π0's
end

# Kinematic fit: 4-momentum constraint only (no π0 mass constraint for this channel)
alg_pi0pi0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_pi0pi0.apply

### Execute on datasets ###
root_files_pi0 = alg_pi0.execute_on([data_3773, incMC_3773, exMC_signal_pi0])
root_files_pi0pi0 = alg_pi0pi0.execute_on([data_3773, incMC_3773, exMC_signal_pi0pi0])