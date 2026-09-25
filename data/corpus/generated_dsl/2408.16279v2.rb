### Dataset description ###
# ψ(3770) @ 3.773 GeV: real data (2.93 fb^-1) + corresponding inclusive MC
psipp_data  = DatasetManager.real_data.find("712_3773")
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the signal process: ψ(3770) → D0 anti-D0 phase space
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC: 200,000 ψ(3770) → D0 anti-D0 phase-space events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psipp_D0D0bar_phsp"
  config.related_dataset = psipp_data       # associated real dataset
  config.events          = 200000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-versus-signal analysis (BOSS) ###
alg_name = "D0StrongPhase"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})   # ψ(3770) CMS energy (GeV)
   .with_decay_card(decay_card_signal)

# Tag side: flavor-tag D0 in the three hadronic modes
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi   # D0 → K-π+, K-π+π0, K-π+π-π+
  # no charm pinned -> both D0 / anti-D0 flavor tags scanned
  # no window declared -> tag mBC / ΔE stored unconditionally, windowed later in ROOT
end

# Signal side: D0 → π+π-π+π- built from the remaining tracks
alg.signal_side do |s|
  s.charged(pip: 2, pim: 2)   # two π+ and two π-
  s.require_charge(0)         # net charge zero; content fixed only by the tag modes + this multiset
end

# 4C kinematic fit: constrain tag+signal four-momentum and the tag D0 to its nominal mass
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# Semileptonic / KL0 tag modes and the ROOT-level corrections belong to separate handling
alg.note(:semileptonic_tag_modes,
  "Semileptonic K- e+ nu_e and K_L0 tag modes use partial reconstruction with " \
  "missing-mass variables analysed in ROOT; they are a separate selection and are " \
  "not part of this hadronic-tag selection.")

# Render the tag spec (takes no Selection argument) and run on the datasets
alg.apply
root_files = alg.execute_on([psipp_data, psipp_incMC, exMC_signal])