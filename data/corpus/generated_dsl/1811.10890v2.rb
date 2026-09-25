# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
data_4180  = DatasetManager.real_data.find("703_4180")     # 4.178 GeV real data (BOSS 703/4180)
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")  # corresponding inclusive MC sample

# Decay card for the signal chain (EvtGen syntax):
#   psi(4260) -> Ds+ Ds*-,  Ds*- -> gamma Ds-,  Ds+ -> mu+ nu_mu
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s+ D_s*- PHSP;
    Enddecay

    Decay D_s*-
    1.000 gamma D_s- PHSP;
    Enddecay

    Decay D_s+
    1.000 mu+ nu_mu PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal process (100k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "DsToMuNu_signal_mc"
  config.related_dataset = data_4180    # tie the signal MC to the 4.178 GeV real data
  config.events          = 100_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) — tag-based analysis ###
alg_name = "DsToMuNuTag"
ds_tag = TagAnalysis.new(alg_name)
ds_tag.set_header(["#{alg_name}Alg/#{alg_name}.h"])
      .set_constant({"ECMS" => [:double, 4.178]})   # 4.178 GeV centre-of-mass energy
      .with_decay_card(decay_card_signal)

# BOSS-side procedures that have no formal DSL construct are preserved as notes.
ds_tag.note(:muon_pid_selection, "signal-side mu+ identified by muon PID using EMC energy in (0.0, 0.3) GeV together with a 2D muon-counter hit-depth requirement; these thresholds live in the generated signal-side PID recipe (fixed v1 defaults) and are not DSL-tunable")
      .note(:background_veto, "events with unused EMC shower energy above 0.4 GeV on the signal side are vetoed to suppress gamma/pi0 and split-off backgrounds")

# Tag side: single-tag Ds- (negative charm), reconstructed in 13 hadronic modes.
# The soft gamma/pi0 from Ds*- and the Ds* mass constraint are handled internally
# by the DTagAlg through the minimum |DeltaE| criterion, so no soft-photon
# reconstruction step is expressed here.
ds_tag.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKKPiPi0, :DstoKsK, :DstoKsKPi0, :DstoKsKsPi,
          :DstoKsKPiPi, :DstoKPiPi, :DstoPiPiPi, :DstoPiEta, :DstoPiPi0Eta,
          :DstoPiEtaPrime, :DstoPiOmega, :DstoKPiPiPi0
  t.charm -1                                  # require the negative-charm Ds- tag
  t.window :mBC,    min: 2.010, max: 2.073    # ST tag M_BC window (GeV/c^2)
  t.window :deltaE, min: -0.05, max: 0.10     # DT signal DeltaE window (GeV)
end

# Signal side: the tracks left over by the tag must give exactly one mu+ (net
# charge +1) and the event must contain one missing massless nu_mu.
ds_tag.signal_side do |s|
  s.charged(mup: 1)      # exactly one mu+ -> any extra good charged track is vetoed
  s.require_charge(1)    # net charge of the signal side = +1
  s.missing :nu_mu       # massless missing muon neutrino
end

# 4C kinematic fit: constrain the final state to the CMS four-momentum and the
# tagged Ds to its nominal mass (the Ds* mass is constrained internally by the tag).
ds_tag.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"D_s-")
  f.chi2_cut 200
end

# Validate/render the tag spec, then apply the same selection chain to data,
# inclusive MC and signal MC.
ds_tag.apply
root_files = ds_tag.execute_on([data_4180, incMC_4180, exMC_signal])