### Dataset description ###
data_3773 = DatasetManager.real_data.find("712_3773")        # ψ(3770) real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")    # corresponding inclusive MC

# Decay card: ψ(3770) → D0 anti-D0, both decaying to the same CP eigenstate K_S0 π0
decay_card_dt = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K_S0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.000 K_S0 pi0 PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive signal MC events for the double-tag channel
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "DT_D0toKsPi0_3773"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_dt
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
# Double-tag: one D0 (charm +1) tag and one anti-D0 (charm -1) tag, both in K_S0 π0.
alg_name = "D0ToKsPi0DT"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(decay_card_dt)

# Tag side 1: D0 (charm +1) → K_S0 π0
alg.tag_side(:D0) do |t|
  t.modes :D0toKsPi0
  t.charm 1
end

# Tag side 2: anti-D0 (charm -1) → K_S0 π0, second tag ranked by invariant mass
alg.tag_side(:D0) do |t|
  t.modes :D0toKsPi0
  t.charm -1
  t.rank_by :inv     # rank the second tag by invariant mass
end

# Kinematic fit: 4C energy-momentum conservation with both tag masses
# constrained to the nominal D0 mass; χ² < 200 (tight cut applied in ROOT).
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# BOSS-side requirement not expressible in the tag DSL: no additional photons
# left on the signal side beyond those used by the two D tags.
alg.note(:extra_photon_veto,
         "require no additional photons on the signal side opposite the two D0 tags " \
         "(no leftover good showers unused by the double-tag reconstruction)")

# Render the spec (no Selection argument for a TagAnalysis) and run on the datasets.
alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_signal])