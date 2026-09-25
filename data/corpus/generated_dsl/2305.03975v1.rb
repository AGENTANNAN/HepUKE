### Dataset preparation ###
# Real data at 3.773 GeV (psi(3770)) and the corresponding inclusive MC
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the correlated signal: psi(3770) -> D0 D0bar,
# signal D0 -> K_S0 pi+ pi- pi0 (K_S0 -> pi+ pi-, pi0 -> gamma gamma),
# opposite-side (tag) D0bar -> K+ pi- flavour mode.
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_S0 pi+ pi- pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive signal MC
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKsPiPiPi0"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "D0CPFracKsPiPiPi0"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})

# Tag side: single tag from the pre-stored D-tag collection (flavour mode D0bar -> K+ pi-).
# Only the K+ pi- flavour tag is available in the current tag-side implementation.
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi
  t.charm -1
end

# Signal side: what the tag did not use.
#   two photons (pi0 -> gamma gamma) + the charged pions of the CP mode
#   (K_S0 -> pi+ pi- plus the D0 daughters pi+ pi-), net charge zero.
alg.signal_side do |s|
  s.photons 2
  s.charged(pip: 2, pim: 2)
  s.require_charge 0
end

# Kinematic fit: 4-momentum constraint + pi0 mass constraint on the two photons.
# (pi+ pi- invariant-mass window for the K_S0 candidate; the K_S0 secondary-vertex
#  fit itself is not expressible here and is captured in the notes below.)
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:pip, :pim).between(0.487, 0.511)
  f.chi2_cut 200
end

# Inexpressible BOSS-side procedures preserved for the systematic-uncertainty step.
alg.note(:ks0_vertex_fit, "K_S0 -> pi+ pi- secondary-vertex fit with decay-length
  significance cut L/sigma_L > 2 and a pi+ pi- invariant-mass window around the
  K_S0 nominal mass; not expressible in the tag-side fit block and applied on the
  reconstructed signal side.")
   .note(:mbc_deltae_stored, "tag-side M_BC and deltaE are stored unconditionally
  (store-not-cut); single-tag and double-tag M_BC / deltaE windows are applied later
  in the ROOT analysis, not in BOSS.")
   .note(:tag_mode_restriction, "only the D0bar -> K+ pi- flavour tag is available in
  the current tag-side implementation; the CP-eigenstate and quasi-CP tag modes
  (Table I) are applied in the ROOT analysis rather than in this BOSS code.")

alg.with_decay_card(decay_card_signal).apply

# Execute on real data, inclusive MC and the exclusive signal MC
root_files = alg.execute_on([data_3773, incMC_3773, exMC_signal])