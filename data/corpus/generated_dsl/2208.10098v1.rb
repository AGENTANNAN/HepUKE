# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")     # psi(3770) real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # tag-based reformed inclusive MC at 3.773 GeV

# Decay card of the signal process: psi(3770) -> D0 D0bar, with the signal D0 -> pi+ pi- pi+ pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi+ pi- pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

### Event selection (BOSS) -- D0 Dbar0 quantum-correlated double-tag analysis ###
# Tag-based analysis: the tag side is a CP-eigenstate D, the signal side is the
# quantum-correlated D -> pi+ pi- pi+ pi-. No Selection chain is used.
alg_name = "D0D0barCP4pi"
d0_cp4pi = TagAnalysis.new(alg_name)
d0_cp4pi.set_header(["#{alg_name}Alg/#{alg_name}.h"])
        .set_constant({"ECMS" => [:double, 3.773]})   # sqrt(s) = 3.773 GeV
        .with_decay_card(decay_card_signal)

# Tag side: one D reconstructed in a CP-eigenstate mode, taken from the pre-stored
# DTagAlg vocabulary. Both tag charges are scanned (neutral species default).
d0_cp4pi.tag_side(:D0) do |t|
  t.modes :D0toKK,        # CP-even tag: K+ K-
          :D0toKsPi0Pi0,  # CP-even tag: K_S0 pi0 pi0
          :D0toKsPi0,     # CP-odd  tag: K_S0 pi0
          :D0toKsEta,     # CP-odd  tag: K_S0 eta
          :D0toKsEtap     # CP-odd  tag: K_S0 eta'
end

# Signal side: the opposite D -> pi+ pi- pi+ pi-
d0_cp4pi.signal_side do |s|
  s.photons 0                 # zero photons on the signal side
  s.charged(pip: 2, pim: 2)   # exactly 2 pi+ and 2 pi- (all PID'd as pions); no extra charged tracks
  s.require_charge 0          # net charge of the four signal pions is zero
end

# Kinematic fit: four-momentum conservation + tag invariant mass constrained to the
# nominal D0 mass, with chi2 < 200.
d0_cp4pi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# BOSS-side requirements that have no formal DSL construct (store-not-cut tag analysis)
d0_cp4pi
  .note(:background_veto, "signal-side K_S0 veto: reject the event if any pi+pi- pair has invariant mass in [0.481, 0.514] GeV/c^2 with L/sigma_L > 2")
  .note(:signal_dE_selection, "signal-side DeltaE required within [-0.026, 0.023] GeV")
  .note(:tag_candidate_ranking, "candidates of each tag mode are ranked by invariant mass (closest to the nominal D0 mass)")
  .note(:tag_modes_not_in_dtagalg, "K_L0 pi0, K_L0 omega, pi0 pi0, K_S0 omega, K_S0 phi, K_L0 pi0 pi0, plus the quasi-CP-even pi+ pi- pi0 and the phase-space-binned K_S,L0 pi+ pi- tags, are not in the DTagAlg vocabulary and are handled in separate algorithms")

d0_cp4pi.apply
d0_cp4pi.execute_on([data_3773, incMC_3773])