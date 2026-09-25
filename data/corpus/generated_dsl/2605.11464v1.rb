### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")      # psi(3770) real data (sqrt(s) = 3.773 GeV)
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")    # corresponding inclusive MC sample

# Decay card for the signal process e+e- -> D+ D-,
# signal side D+ -> K_S0 K_L0 pi+ with K_S0 -> pi+ pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K_S0 K_L0 pi+ PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal channel: 500k events of D+ -> K_S0 K_L0 pi+, K_S0 -> pi+ pi-
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dp_KsKlpi"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (TagAnalysis) ###
alg_name = "DpTagKsKlpi"
tag_alg = TagAnalysis.new(alg_name)
tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.773]})
       .with_decay_card(decay_card_signal)

# Tag side: single-tag D- reconstructed in the three hadronic modes
#   D- -> K+ pi- pi-, K_S0 pi-, K_S0 pi+ pi- pi-
tag_alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKsPiPiPi   # DTagAlg channel names (charge-conjugate of D+ modes)
  t.charm -1                                     # pin the reconstructed (anti-)D state to D-
  t.window :mBC, min: 1.863, max: 1.877          # beam-constrained mass window (GeV/c^2)
end

# Signal side: D+ -> K_S0 K_L0 pi+ ; K_S0 -> pi+ pi- ; K_L0 is the one missing particle
tag_alg.signal_side do |s|
  s.charged(pip: 2, pim: 1)    # exactly 2 pi+ and 1 pi- (pi+ from D+ plus the K_S0 daughters)
  s.require_charge 1           # net charge +1
  s.missing :K_L0              # K_L0 treated as missing (massive missing-particle form)
  s.min_photon_angle 10.0      # minimum photon angle to charged tracks (degrees)
end

# Kinematic fit: four-momentum conservation + K_S0 mass constraint on the pi+pi- pair
tag_alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut 100
end

# BOSS-side procedures with no dedicated DSL construct
tag_alg
  .note(:background_veto,
        "reconstructed pi0 and eta -> gamma gamma candidates are vetoed to suppress the "
        "D+ -> K_S0 pi+ eta and D+ -> K_S0 K_S0 pi+ backgrounds")
  .note(:second_vertex,
        "K_S0 reconstructed from an oppositely charged pion pair with the standard BESIII "
        "secondary-vertex fit prior to the kinematic fit")
  .note(:tag_best_candidate,
        "best single-tag D- candidate per event chosen by minimum |Delta E|")

tag_alg.apply
tag_alg.execute_on([data_3773, incMC_3773, exMC_signal])