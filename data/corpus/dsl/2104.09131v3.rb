# Amplitude analysis of D+ -> K+ K_S0 pi0 at psi(3770)
# Paper: 2104.09131v3
# Tag-based double-tag analysis (DTagAlg)
# Signal: D+ -> K+ K_S0 pi0, K_S0 -> pi+ pi-, pi0 -> gamma gamma
# Tag: D- with 6 hadronic tag modes

alg = TagAnalysis.new('Dp_to_KpKsPi0_DT')
alg.set_header(['EmcRecEventTag/EmcRecEventTag.h'])
alg.set_constant(ECMS: 3.773)

alg.with_decay_card <<~DECAY
  Decay psi(3770)
  1.000 D+ D- PHSP;
  Enddecay
  Decay D+
  1.000 K+ K_S0 pi0 PHSP;
  Enddecay
  Decay D-
  1.000 tag hadronic;
  Enddecay
  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
DECAY

# Tag side 1: D- with 6 hadronic modes
alg.tag_side(:Dm) do |t|
  t.modes :DmtoKpPimPim, :DmtoKpPimPimPi0, :DmtoKsPim,
          :DmtoKsPimPi0, :DmtoKsPimPimPip, :DmtoKpKmPim
  t.charm(-1)
end

# Tag side 2: D+ with 1 hadronic tag mode (D+ -> K- pi+ pi+)
alg.tag_side(:Dp) do |t|
  t.modes :DptoKmPipPip
  t.charm(1)
  t.rank_by :inv
end

# Signal side: K+ K_S0 pi0
# K_S0 -> pi+ pi- reconstructed; pi0 -> gamma gamma
alg.signal_side do |s|
  s.charged(kp: 1)
  s.photons 2
  s.min_photon_angle 10.0
  s.require_charge 1      # K+ has +1
end

# Kinematic fit: 4C + pi0 mass constraint + tag-side D mass constraints (6C)
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Dplus)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:Dplus)
  f.chi2_cut 200
end

alg.apply

# Exclusive MC
DatasetManager.create_exclusive_mc do |mc|
  mc.decay_card <<~DECAY
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay
    Decay D+
    1.000 K+ K_S0 pi0 PHSP;
    Enddecay
    Decay D-
    1.000 tag hadronic;
    Enddecay
    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
  DECAY
  mc.related_dataset DatasetManager.load_real_data.find('712_3773')
  mc.n_events 2_000_000
end

# Execute on psi(3770) data
alg.execute_on([
  DatasetManager.load_real_data.find('712_3773'),
  DatasetManager.load_inclusive_mc.find('712_3773')
])

alg.note(:post_kinematic_cuts,
  'K_S0 mass window, pi0 mass window, flight distance significance cuts, ' \
  'and final kinematic fit chi2 optimization are applied in ROOT analysis.')
alg.note(:amplitude_analysis,
  'The amplitude analysis of D+ -> K+ K_S0 pi0 is performed in ROOT stage, ' \
  'including the Dalitz plot fit and systematic uncertainty evaluation.')