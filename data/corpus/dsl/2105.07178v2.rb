# Ds+ -> tau+ nu_tau at 6 energy points
# Paper: 2105.07178v2
# Tag-based ST+missing analysis
# 14 Ds- hadronic tag modes, signal: tau+ -> pi+ pi0 anti-nu_tau

alg = TagAnalysis.new('Dsp_to_taunu_ST')
alg.set_header(['EmcRecEventTag/EmcRecEventTag.h'])
alg.set_constant(ECMS: 4.178)

alg.with_decay_card <<~DECAY
  Decay anti-D_s-
  1.000 tag hadronic;
  Enddecay
  Decay D_s+
  1.000 tau+ nu_tau PHSP;
  Enddecay
  Decay tau+
  1.000 pi+ pi0 anti-nu_tau PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
DECAY

# Tag side: Ds- with 14 hadronic modes
alg.tag_side(:Ds) do |t|
  t.mode_group :hadronic
  t.charm(-1)
end

# Signal side: pi+ pi0 + missing nu_tau
alg.signal_side do |s|
  s.charged(pip: 1)
  s.photons 2
  s.min_photon_angle 10.0
  s.missing :nu_tau         # massless (tau neutrino)
  s.require_charge 1        # pi+ has +1
end

# Kinematic fit: 4C + pi0 mass constraint
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg.apply

# Exclusive MC
DatasetManager.create_exclusive_mc do |mc|
  mc.decay_card <<~DECAY
    Decay anti-D_s-
    1.000 tag hadronic;
    Enddecay
    Decay D_s+
    1.000 tau+ nu_tau PHSP;
    Enddecay
    Decay tau+
    1.000 pi+ pi0 anti-nu_tau PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
  DECAY
  mc.related_dataset DatasetManager.load_real_data.find('703_4180')
  mc.n_events 2_000_000
end

# Execute on all 6 energy points
alg.execute_on([
  DatasetManager.load_real_data.find('703_4180'),   # 4.178 GeV, 3189 pb^-1
  DatasetManager.load_real_data.find('703_4190'),   # 4.189 GeV
  DatasetManager.load_real_data.find('703_4200'),   # 4.199 GeV
  DatasetManager.load_real_data.find('703_4210'),   # 4.209 GeV
  DatasetManager.load_real_data.find('703_4220'),   # 4.219 GeV
  DatasetManager.load_real_data.find('703_4230'),   # 4.226 GeV
  DatasetManager.load_inclusive_mc.find('703_4180')
])

alg.note(:tag_modes,
  '14 Ds- hadronic tag modes: Ds- -> K_S0 K-, K_S0 K- pi+ pi-, ' \
  'K+ K- pi-, K+ K- pi- pi0, K_S0 K+ pi- pi-, pi+ pi+ pi- pi- pi-, ' \
  'K+ K+ K- pi- pi-, K_S0 K- pi0, K+ K- pi- pi+ pi-, ' \
  'K_S0 K+ pi- pi- pi0, and others.')
alg.note(:post_kinematic_cuts,
  'MM^2 signal/background discrimination, E_extra_gamma cuts, ' \
  'pi0 mass window, and final fit chi2 optimization are applied in ROOT analysis.')
alg.note(:multi_energy,
  'Combined analysis over 6 energy points. Branching fractions measured ' \
  'by fitting the MM^2 distributions at each energy point simultaneously.')