# DSL auto-generated from 2505.12086v1
# Paper: Observation of D+ → π+ηη at BESIII
# TagAnalysis: D+ tagged (ST pattern at ψ(3770)) with 6 D- tag modes
# Signal side: π+ + 2η(→γγ each) = 1 charged pion + 4 photons

### Dataset preparation ###
data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

decay_card = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 pi+ eta eta PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Dplus_pi_eta_eta"
  config.related_dataset = data_3773
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### TagAnalysis: tag D- (anti-D+) via 6 modes, fully reconstruct D+ → π+ηη on signal side ###
alg = TagAnalysis.new("DplusPiEtaEta")
alg.set_header(["DplusPiEtaEtaAlg/DplusPiEtaEta.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })
    .note(:bdts, "BDTG classifier trained on kinematic variables used for background suppression; applied in ROOT analysis")
    .note(:tag_windows, "Tag ST windows: M_BC in [1.860, 1.880] GeV, |ΔE| < 0.040 GeV; applied in ROOT")
    .note(:dplus_mass_constraint, "M(π+ηη) constrained to D+ nominal mass in kinematic fit; not expressible in TagAnalysis DSL fit block (eta intermediate not in derived participant set)")
    .with_decay_card(decay_card)

# Tag side: D- (anti-D+, charm -1) with 6 hadronic modes
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

# Signal side: 1 π+ + 4 photons (2η → 4γ)
alg.signal_side do |s|
  s.photons 4
  s.charged pip: 1
  s.min_photon_angle 10.0
end

# Kinematic fit: 4C + two η mass constraints
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_signal])