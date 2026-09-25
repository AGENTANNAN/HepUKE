### Dataset preparation ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data at 3.773 GeV
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC sample

# Decay card for the double-tag signal channel:
#   ψ(3770) -> D+ D-;  tag D- -> K+ π- π-;  signal D+ -> K1(1270)0 e+ ν_e (ISGW2),
#   with the K1(1270)0 described by a relativistic Breit-Wigner in the K- π+ π0 channel.
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay D+
    1.0000 K1(1270)0 e+ nu_e ISGW2;
    Enddecay

    Decay K1(1270)0
    1.0000 K- pi+ pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 200k exclusive signal-MC events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_K1enu_DT"
  config.related_dataset = psi3770_data
  config.events          = 200_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) — tag-based, double-tag technique ###
alg = TagAnalysis.new("K1enuDTag")
alg.set_header(["K1enuDTagAlg/K1enuDTag.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(decay_card_signal)
   # BOSS-side procedures with no dedicated DSL construct
   .note(:positron_selection, "positron identified by CLe > 0.001, CLe/(CLe+CLpi+CLK) > 0.8 and E/p > 0.8")
   .note(:pi0_selection, "pi0 -> gamma gamma candidate required to have momentum > 0.15 GeV/c and |cos(theta_decay)| < 0.8")
   .note(:fsr_recovery, "FSR / bremsstrahlung photons within 5 degrees of the positron are recovered and included in the kinematic fit")

# Tag side: the D- is reconstructed in six hadronic modes (charm -1 pins the D-)
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

# Signal side: everything the tag did not use
alg.signal_side do |s|
  s.photons 2                      # two photons for pi0 -> gamma gamma
  s.min_photon_angle 10.0          # photon opening angle > 10 degrees
  s.charged(km: 1, pip: 1, ep: 1)  # K-, pi+, e+
  s.missing :nu_e                  # missing neutrino (massless)
end

# 5C kinematic fit: 4-momentum conservation + gamma gamma mass constrained to pi0
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)               # pi0 mass window
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # 4C + 1C = 5C
  f.chi2_cut 200
end

alg.apply
alg.execute_on([psi3770_data, psi3770_incMC, exMC_signal])