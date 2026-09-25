# Core DSL classes and dependencies are loaded automatically at execution.
### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")      # psi(3770) real data (BOSS 712, 3.773 GeV)
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC sample

# Decay card for psi(3770) -> D0 anti-D0, with D0 -> a0(980)- e+ nu_e (signal)
# and anti-D0 -> K+ pi- (tag side).
decay_card_d0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 a0(980)- e+ nu_e PHSP;
    Enddecay

    Decay a0(980)-
    1.000 eta pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for psi(3770) -> D+ D-, with D+ -> a0(980)0 e+ nu_e (signal)
# and D- -> K+ pi- pi- (tag side).
decay_card_dp = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 a0(980)0 e+ nu_e PHSP;
    Enddecay

    Decay a0(980)0
    1.000 eta pi0 PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples, 200k events per signal mode
exMC_d0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0_a0ev"
  config.related_dataset = data_3773
  config.events          = 200_000
  config.decay_card      = decay_card_d0
  config.cross_section   = :default
end

exMC_dp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dplus_a0ev"
  config.related_dataset = data_3773
  config.events          = 200_000
  config.decay_card      = decay_card_dp
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
# ------------------------------------------------------------------
# Mode I: D0 -> a0(980)- e+ nu_e  (tag the opposite-charm anti-D0)
# ------------------------------------------------------------------
alg_d0 = TagAnalysis.new("D0A0EvTag")
alg_d0.set_header(["D0A0EvTagAlg/D0A0EvTag.h"])
      .set_constant({"ECMS" => [:double, 3.773]})       # 3.773 GeV
      .with_decay_card(decay_card_d0)

# Tag side: anti-D0 -> K+ pi-, K+ pi- pi0, K+ pi- pi+ pi- ; DeltaE in [-64, +35] MeV
alg_d0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.window :deltaE, min: -0.064, max: 0.035
end

# Signal side: one pi- + one e+ + two photons (from eta), net charge 0, missing nu_e
alg_d0.signal_side do |s|
  s.photons 2
  s.charged(pim: 1, ep: 1)
  s.require_charge 0
  s.missing :nu_e
end

# Kinematic fit: 4C + constrain the photon pair to the nominal eta mass, chi2 < 200
alg_d0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_d0.note(:pid_correction_method, "Electron PID uses the tag-framework fixed v1 thresholds; " \
                                    "intended criteria E/pc > 0.8 and combined e-PID ratio > 0.8 " \
                                    "are not DSL-tunable. pi/K separation is implicit in the tag modes.")
      .note(:background_veto, "K_L0 veto via lateral-moment cut (0, 0.35) on the higher-energy " \
                              "photon from eta.")
      .note(:tag_framework_defaults, "Explicit |cos(theta)|, |Vz|, Vr and photon-energy thresholds " \
                                     "follow the standard tag-framework selections.")

alg_d0.apply

# ------------------------------------------------------------------
# Mode II: D+ -> a0(980)0 e+ nu_e  (tag the opposite-charm D-)
# ------------------------------------------------------------------
alg_dp = TagAnalysis.new("DplusA0EvTag")
alg_dp.set_header(["DplusA0EvTagAlg/DplusA0EvTag.h"])
      .set_constant({"ECMS" => [:double, 3.773]})       # 3.773 GeV
      .with_decay_card(decay_card_dp)

# Tag side: D- -> K+ pi- pi-, K+ pi- pi+ pi0, K_S0 pi-, K_S0 pi- pi0, K_S0 pi- pi+ pi-,
# K+ K- pi- ; DeltaE in [-60, +34] MeV
alg_dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.window :deltaE, min: -0.060, max: 0.034
end

# Signal side: one e+ + four photons (eta -> gamma gamma, pi0 -> gamma gamma),
# net charge 1, missing nu_e
alg_dp.signal_side do |s|
  s.photons 4
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

# Kinematic fit: 4C + constrain the photon pair to the nominal eta mass and an
# additional photon pair to the nominal pi0 mass, chi2 < 200
alg_dp.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_dp.note(:pid_correction_method, "Electron PID uses the tag-framework fixed v1 thresholds; " \
                                    "intended criteria E/pc > 0.8 and combined e-PID ratio > 0.8 " \
                                    "are not DSL-tunable. pi/K separation is implicit in the tag modes.")
      .note(:background_veto, "K_L0 veto via lateral-moment cut (0, 0.35) on the higher-energy " \
                              "photon from eta; the best a0 combination is chosen by the smallest " \
                              "sum of the pi0 and eta 1C chi2 in ROOT.")
      .note(:tag_framework_defaults, "Explicit |cos(theta)|, |Vz|, Vr and photon-energy thresholds " \
                                     "follow the standard tag-framework selections.")

alg_dp.apply

# Execute the tag analyses on data, inclusive MC and the corresponding signal MC
root_files_d0 = alg_d0.execute_on([data_3773, incMC_3773, exMC_d0])
root_files_dp = alg_dp.execute_on([data_3773, incMC_3773, exMC_dp])