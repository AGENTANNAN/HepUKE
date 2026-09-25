# BESIII search for D+ → γρ+ and D+ → γK*+ using DT method
# ArXiv: 2410.06500v1
# Dataset: 20.3 fb^-1 at √s = 3.773 GeV (ψ(3770))
# Two signal modes share the tag side; separate TagAnalysis objects per mode

### Dataset preparation ###
psip3770_data = DatasetManager.real_data.find("712_3773")
psip3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card: ψ(3770) → D+ D-; signal D+ → γρ+, ρ+ → π+π0
decay_card_rho = <<~DECAYCARD
    Decay psi(3770)
    1.0000  D+  D-  PHSP;
    Enddecay

    Decay D+
    1.0000  gamma  rho+  PHSP;
    Enddecay

    Decay rho+
    1.0000  pi+  pi0  VSS;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma  PHSP;
    Enddecay

    Decay D-
    1.0000  K+  pi-  pi-  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: ψ(3770) → D+ D-; signal D+ → γK*+, K*+ → K+π0
decay_card_kstar = <<~DECAYCARD
    Decay psi(3770)
    1.0000  D+  D-  PHSP;
    Enddecay

    Decay D+
    1.0000  gamma  K*+  PHSP;
    Enddecay

    Decay K*+
    1.0000  K+  pi0  VSS;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma  PHSP;
    Enddecay

    Decay D-
    1.0000  K+  pi-  pi-  PHSP;
    Enddecay

    End
DECAYCARD

exMC_rho = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_gamma_rho"
  config.related_dataset = psip3770_data
  config.events = 1_000_000
  config.decay_card = decay_card_rho
  config.cross_section = :default
end

exMC_kstar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_gamma_kstar"
  config.related_dataset = psip3770_data
  config.events = 1_000_000
  config.decay_card = decay_card_kstar
  config.cross_section = :default
end

### Mode I: D+ → γρ+, ρ+ → π+π0 ###
alg_rho = TagAnalysis.new("Dp2GammaRho")
alg_rho.set_header(["Dp2GammaRhoAlg/Dp2GammaRho.h"])
  .set_constant(ECMS: 3.773)
  .with_decay_card(decay_card_rho)

alg_rho.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKKPi,
          :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm -1
end

alg_rho.signal_side do |s|
  s.photons 3
  s.min_photon_angle 10.0
  s.charged(pip: 1)
end

alg_rho.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_rho.apply
alg_rho.execute_on([psip3770_data, psip3770_incMC, exMC_rho])

### Mode II: D+ → γK*+, K*+ → K+π0 ###
alg_kstar = TagAnalysis.new("Dp2GammaKstar")
alg_kstar.set_header(["Dp2GammaKstarAlg/Dp2GammaKstar.h"])
  .set_constant(ECMS: 3.773)
  .with_decay_card(decay_card_kstar)

alg_kstar.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKKPi,
          :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm -1
end

alg_kstar.signal_side do |s|
  s.photons 3
  s.min_photon_angle 10.0
  s.charged(kp: 1)
end

alg_kstar.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_kstar.apply
alg_kstar.execute_on([psip3770_data, psip3770_incMC, exMC_kstar])