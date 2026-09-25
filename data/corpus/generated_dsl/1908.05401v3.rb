### Dataset preparation ###
psi3770_data  = DatasetManager.real_data.find("712_3773")      # psi(3770) real data
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC

# Decay card for the representative charged signal:
#   psi(3770) -> D+ D- ,  D+ -> phi pi+ ,  phi -> K+ K-
decay_card_dplus = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 phi pi+ PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Decay card for the representative neutral signal:
#   psi(3770) -> D0 anti-D0 ,  D0 -> phi K- pi+ ,  phi -> K+ K-
decay_card_d0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 phi K- pi+ PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for each representative signal (efficiency reference for inclusive phi X)
exMC_dplus = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DplusToPhiPi"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_dplus
  config.cross_section   = :default
end

exMC_d0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0ToPhiKPi"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_d0
  config.cross_section   = :default
end

### Tag analysis 1 — D+ -> phi X : tag the D- hadronically ###
alg_dplus_name = "DplusPhiXTag"
alg_dplus = TagAnalysis.new(alg_dplus_name)
alg_dplus.set_header(["#{alg_dplus_name}Alg/#{alg_dplus_name}.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_dplus)

# Tag the D- in the hadronic modes K+pi-pi-, K+pi-pi-pi0, K_S pi-, K_S pi- pi0, K_S pi- pi- pi-
alg_dplus.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm -1        # pin the tagged side to the D-
end

# Signal side: the inclusive phi X  ->  require exactly the phi -> K+ K- pair
alg_dplus.signal_side do |s|
  s.charged(kp: 1, km: 1)   # exactly one K+ and one K-
  s.require_charge 0        # net charge zero
end

# Four-momentum-constrained kinematic fit over tag + signal
alg_dplus.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_dplus.note(:quantum_correlation_correction,
  "quantum-correlation correction factor taken as unity for the charged D+ (D-) tag sample")

alg_dplus.apply    # tag analysis takes no Selection argument
root_files_dplus = alg_dplus.execute_on([psi3770_data, psi3770_incMC, exMC_dplus])

### Tag analysis 2 — D0 -> phi X : tag the Dbar0 hadronically ###
alg_d0_name = "D0PhiXTag"
alg_d0 = TagAnalysis.new(alg_d0_name)
alg_d0.set_header(["#{alg_d0_name}Alg/#{alg_d0_name}.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_d0)

# Tag the Dbar0 in the hadronic modes K+pi-, K+pi-pi0, K+pi-pi-pi+
alg_d0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1        # pin the tagged side to the Dbar0
end

# Signal side: the inclusive phi X  ->  require exactly the phi -> K+ K- pair
alg_d0.signal_side do |s|
  s.charged(kp: 1, km: 1)   # exactly one K+ and one K-
  s.require_charge 0        # net charge zero
end

# Four-momentum-constrained kinematic fit over tag + signal
alg_d0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0.note(:quantum_correlation_correction,
  "quantum-correlation correction factor applied for the neutral D0 (Dbar0) tag sample, taken from external inputs; carried outside the BOSS selection")

alg_d0.apply    # tag analysis takes no Selection argument
root_files_d0 = alg_d0.execute_on([psi3770_data, psi3770_incMC, exMC_d0])