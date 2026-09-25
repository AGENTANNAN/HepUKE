# Core DSL classes and dependencies are loaded automatically at execution.

### Datasets: √s = 3.773 GeV (ψ(3770)) ###
data_3773  = DatasetManager.real_data.find("712_3773")      # 3.773 GeV real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")    # corresponding inclusive MC

### Decay cards for the four representative signal modes (EvtGen syntax) ###
# The tagged Dbar (other side) is generated generically via a representative tag mode.
decay_card_d0_kkpi0pi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K+ K- pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- pi0 PHSP;
    Enddecay

    End
DECAYCARD

decay_card_d0_kskspipi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_S0 K_S0 pi+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_dp_kkpipi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 K+ K- pi+ pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_dp_kskspipi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 K_S0 K_S0 pi+ pi0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive signal MC: 500k events per representative mode ###
exMC_d0_kkpi0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_d0_kkpi0pi0"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_d0_kkpi0pi0
  config.cross_section   = :default
end

exMC_d0_kskspipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_d0_kskspipi"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_d0_kskspipi
  config.cross_section   = :default
end

exMC_dp_kkpipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_dp_kkpipi0"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_dp_kkpipi0
  config.cross_section   = :default
end

exMC_dp_kskspipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_dp_kskspipi0"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_dp_kskspipi0
  config.cross_section   = :default
end

### ------------------------------------------------------------------ ###
### Double-tag analysis: tag D0bar; signal D0 → K+K-π0π0            ###
### ------------------------------------------------------------------ ###
alg_d0_kkpi0pi0 = TagAnalysis.new("D0TagKKPi0Pi0")
alg_d0_kkpi0pi0.set_header(["D0TagKKPi0Pi0Alg/D0TagKKPi0Pi0.h"])
               .set_constant({"ECMS" => [:double, 3.773]})
               .with_decay_card(decay_card_d0_kkpi0pi0)

# Tag side: D0bar reconstructed from the pre-stored DTag collection
alg_d0_kkpi0pi0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi   # K+π-, K+π-π0, K+π-π+π-
  t.charm -1                                     # tag the D0bar (opposite flavour to signal D0)
end

# Signal side: everything the tag did not use
alg_d0_kkpi0pi0.signal_side do |s|
  s.charged(kp: 1, km: 1)     # K+ and K-
  s.photons 4                 # two π0 → 4γ
  s.require_charge 0          # net charge 0
end

# Kinematic fit: 4C + two π0 mass constraints, loose χ² < 200
alg_d0_kkpi0pi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # π0 mass constraint #1
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # π0 mass constraint #2
  f.chi2_cut 200
end

alg_d0_kkpi0pi0.apply
root_files_d0_kkpi0pi0 = alg_d0_kkpi0pi0.execute_on([data_3773, incMC_3773, exMC_d0_kkpi0pi0])

### ------------------------------------------------------------------ ###
### Double-tag analysis: tag D0bar; signal D0 → K_S0 K_S0 π+π-       ###
### ------------------------------------------------------------------ ###
alg_d0_kskspipi = TagAnalysis.new("D0TagKsKsPiPi")
alg_d0_kskspipi.set_header(["D0TagKsKsPiPiAlg/D0TagKsKsPiPi.h"])
               .set_constant({"ECMS" => [:double, 3.773]})
               .with_decay_card(decay_card_d0_kskspipi)

alg_d0_kskspipi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

alg_d0_kskspipi.signal_side do |s|
  s.charged(pip: 3, pim: 3)   # two K_S0 → π+π- plus the prompt π+π-
  s.require_charge 0          # net charge 0
end

# The two K_S0 are reconstructed with a secondary-vertex fit on the remaining
# π+π- pairs (mass window and decay-length cut) — not expressible on a tag-signal side.
alg_d0_kskspipi.note(:ks0_secondary_vertex,
  "The two K_S0 → π+π- candidates on the signal side are reconstructed from the " \
  "remaining π+π- pairs with a secondary-vertex fit; K_S0 mass window is " \
  "0.486–0.510 GeV/c² and the decay length is required to be greater than 2σ " \
  "of the vertex resolution. Not expressible in the TagAnalysis signal-side primitives.")

# Kinematic fit: 4C only (no extra mass constraint for this channel), χ² < 200
alg_d0_kskspipi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_kskspipi.apply
root_files_d0_kskspipi = alg_d0_kskspipi.execute_on([data_3773, incMC_3773, exMC_d0_kskspipi])

### ------------------------------------------------------------------ ###
### Double-tag analysis: tag D-; signal D+ → K+K-π+π0                ###
### ------------------------------------------------------------------ ###
alg_dp_kkpipi0 = TagAnalysis.new("DpTagKKPiPi0")
alg_dp_kkpipi0.set_header(["DpTagKKPiPi0Alg/DpTagKKPiPi0.h"])
              .set_constant({"ECMS" => [:double, 3.773]})
              .with_decay_card(decay_card_dp_kkpipi0)

alg_dp_kkpipi0.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0   # K+π-π-, K_S0π-, K+π-π-π0
  t.charm -1                                       # tag the D-
end

alg_dp_kkpipi0.signal_side do |s|
  s.charged(kp: 1, km: 1, pip: 1)   # K+, K-, π+
  s.photons 2                       # π0 → 2γ
  s.require_charge 1                # net charge +1
end

alg_dp_kkpipi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # π0 mass constraint
  f.chi2_cut 200
end

alg_dp_kkpipi0.apply
root_files_dp_kkpipi0 = alg_dp_kkpipi0.execute_on([data_3773, incMC_3773, exMC_dp_kkpipi0])

### ------------------------------------------------------------------ ###
### Double-tag analysis: tag D-; signal D+ → K_S0 K_S0 π+π0          ###
### ------------------------------------------------------------------ ###
alg_dp_kskspipi0 = TagAnalysis.new("DpTagKsKsPiPi0")
alg_dp_kskspipi0.set_header(["DpTagKsKsPiPi0Alg/DpTagKsKsPiPi0.h"])
                .set_constant({"ECMS" => [:double, 3.773]})
                .with_decay_card(decay_card_dp_kskspipi0)

alg_dp_kskspipi0.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0
  t.charm -1
end

alg_dp_kskspipi0.signal_side do |s|
  s.charged(pip: 3, pim: 2)   # two K_S0 → π+π- (2π+2π-) plus the prompt π+
  s.photons 2                 # π0 → 2γ
  s.require_charge 1          # net charge +1
end

# Same K_S0 secondary-vertex reconstruction as the D0 → K_S0 K_S0 π+π- channel.
alg_dp_kskspipi0.note(:ks0_secondary_vertex,
  "The two K_S0 → π+π- candidates on the signal side are reconstructed from the " \
  "remaining π+π- pairs with a secondary-vertex fit; K_S0 mass window is " \
  "0.486–0.510 GeV/c² and the decay length is required to be greater than 2σ " \
  "of the vertex resolution. Not expressible in the TagAnalysis signal-side primitives.")

alg_dp_kskspipi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # π0 mass constraint
  f.chi2_cut 200
end

alg_dp_kskspipi0.apply
root_files_dp_kskspipi0 = alg_dp_kskspipi0.execute_on([data_3773, incMC_3773, exMC_dp_kskspipi0])