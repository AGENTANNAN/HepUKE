# Core DSL classes and dependencies are loaded automatically at execution.

### Datasets ###
data_3773  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data at √s = 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # Corresponding inclusive MC sample

### Decay cards (EvtGen) — one per signal mode ###
# Signal mode I: D0 → K_S0 K+ K-  (K_S0 → π+π-)
decay_card_d0_ks_kk = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_S0 K+ K- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Signal mode II: D0 → K_L0 K+ K-
decay_card_d0_kl_kk = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_L0 K+ K- PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC — 1M events per signal mode ###
exMC_d0_ks_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_d0_ks_kk"   # D0 → K_S0 K+ K- signal MC
  config.related_dataset = data_3773
  config.events          = 1000000
  config.decay_card      = decay_card_d0_ks_kk
  config.cross_section   = :default
end

exMC_d0_kl_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_d0_kl_kk"   # D0 → K_L0 K+ K- signal MC
  config.related_dataset = data_3773
  config.events          = 1000000
  config.decay_card      = decay_card_d0_kl_kk
  config.cross_section   = :default
end

# =========================================================================
### Mode I — Dbar0 (tag) + D0 → K_S0 K+ K-  (double-tag) ###
# =========================================================================
alg_ks_kk = TagAnalysis.new("D0TagKsKK")
alg_ks_kk.set_header(["D0TagKsKKAlg/D0TagKsKK.h"])
         .set_constant({"ECMS" => [:double, 3.773]})   # √s = 3.773 GeV
         .with_decay_card(decay_card_d0_ks_kk)

# Tag side: D̄0 tagged with charm = -1 through the seven hadronic modes
alg_ks_kk.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKK, :D0toPiPi,
          :D0toKsPi0, :D0toPiPiPi0, :D0toKsPiPi
  t.charm -1
end

# Signal side: fully reconstructed K_S0 K+ K- from kp, km and a π+π- pair (net charge 0)
alg_ks_kk.signal_side do |s|
  s.charged(kp: 1, km: 1, pip: 1, pim: 1)
  s.require_charge 0
end

# Fit: 4-momentum constraint + K_S0 mass constraint on the π+π- pair
alg_ks_kk.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut 200
end

alg_ks_kk.apply    # no Selection argument for the tag layer
alg_ks_kk.execute_on([data_3773, incMC_3773, exMC_d0_ks_kk])

# =========================================================================
### Mode II — Dbar0 (tag) + D0 → K_L0 K+ K-  (double-tag, K_L0 missing) ###
# =========================================================================
alg_kl_kk = TagAnalysis.new("D0TagKlKK")
alg_kl_kk.set_header(["D0TagKlKKAlg/D0TagKlKK.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .note(:kl_emc_shower_momentum,
               "the signal-side K_L0 momentum is obtained from its associated EMC shower; "
               "in the DSL fit the K_L0 is treated as a massive missing particle with the "
               "nominal K_L0 mass and floated momentum, instead of using the EMC-shower "
               "momentum direction as the seed")
         .with_decay_card(decay_card_d0_kl_kk)

# Tag side: same D̄0 (charm = -1) hadronic tag modes
alg_kl_kk.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKK, :D0toPiPi,
          :D0toKsPi0, :D0toPiPiPi0, :D0toKsPiPi
  t.charm -1
end

# Signal side: one K+, one K- (net charge 0) + missing massive K_L0
alg_kl_kk.signal_side do |s|
  s.charged(kp: 1, km: 1)
  s.require_charge 0
  s.missing :K_L0      # massive missing particle, mass fixed to nominal K_L0
end

# Fit: 4-momentum constraint with the missing K_L0
alg_kl_kk.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_kl_kk.apply
alg_kl_kk.execute_on([data_3773, incMC_3773, exMC_d0_kl_kk])