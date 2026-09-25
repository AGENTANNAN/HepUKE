# =====================================================================
# Dataset preparation
# =====================================================================
data_3773  = DatasetManager.real_data.find("712_3773")       # ψ(3770) real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")    # Corresponding inclusive MC

# Signal decay card: e+e- -> D+ D-, signal D+ -> K+ K_S0 pi0 (EvtGen format)
decay_card_signal = <<~DECAYCARD
  Decay psi(3770)
  1.000 D+ D- PHSP;
  Enddecay

  Decay D+
  1.000 K+ K_S0 pi0 PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive signal MC: 2,000,000 events with the same decay card
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpKsKPi0"
  config.related_dataset = data_3773           # matched to the ψ(3770) data sample
  config.events          = 2_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# =====================================================================
# Tag analysis (tag-based: D-tag)
# =====================================================================
alg = TagAnalysis.new("DpKsKPi0Tag")
alg.set_header(["DpKsKPi0TagAlg/DpKsKPi0Tag.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })   # √s = 3.773 GeV
   .with_decay_card(decay_card_signal)

# --- Tag side 1: D- through six hadronic modes ---
alg.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi,        # D- -> K+ pi- pi-
          :DptoKPiPiPi0,     # D- -> K+ pi- pi- pi0
          :DptoKsPi,         # D- -> K_S0 pi-
          :DptoKsPiPi0,      # D- -> K_S0 pi- pi0
          :DptoKsPiPiPi,     # D- -> K_S0 pi- pi- pi+
          :DptoKKPi          # D- -> K+ K- pi-
  t.charm(-1)                # pin the tagged side to D-
end

# --- Tag side 2: D+ -> K- pi+ pi+, double tag ranked by invariant mass ---
alg.tag_side(:Dp) do |t|
  t.modes :DptoKPiPi         # D+ -> K- pi+ pi+
  t.charm(1)                 # pin the tagged side to D+
  t.rank_by(:inv)            # rank the double-tag candidate by invariant mass
end

# --- Signal side: exactly one charged K+, exactly two photons ---
alg.signal_side do |s|
  s.charged(kp: 1)           # exactly one charged K+
  s.photons 2                # exactly two photons (pi0 -> gamma gamma)
  s.min_photon_angle 10.0    # photon separation angle > 10 degrees
  s.require_charge 1         # net charge +1
end

# --- Kinematic fit ---
alg.fit do |f|
  f.constrain_four_momentum                                                 # 4-momentum conservation
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)    # gamma-gamma mass -> m(pi0)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Dplus)           # nominal D mass on tag candidate 1
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:Dplus)           # nominal D mass on tag candidate 2
  f.chi2_cut 200                                                            # accept chi2 < 200
end

alg.apply   # no Selection argument for tag-based analyses
alg.execute_on([data_3773, incMC_3773, exMC_signal])