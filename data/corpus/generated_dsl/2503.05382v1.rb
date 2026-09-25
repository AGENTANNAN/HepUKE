### Dataset description ###
# ψ(3770) at 3.773 GeV — 20.3 fb^-1 real data + matching inclusive MC
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the double-tag signal: ψ(3770) → D+ D-, D+ → K+ K- π+ π+ π- (phase space)
decay_card_Dp_KKpipipi = <<~DECAYCARD
    Decay psi(3770)
    1.0 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0 K+ K- pi+ pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the φ channel: D+ → φ π+ π+ π-, φ → K+ K- (same final state as above)
decay_card_Dp_phi = <<~DECAYCARD
    Decay psi(3770)
    1.0 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0 phi pi+ pi+ pi- PHSP;
    Enddecay

    Decay phi
    1.0 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Exclusive MC (1,000,000 events each) for the two commissioned signal modes
exMC_Dp_KKpipipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dp_KKpipipi"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_Dp_KKpipipi
  config.cross_section   = :default
end

exMC_Dp_phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dp_phi3pi"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_Dp_phi
  config.cross_section   = :default
end

### Event selection — tag-based double tag ###
alg = TagAnalysis.new("DpKKpipipiDTag")
alg.set_header(["DpKKpipipiDTagAlg/DpKKpipipiDTag.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })
   .with_decay_card(decay_card_Dp_KKpipipi)
   # ΔE windows are mode-dependent (Table 2); only the nominal ±60 MeV window is applied at BOSS level,
   # the final per-mode windows and the M_BC signal region are handled in ROOT.
   .note(:efficiency_curve, "tag-side ΔE windows are mode-dependent (Table 2, nominally ±60 MeV); a single nominal ±60 MeV window is applied at BOSS level and the final per-mode windows are applied in ROOT")
   # Modes not commissioned because their tag modes are absent from the mode vocabulary.
   .note(:uncommissioned_modes, "D+ -> K_S0 K+ pi+ pi- pi0, D+ -> K_S0 K+ eta and D+ -> K_S0 K+ omega are not implemented: the required tag modes are absent from the available DTagAlg tag-mode vocabulary, so only D+ -> K+ K- pi+ pi+ pi- and D+ -> phi pi+ pi+ pi- are commissioned")

# Tag side: D- reconstructed through the six hadronic modes (charm -1 pins the D- side)
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :deltaE, abs: 0.060
end

# Signal side: D+ → K+ K- π+ π+ π- built from the tracks left by the tag;
# pure double tag, no additional (photon / missing) signal-side particles
alg.signal_side do |s|
  s.charged(kp: 1, km: 1, pip: 2, pim: 1)
  s.require_charge 1
end

# 4C kinematic fit: constrain the total four-momentum, χ² < 200
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_Dp_KKpipipi, exMC_Dp_phi])