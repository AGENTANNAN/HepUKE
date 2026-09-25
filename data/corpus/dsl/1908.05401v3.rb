### Dataset preparation ###
# psi(3770) data (2.93 fb^-1) and inclusive MC
psip3770_data = DatasetManager.real_data.find("712_3773")
psip3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for D+ -> phi X (representative channel: D+ -> phi pi+)
# At psi(3770), KKMC generates psi(3770) -> D+ D-
decay_card_Dp_phiX = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay
    Decay D+
    1.000 phi pi+ PHSP;
    Enddecay
    Decay phi
    1.000 K+ K- VSS;
    Enddecay
    End
DECAYCARD

# Decay card for D0 -> phi X (representative channel: D0 -> phi K- pi+)
decay_card_D0_phiX = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay
    Decay D0
    1.000 phi K- pi+ PHSP;
    Enddecay
    Decay phi
    1.000 K+ K- VSS;
    Enddecay
    End
DECAYCARD

# Exclusive MC samples
exMC_Dp_phiX = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_to_phi_X_inclusive"
  config.related_dataset = psip3770_data
  config.events = 100000
  config.decay_card = decay_card_Dp_phiX
  config.cross_section = :default
end

exMC_D0_phiX = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0_to_phi_X_inclusive"
  config.related_dataset = psip3770_data
  config.events = 100000
  config.decay_card = decay_card_D0_phiX
  config.cross_section = :default
end

### D+ -> phi X (tag D- in hadronic modes, signal phi -> K+ K- on D+ side) ###

alg_Dp_phiX = TagAnalysis.new("DpToPhiX")
alg_Dp_phiX.set_header(["DpToPhiXAlg/DpToPhiX.h"])
           .set_constant({"ECMS" => [:double, 3.773]})
           .note(:qc_correction, "QC correction factor f_QC not applied for charged D tags (taken as unity)")
           .note(:inclusive_signal, "Signal side is inclusive phi X; all known D+ -> phi X decays included in exclusive MC for efficiency determination")

alg_Dp_phiX.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm -1
end

alg_Dp_phiX.signal_side do |s|
  s.charged(kp: 1, km: 1)
  s.require_charge 0
end

alg_Dp_phiX.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_Dp_phiX.with_decay_card(decay_card_Dp_phiX)
alg_Dp_phiX.apply
alg_Dp_phiX.execute_on([psip3770_data, psip3770_incMC, exMC_Dp_phiX])

### D0 -> phi X (tag anti-D0 in hadronic modes, signal phi -> K+ K- on D0 side) ###

alg_D0_phiX = TagAnalysis.new("D0ToPhiX")
alg_D0_phiX.set_header(["D0ToPhiXAlg/D0ToPhiX.h"])
           .set_constant({"ECMS" => [:double, 3.773]})
           .note(:qc_correction, "Quantum correlation (QC) correction factor f_QC applied for neutral D tags following Refs. [23, 24]; determined from external inputs, not expressible in DSL")
           .note(:inclusive_signal, "Signal side is inclusive phi X; all known D0 -> phi X decays included in exclusive MC for efficiency determination")

alg_D0_phiX.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

alg_D0_phiX.signal_side do |s|
  s.charged(kp: 1, km: 1)
  s.require_charge 0
end

alg_D0_phiX.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_D0_phiX.with_decay_card(decay_card_D0_phiX)
alg_D0_phiX.apply
alg_D0_phiX.execute_on([psip3770_data, psip3770_incMC, exMC_D0_phiX])