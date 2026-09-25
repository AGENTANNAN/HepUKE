### Dataset description ###
# psi(3770) data: 20.3 fb^-1 at sqrt(s) = 3.773 GeV
# Double-tag (DT) method used to measure D0/D+ -> KKKpi absolute BFs
psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for psi(3770) -> D0 D0bar (KKMC + EvtGen)
decay_card_D0D0bar = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for psi(3770) -> D+ D-
decay_card_DpDm = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K- pi+ pi+ PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_D0D0bar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exclusive_psi3770_D0D0bar"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_D0D0bar
  config.cross_section = :default
end

exMC_DpDm = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exclusive_psi3770_DpDm"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_DpDm
  config.cross_section = :default
end

### Tag-based event selection (BOSS) ###

# ============================================================
# D0 DT analysis at psi(3770)
# Tag side 1: anti-D0 tagged via D0toKPi, D0toKPiPi0, D0toKPiPiPi (charm -1)
# Tag side 2: D0 tagged via same modes (charm +1)
# Signal D0 decays: KS0 K+ K- pi0, KS0 KS0 K- pi+, KS0 KS0 K+ pi-, K+ K- K- pi+
#   plus phi-submodes: phi KS0 pi0, phi K- pi+
# Signal D0 reconstructed from remaining tracks not used by tags;
# identified via DeltaE_sig and M_BC^sig in ROOT analysis.
# ============================================================

alg_D0 = TagAnalysis.new("D0toKKKPi_DT")
alg_D0.set_header(["D0toKKKPi_DTAlg/D0toKKKPi_DT.h"])
        .set_constant({"ECMS" => [:double, 3.773]})

alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1   # tag anti-D0
end

alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1    # tag D0 (signal side)
  t.rank_by :inv
end

alg_D0.signal_side do |s|
  # No extra particles beyond the D Dbar pair at psi(3770)
  # Signal D0 decay products are the remaining tracks after both tags
end

alg_D0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_D0
  .note(:dt_method, "Double-tag method at psi(3770): ST anti-D0 tags (K+pi-, K+pi-pi0, K+pi-pi-pi+) select events; signal D0 decays identified among remaining tracks via DeltaE_sig and M_BC^sig; 5 inclusive KKKpi signal modes + 3 phi-dominated sub-modes measured. Signal D0 decay modes: KS0 K+ K- pi0, KS0 KS0 K- pi+, KS0 KS0 K+ pi-, K+ K- K- pi+, and phi sub-modes phi KS0 pi0, phi K- pi+.")
  .note(:st_yield, "total ST D0 yield = (1688.7 +/- 0.5) x 10^4 from M_BC fits; DT yields extracted from fits to M_BC^sig distributions")
  .note(:detection_efficiency, "efficiencies determined from exclusive MC with psi(3770) -> D0 D0bar; detection efficiencies ~1-5% depending on signal mode; efficiency averaged over tag modes weighted by ST yields")
  .with_decay_card(decay_card_D0D0bar)
  .apply

# ============================================================
# D+ DT analysis at psi(3770)
# Tag side 1: D- tagged via 6 hadronic modes (charm -1)
# Tag side 2: D+ tagged via same modes (charm +1)
# Signal D+ decays: KS0 K+ K- pi+, phi KS0 pi+
# ============================================================

alg_Dp = TagAnalysis.new("DptoKKKPi_DT")
alg_Dp.set_header(["DptoKKKPi_DTAlg/DptoKKKPi_DT.h"])
        .set_constant({"ECMS" => [:double, 3.773]})

alg_Dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1   # tag D-
end

alg_Dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm 1    # tag D+ (signal side)
  t.rank_by :inv
end

alg_Dp.signal_side do |s|
  # No extra particles beyond D+ D- pair at psi(3770)
end

alg_Dp.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Dplus)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:Dplus)
  f.chi2_cut 200
end

alg_Dp
  .note(:dt_method, "Double-tag method at psi(3770): ST D- tags (K+pi-pi-, KS0 pi-, K+pi-pi-pi0, KS0 pi-pi0, KS0 pi-pi-pi+, K+K-pi-) select events; signal D+ decays identified among remaining tracks. Signal modes: KS0 K+ K- pi+ and phi KS0 pi+ (phi -> K+ K-).")
  .note(:st_yield, "total ST D- yield = (1096.0 +/- 0.4) x 10^4 from M_BC fits")
  .with_decay_card(decay_card_DpDm)
  .apply

# Execute the algorithms
root_files_D0 = alg_D0.execute_on([psi3770_data, psi3770_incMC, exMC_D0D0bar])
root_files_Dp = alg_Dp.execute_on([psi3770_data, psi3770_incMC, exMC_DpDm])