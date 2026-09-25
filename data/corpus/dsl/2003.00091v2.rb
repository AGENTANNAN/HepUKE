# DSL for 2003.00091v2: Quantum-correlated D meson pairs at psi(3770)
# Strong-phase parameters c_i, s_i, c_i', s_i' from D->K_S/L0 pi+ pi-
# DT method with flavor, CP, mixed-CP, and K_L0 tags

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for psi(3770) -> D0 D0bar (KKMC convention)
decay_card_psipp = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay
  End
DECAYCARD

exMC_psipp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_DDbar_3773"
  config.related_dataset = psi3770_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_psipp
  config.cross_section   = :default
end

# TagAnalysis for DT: tag side 1 = hadronic/CP tags, tag side 2 = K_S0 pi+ pi- signal
alg = TagAnalysis.new("D0StrongPhaseDT")
alg.set_header(["D0StrongPhaseDTAlg/D0StrongPhaseDT.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })
    .with_decay_card(decay_card_psipp)

# Tag side 1: flavor tags, CP-even, CP-odd, mixed-CP tags
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi,
          :D0toKK, :D0toPiPi, :D0toKsPi0, :D0toPiPiPi0,
          :D0toKsEta, :D0toKsOmega,
          :D0toKsPiPi
  # Note: semileptonic tag D0->KeNu, K_L0 tags (K_L0 pi0, K_L0 pi0 pi0),
  # and some CP-odd modes (K_S0 eta'(pi+pi-eta), K_S0 eta'(gamma pi+pi-))
  # are not available in DTagAlg. Their ST yields are estimated from known BFs.
end

# Tag side 2: signal D0 -> K_S0 pi+ pi-
alg.tag_side(:D0) do |t|
  t.modes :D0toKsPiPi
  t.rank_by :inv
end

alg.signal_side do |s|
  s.photons 0
  # No leftover tracks: both D mesons are fully reconstructed in DT
end

alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg.note(:tag_mode_unavailable, "Semileptonic tag D0->KeNu, K_L0 pi0, K_L0 pi0 pi0 not in DTagAlg; ST yields from known BFs and N_DDbar.")
    .note(:partial_reconstruction, "K_S0 pi+ pi_miss and K_S0(pi0 pi0_miss) pi+ pi- partial reconstruction for signal side not expressible in TagAnalysis DSL.")
    .note(:kl_reconstruction, "K_L0 pi+ pi- signal reconstructed via missing-mass-squared technique; K_L0 not directly detectable.")
    .note(:dalitz_binning, "Dalitz plot binning (equal, optimal, modified-optimal) and efficiency-matrix correction applied at ROOT level.")
    .note(:peaking_background, "Peaking backgrounds from DCS decays and K_S0->pi0pi0 estimated from MC and subtracted; ROOT-level fit.")

alg.apply
alg.execute_on([psi3770_data, psi3770_incMC, exMC_psipp])