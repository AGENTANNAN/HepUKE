# 2508.17819v2: Search for CP-forbidden D0 anti-D0 -> (K_S0 pi0)(K_S0 pi0) at psi(3770)
# Double-tag method, upper limit, 20.28 fb^-1 at 3.773 GeV
# Both D0 mesons reconstructed via DTagAlg double-tag (DT) pattern

DatasetManager.load_real_data("config/BES3_dataset.md")
DatasetManager.load_inclusive_mc("config/BES3_incMC.md")

data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

decay_card = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.000 K_S0 pi0 PHSP;
  Enddecay
  Decay anti-D0
  1.000 K_S0 pi0 PHSP;
  Enddecay
  End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_D0D0_KsPi0_KsPi0"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

alg = TagAnalysis.new("D0D0CpSearch")
alg.set_header(["D0D0CpSearchAlg/D0D0CpSearch.h"])
    .set_constant({ "ECMS" => [:double, 3.773] })
    .with_decay_card(decay_card)

alg.tag_side(:D0) do |t|
  t.modes :D0toKsPi0
  t.charm 1
end

alg.tag_side(:D0) do |t|
  t.modes :D0toKsPi0
  t.charm -1
  t.rank_by :inv
end

alg.signal_side do |s|
  s.photons 0
end

alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg.note(:signal_sideband_method, "Three regions in 2D M_BC plane: region I (both D0 in signal [1.859,1.873] GeV), region II (one signal one sideband [1.830,1.850] or [1.878,1.900] GeV), region III (both sideband). N_sig = N_I - N_II/2 + N_III/4. Applied in ROOT.")
alg.note(:upper_limit, "No significant signal observed (N_obs = -19 +/- 10). Upper limit at 90% C.L. set in ROOT analysis.")
alg.note(:cp_forbidden, "CP-forbidden process: both D0 and anti-D0 decay to the same CP eigenstate K_S0 pi0.")

alg.apply
alg.execute_on([data_3773, incMC_3773, sig_mc])