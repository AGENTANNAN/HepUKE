### Dataset description ###
# psi(3770) real data (3.773 GeV, ~7.9 fb^-1) and the matching inclusive MC sample
data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for the signal process: psi(3770) -> D0 anti-D0,
#   D0 -> K0 pi- mu+ nu_mu  (K0 -> K_S0 -> pi+ pi-),  anti-D0 -> K+ pi-
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K0 pi- mu+ nu_mu PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay K0
    1.000 K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 2M exclusive signal MC events matched to the 3.773 GeV dataset
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_D0_K0pimumu"
  config.related_dataset = data_3773
  config.events         = 2_000_000
  config.decay_card     = decay_card_signal
  config.cross_section  = :default
end

### Event selection (BOSS) — D-tag single-tag analysis ###
# Tag-based analysis: the tag side is the pre-stored anti-D0 (charm -1);
# the signal side is everything the tag did not use (D0 -> K0 pi- mu+ nu_mu).
alg = TagAnalysis.new("D0TagKstMuNu")
alg.set_header(["D0TagKstMuNuAlg/D0TagKstMuNu.h"])
   .set_constant({"ECMS" => [:double, 3.773]})   # lab-frame c.m. energy for the 4C fit
   .with_decay_card(decay_card_signal)
# The paper's sixth ST entry duplicates the fifth final state and cannot be separated;
# only the five separable hadronic ST modes are used in the tag reconstruction.
alg.note(:tag_mode_ambiguity,
  "the paper's sixth ST mode shares the same final state as K+ pi- pi- pi+ pi0 and cannot be separated; only the five listed hadronic ST modes are used")

# Tag side: anti-D0 (charm -1) reconstructed from the five hadronic ST modes.
# mBC / deltaE are stored unconditionally (store-not-cut) — no window declared.
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPiPi, :D0toKPiPi0, :D0toKsPiPi, :D0toKPiPiPiPi0
  t.charm -1
end

# Signal side: D0 -> K0 pi- mu+ nu_mu, with K0 -> K_S0 -> pi+ pi-
#   charged content: pi+ (from K_S0) = 1, pi- (K* + K_S0) = 2, mu+ = 1 ; net charge 0
#   one massless nu_mu missing (semileptonic tag)
alg.signal_side do |s|
  s.charged(pip: 1, pim: 2, mup: 1)
  s.require_charge(0)
  s.missing :nu_mu
end

# 4C kinematic fit: tag + signal + missing neutrino four-momenta constrained to the
# lab-frame c.m. energy (measured CMS 4-vector); chi2 < 200
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# Render the tag analysis (no Selection argument) and run on data + inclusive MC + signal MC
alg.apply
root_files = alg.execute_on([data_3773, incMC_3773, exMC_signal])