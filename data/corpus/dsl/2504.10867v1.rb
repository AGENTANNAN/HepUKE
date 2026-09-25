# Paper: 2504.10867v1 — Precise measurement of the form factors in D0 -> K*(892)- mu+ nu_mu
# Tag-based analysis: D-tag ST method at psi(3770) (3.773 GeV, 7.9 fb^-1)
# Signal: D0 -> K0 pi- mu+ nu_mu (K0 -> K_S0 -> pi+ pi-)
# 6 ST modes for anti-D0 tag

### Dataset preparation ###
data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Exclusive MC: D0 -> K0 pi- mu+ nu_mu, anti-D0 -> inclusive hadronic tag
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

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_D0_K0pimumu_nu"
  config.related_dataset = data_3773
  config.events = 2_000_000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

### TagAnalysis (D-tag ST method) ###
alg = TagAnalysis.new("D0KstarMuNu")
alg.set_header(["D0KstarMuNuAlg/D0KstarMuNu.h"])
    .set_constant({"ECMS" => [:double, 3.773]})
    .with_decay_card(decay_card_signal)

# Single tag: anti-D0 reconstructed from 5 hadronic tag modes (6th mode unavailable)
# Modes: K+ pi-, K+ pi- pi- pi+, K+ pi- pi0, K_S0 pi+ pi-, K+ pi- pi- pi+ pi0
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPiPi, :D0toKPiPi0, :D0toKsPiPi, :D0toKPiPiPiPi0
  t.charm -1   # anti-D0 tag
end

# Signal side: K_S0 -> pi+ pi-, D0 -> pi- mu+ nu_mu
# Charged tracks: 1 pi+ (from K_S0), 2 pi- (from K_S0 + D0), 1 mu+ (from D0)
alg.signal_side do |s|
  s.charged(pip: 1, pim: 2, mup: 1)
  s.require_charge 0       # +1 -2 +1 = 0
  s.missing :nu_mu         # massless neutrino
end

# 4C kinematic fit: tag + signal + neutrino = ecms_lab
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# Paper Table I has 6 entries but only 5 distinct DTagAlg modes.
# Entries 5 (K+ pi- pi- pi+ pi0) and 6 (K+ pi- pi+ pi- pi0) both
# correspond to :D0toKPiPiPiPi0 — the same final state (K + 3pi + pi0).
# The paper applies different DeltaE windows to these sub-topologies,
# but DTagAlg treats them as one mode; the DSL cannot separate them.
alg.note(:tag_mode_limitation,
  "Paper Table I has 6 ST entries but only 5 distinct DTagAlg modes; entries 5-6 both map to :D0toKPiPiPiPi0. The paper's separate DeltaE windows for these sub-channels cannot be expressed at the DTagAlg level.")

# Post-BOSS cuts applied in ROOT:
# - M(K0 pi- mu+ (pi0)) < 1.60 GeV/c^2 (anti-hadronic background veto)
# - E_gamma_max < 0.15 GeV (extra pi0 suppression)
# - U_miss fit and window cut
alg.note(:background_veto,
  "M(K0 pi- mu+ (pi0)) < 1.60 GeV/c^2 veto against D0 -> K0 pi+ pi- (pi0) background; E_gamma_max < 0.15 GeV to suppress extra pi0 events; U_miss windows applied in ROOT")

alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_signal])