# DSL for arxiv:1811.11349v2
# D0 -> Kbar0 pi- e+ nu_e at sqrt(s) = 3.773 GeV (BOSS 712, 3773)
# Tag-based: ST D0bar in 3 modes, DT signal: K_S0 + pi- + e+ + nu_e

decay_card_d0_semilep = <<~DECAYCARD
  Decay psi(3770)
  1 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1 anti-K0 pi- e+ nu_e PHSP;
  Enddecay
  Decay anti-D0
  1 K+ pi- PHSP;
  Enddecay
  End
DECAYCARD

data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

sig_d0_semilep = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_d0_kspienu"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_d0_semilep
  config.cross_section   = :default
end

alg = TagAnalysis.new("D0TagKsPiENu")
alg.set_header(["D0TagKsPiENuAlg/D0TagKsPiENu.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })
  .with_decay_card(decay_card_d0_semilep)

# ST tag side: anti-D0 reconstructed in 3 hadronic modes
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1   # anti-D0 (charm +1 = D0bar)
end

# Signal side: D0 -> K_S0 + pi- + e+ + nu_e
# K_S0 -> pi+ pi- (vertex fit handled externally, noted below)
# Total charged: 1 pi+ (from K_S0) + 2 pi- (one from K_S0, one from D0) + 1 e+
alg.signal_side do |s|
  s.charged(pip: 1, pim: 2, ep: 1)
  s.require_charge 0        # +1 + 2*(-1) + (+1) = 0
  s.missing :nu_e            # massless neutrino
end

# Kinematic fit: 4C + D mass constraint
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.note(:ks_reconstruction, "K_S0 -> pi+ pi- with secondary vertex fit. M(pi+pi-) within ±12 MeV of nominal K_S0 mass. Decay length significance > 2 sigma.")
alg.note(:ks_mass_constraint, "K_S0 mass constraint via secondary vertex fit of pi+ pi- pair from signal-side tracks — not expressible in tag_fit.invariant_mass_of")
alg.note(:pion_pid, "Pion PID: CL_pi > CL_K; electron PID: E/p > 0.8 and EMC shower shape")
alg.note(:electron_pid, "e+ candidate: E/p ratio and EMC energy deposit consistent with electron hypothesis")
alg.note(:umiss_fit, "Signal yield from fit to U_miss = E_miss - |p_miss| distribution; form factor parameters from 5D fit — ROOT-level analysis")
alg.note(:tag_deltae_windows, "Tag-mode-dependent DeltaE requirements applied in ST selection (see Table I in paper)")
alg.note(:ks_veto, "For D0 -> K_S0 pi- e+ nu_e signal, veto D0 -> K_S0 pi+ pi- background")

alg.apply
alg.execute_on([data_3773, incMC_3773, sig_d0_semilep])