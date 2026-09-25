# === Shared datasets ===
psip_data = DatasetManager.real_data.find("712_3773")
psip_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================================
# Mode I: D0 -> omega gamma_prime, omega -> pi+ pi- pi0, pi0 -> gamma gamma
# ============================================================================

decay_card_mode1 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 omega gamma_prime PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_D0_omega_gammap"
  config.related_dataset = psip_data
  config.events = 100_000
  config.decay_card = decay_card_mode1
  config.cross_section = :default
end

alg_mode1 = TagAnalysis.new("D0ToOmegaGammap")
alg_mode1.set_header(["D0ToOmegaGammapAlg/D0ToOmegaGammap.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })
         .with_decay_card(decay_card_mode1)

alg_mode1.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

alg_mode1.signal_side do |s|
  s.photons 2
  s.charged(pip: 1, pim: 1)
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :gamma_prime, mass: nil
end

alg_mode1.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)
  f.invariant_mass_of(:pip, :pim, :gamma, :gamma).between(0.700, 0.850)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:pip, :pim, :gamma, :gamma).constrain_to_nominal_mass_of(:omega)
  f.chi2_cut 200
end

alg_mode1
  .note(:background_veto, "E_oth_gamma_tot < 0.1 GeV veto on extra photon energy; |cos_theta_recoil| < 0.7 to suppress endcap background")
  .apply

alg_mode1.execute_on([psip_data, psip_incMC, exMC_mode1])

# ============================================================================
# Mode II: D0 -> gamma gamma_prime
# ============================================================================

decay_card_mode2 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 gamma gamma_prime PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_D0_gamma_gammap"
  config.related_dataset = psip_data
  config.events = 100_000
  config.decay_card = decay_card_mode2
  config.cross_section = :default
end

alg_mode2 = TagAnalysis.new("D0ToGammaGammap")
alg_mode2.set_header(["D0ToGammaGammapAlg/D0ToGammaGammap.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })
         .with_decay_card(decay_card_mode2)

alg_mode2.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

alg_mode2.signal_side do |s|
  s.photons 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.5
  s.missing :gamma_prime, mass: nil
end

alg_mode2.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mode2
  .note(:background_veto, "E_oth_gamma_tot < 0.1 GeV veto on extra photon energy; |cos_theta_recoil| < 0.7 to suppress endcap background")
  .apply

alg_mode2.execute_on([psip_data, psip_incMC, exMC_mode2])