# BOSS Ruby DSL for arXiv:2007.07674v2
# "Observation of the Doubly Cabibbo-Suppressed Decay D+ -> K+ pi+ pi- pi0
#  and Evidence for D+ -> K+ omega"
#
# Analysis: e+e- -> psi(3770) -> D+ D- at sqrt(s) = 3.773 GeV
#   Double-tag method: tag D- via CF decays, reconstruct signal D+ from remaining tracks.
#   Signal channels: D+ -> K+ pi+ pi- pi0 (DCS, with pi0 -> gamma gamma)
#                    D+ -> K+ omega, omega -> pi+ pi- pi0 (sub-resonance)
#   D- tag modes: D- -> K+ pi- pi-, D- -> K_S0 pi-, D- -> K+ pi- pi- pi0
#   BF measured using DT yields vs ST yields. CP asymmetry also measured.
#   D+ -> K+ omega selected by additional ROOT-level cuts on M(3pi), cos(theta_omega), lambda/lambda_max.

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for signal D+ -> K+ pi+ pi- pi0 (includes K+ eta, K+ omega, K+ phi contributions)
# with D- -> CF tag modes
decay_card_k3pi_pi0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K+ pi+ pi- pi0 PHSP;
    Enddecay

    Decay D-
    0.340 K+ pi- pi- PHSP;
    0.270 K_S0 pi- PHSP;
    0.390 K+ pi- pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay
End
DECAYCARD

# Decay card for D+ -> K+ omega, omega -> pi+ pi- pi0
decay_card_komega = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K+ omega PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay D-
    0.340 K+ pi- pi- PHSP;
    0.270 K_S0 pi- PHSP;
    0.390 K+ pi- pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay
End
DECAYCARD

# Exclusive MC for D+ -> K+ pi+ pi- pi0 (inclusive of resonances)
exMC_k3pi_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "dplus_k3pi_pi0_signal_mc"
  config.related_dataset = psi3770_data
  config.events = 1_000_000
  config.decay_card = decay_card_k3pi_pi0
  config.cross_section = :default
end

# Exclusive MC for D+ -> K+ omega
exMC_komega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "dplus_komega_signal_mc"
  config.related_dataset = psi3770_data
  config.events = 500_000
  config.decay_card = decay_card_komega
  config.cross_section = :default
end

# ===========================================================================
# TagAnalysis: D- tag + signal D+ -> K+ pi+ pi- pi0
#   Both D+ -> K+ pi+ pi- pi0 (inclusive) and D+ -> K+ omega are selected here;
#   omega sub-sample separated in ROOT via M(pi+pi-pi0) and helicity cuts.
# ===========================================================================
alg = TagAnalysis.new("DplusToK3PiPi0")
alg.set_header(["DplusToK3PiPi0Alg/DplusToK3PiPi0.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(decay_card_k3pi_pi0)
   .note(:ks_veto, "Signal side: |M(pi+pi-) - M(K_S0)| > 20 MeV/c2 to reject D+ -> K_S0 K+ pi0 peaking background; applied in ROOT after reconstruction")
   .note(:opening_angle_cut, "Opening angle between D+ and D- candidates > 160 degrees to suppress non-DDbar events; applied in ROOT")
   .note(:tag_deltae_windows, "Tag deltaE window: (-25, 25) MeV for modes without pi0, (-55, 40) MeV for modes with pi0; applied in ROOT on mBC_tag vs mBC_sig 2D fit")
   .note(:omega_selection, "D+ -> K+ omega selected in ROOT by |M(pi+pi-pi0) - M_omega| < 40 MeV/c2, |cos(theta_omega)| > 0.57, lambda/lambda_max > 0.21")

# Tag side: D- reconstructed via three CF decay modes
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0
  t.charm -1
end

# Signal side: D+ -> K+ pi+ pi- pi0
# Content: 1 K+, 2 pi+, 1 pi-, pi0 -> gamma gamma (2 photons)
alg.signal_side do |s|
  s.photons 2
  s.charged(kp: 1, pip: 2, pim: 1)
  s.require_charge 1
end

# Kinematic fit: 4-momentum conservation + pi0 mass constraint
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg.apply
alg.execute_on([psi3770_data, psi3770_incMC, exMC_k3pi_pi0, exMC_komega])