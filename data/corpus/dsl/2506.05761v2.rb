### TagAnalysis — semileptonic double-tag at psi(3770) ###
# Paper: Observation of D+ -> K_S0 pi0 mu+ nu_mu, test of LFU,
#   first angular analysis of D+ -> K*0 l+ nu_l
# Uses 20.3 fb^-1 at 3.773 GeV

psip3770_data = DatasetManager.real_data.find("712_3773")
psip3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

decay_card = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- VSS;
    Enddecay
    Decay D+
    1.0000 K_S0 pi0 mu+ nu_mu PHSP;
    Enddecay
    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay
    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
End
DECAYCARD

exMC_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_KsPi0munu_exclusive_mc"
  config.related_dataset = psip3770_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

alg = TagAnalysis.new("DpKsPi0LNu")
alg.set_header(["DpKsPi0LNuAlg/DpKsPi0LNu.h"])
  .set_constant({"ECMS" => [:double, 3.773]})

# Tag side: D- reconstructed via 6 hadronic decay modes (ST)
# Paper modes: K+pi-pi-, K_S0pi-, K+pi-pi-pi0, K_S0pi-pi0, K_S0pi+pi-pi-, K+K-pi-
# Mapped to Dplus mode symbols (DTagAlg handles CP conjugation)
alg.tag_side(:Dm) do
  modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0,
        :DptoKsPiPiPi, :DptoKKPi
end

# Signal side: D+ -> K_S0 pi0 l+ nu_l
# e+ channel: D+ -> K_S0 pi0 e+ nu_e
# mu+ channel: D+ -> K_S0 pi0 mu+ nu_mu (first observation)
alg.signal_side do
  charged(pip: 0, pim: 2, ep: 1, at_least: true)
  photons 2
  missing :nu_e
end

alg.fit do
  constrain_four_momentum
  chi2_cut 200
end

alg.with_decay_card(decay_card).apply
  .note(:two_lepton_channels, "both e+ and mu+ channels analysed separately; this spec covers the e+ channel; mu+ channel uses mu+ PID with EMC energy < 0.3 GeV")
  .note(:angular_analysis, "first full angular analysis of D+ -> K*0 l+ nu_l with CP asymmetries A_i and averaged observables S_i (i=2-9) measured in 5 q^2 bins; form factor ratios r_V=1.42, r_2=0.75 extracted from helicity amplitude fit")
  .note(:s_wave, "S-wave (K_S0 pi0)_S-wave component included via LASS parametrization; fraction ~6-7% of total decay rate")
  .note(:signal_mu_channel, "mu+ channel has additional PID requirement: L(mu) > L(e), L(mu) > L(K), L(mu) > L(pi); E_EMC^mu < 0.3 GeV; max extra photon energy < 0.2 GeV")
  .execute_on([psip3770_data, psip3770_incMC, exMC_mu])