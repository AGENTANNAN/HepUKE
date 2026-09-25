### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")       # ψ(3770) real data at 3.773 GeV (20.3 fb^-1)
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")    # corresponding inclusive MC sample

# Decay card for the signal process (e+ channel): ψ(3770) -> D+ D-, with the tag
# D- -> K+ pi- pi- and the signal D+ -> K_S0 pi0 e+ nu_e (K_S0 -> pi+ pi-, pi0 -> gamma gamma)
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay D+
    1.0000 K_S0 pi0 e+ nu_e PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dp_KsPi0Enu"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS, tag-based) ###
alg_name = "DpToKsPi0EnuTag"
tag_alg = TagAnalysis.new(alg_name)
tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.773]})
       .with_decay_card(decay_card_signal)

# Tag side: D- reconstructed in six hadronic modes (declared as their D+ DTagAlg
# channel names and pinned to the D- side with charm -1)
tag_alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,      # D- -> K+ pi- pi-
          :DptoKsPi,       # D- -> K_S0 pi-
          :DptoKPiPiPi0,   # D- -> K+ pi- pi- pi0
          :DptoKsPiPi0,    # D- -> K_S0 pi- pi0
          :DptoKsPiPiPi,   # D- -> K_S0 pi+ pi- pi-
          :DptoKKPi        # D- -> K+ K- pi-
  t.charm -1
end

# Signal side: D+ -> K_S0 pi0 e+ nu_e (everything the tag did not use)
tag_alg.signal_side do |s|
  s.photons 2                                       # exactly two photons from pi0 -> gamma gamma
  s.charged(ep: 1, pip: 1, pim: 1, at_least: true)  # at least one e+ plus the K_S0 -> pi+ pi- pair
  s.missing :nu_e                                   # missing (massless) neutrino
end

# Four-momentum constraint kinematic fit with chi2 < 200
tag_alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

tag_alg.note(:separate_channel_spec, "The muon channel D+ -> K_S0 pi0 mu+ nu is a separate "
  "specification; its additional selection requires L(mu) > L(e), L(mu) > L(K), L(mu) > L(pi), "
  "EMC energy < 0.3 GeV and maximum extra photon energy < 0.2 GeV.")

tag_alg.apply
root_files = tag_alg.execute_on([data_3773, incMC_3773, exMC_signal])