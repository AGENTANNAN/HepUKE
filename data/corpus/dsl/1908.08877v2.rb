### Dataset preparation ###
# psi(3770) data (2.93 fb^-1) and inclusive MC
psip3770_data = DatasetManager.real_data.find("712_3773")
psip3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for D+ -> tau+ nu_tau, tau+ -> pi+ anti-nu_tau (signal)
# At psi(3770), KKMC generates psi(3770) -> D+ D-
decay_card_Dp_taunu = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay
    Decay D+
    1.000 tau+ nu_tau PHSP;
    Enddecay
    Decay tau+
    1.000 pi+ anti-nu_tau PHSP;
    Enddecay
    End
DECAYCARD

# Decay card for D+ -> mu+ nu_mu (cross-check / normalization channel)
decay_card_Dp_munu = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay
    Decay D+
    1.000 mu+ nu_mu PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC samples
exMC_Dp_taunu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_to_tau_nu"
  config.related_dataset = psip3770_data
  config.events = 100000
  config.decay_card = decay_card_Dp_taunu
  config.cross_section = :default
end

exMC_Dp_munu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_to_mu_nu"
  config.related_dataset = psip3770_data
  config.events = 100000
  config.decay_card = decay_card_Dp_munu
  config.cross_section = :default
end

### D+ -> tau+ nu_tau (ST + missing: tag D- in hadronic modes, signal pi+ from tau->pi+ anti-nu_tau) ###

alg_Dp_taunu = TagAnalysis.new("DpToTauNu")
alg_Dp_taunu.set_header(["DpToTauNuAlg/DpToTauNu.h"])
             .set_constant({"ECMS" => [:double, 3.773]})
             .note(:signal_extraction, "MM2 (missing mass squared) method used for signal extraction; mu-like (E_EMC <= 300 MeV) and pi-like (E_EMC > 300 MeV) samples fitted simultaneously")
             .note(:background_veto, "Post-selection requirements not expressible in DSL: E_EMC/|pc| < 0.95 for pi-like sample; E_max < 300 MeV for extra showers; |cos(theta_missing)| < 0.95(mu)/0.75(pi); opening angle alpha > 25deg(mu)/45deg(pi); applied in BOSS code")
             .note(:peaking_backgrounds, "Peaking bkg from D+->pi0 pi+ and D+->K_L0 pi+ estimated from data control samples with kernel estimation method")

alg_Dp_taunu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm -1
end

alg_Dp_taunu.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge 1
  s.missing :nu_e
end

alg_Dp_taunu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_Dp_taunu.with_decay_card(decay_card_Dp_taunu)
alg_Dp_taunu.apply
alg_Dp_taunu.execute_on([psip3770_data, psip3770_incMC, exMC_Dp_taunu])

### D+ -> mu+ nu_mu (normalization channel, same tag scheme) ###

alg_Dp_munu = TagAnalysis.new("DpToMuNu")
alg_Dp_munu.set_header(["DpToMuNuAlg/DpToMuNu.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            .note(:tag_mode_unavailable, "DptoKKPi tag mode not available in authoritative mode list; omitted from tag side modes")

alg_Dp_munu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm -1
end

alg_Dp_munu.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge 1
  s.missing :nu_e
end

alg_Dp_munu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_Dp_munu.with_decay_card(decay_card_Dp_munu)
alg_Dp_munu.apply
alg_Dp_munu.execute_on([psip3770_data, psip3770_incMC, exMC_Dp_munu])