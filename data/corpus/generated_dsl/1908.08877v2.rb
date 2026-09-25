# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # psi(3770) real data, 2.93 fb^-1 (rounds 03+04)
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # tag-based reformed inclusive MC at 3.773 GeV

# --- Decay cards (EvtGen syntax) -------------------------------------------
# The tag side (D-) is not decayed explicitly in these cards: it is left to the
# default EvtGen decay table so that all five hadronic tag modes are populated
# in the generated sample.
decay_card_tau = <<~DECAYCARD
    Decay psi(3770)
    1.00000 D+ D-
    Enddecay

    Decay D+
    1.00000 tau+ nu_tau        PHOTOS PHSP;
    Enddecay

    Decay tau+
    1.00000 pi+ anti-nu_tau    PHSP;
    Enddecay

    End
DECAYCARD

decay_card_mu = <<~DECAYCARD
    Decay psi(3770)
    1.00000 D+ D-
    Enddecay

    Decay D+
    1.00000 mu+ nu_mu          PHOTOS PHSP;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC samples (100k events each) -------------------------------
# Signal channel: D+ -> tau+ nu_tau, tau+ -> pi+ anti-nu_tau
exMC_tau = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_DtoTauNu"
  config.related_dataset = psi3770_data
  config.events          = 100_000
  config.decay_card      = decay_card_tau
  config.cross_section   = :default
end

# Normalization channel: D+ -> mu+ nu_mu
exMC_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_DtoMuNu"
  config.related_dataset = psi3770_data
  config.events          = 100_000
  config.decay_card      = decay_card_mu
  config.cross_section   = :default
end

### Event selection (BOSS) — single-tag D- analysis ###
# Signal channel: D+ -> tau+ nu_tau, tau+ -> pi+ anti-nu_tau
alg_name_tau = "STDmTauNu"
alg_tau = TagAnalysis.new(alg_name_tau)
alg_tau.set_header(["#{alg_name_tau}Alg/#{alg_name_tau}.h"])
       .set_constant({"ECMS" => [:double, 3.773]})

# Tag side: the D- is taken from the pre-stored DTag candidates in the five
# hadronic modes; charm = -1 pins the reconstructed D- (single tag).
alg_tau.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm -1
end

# Signal side: exactly one pi+ recoiling against the tag, plus the missing
# (tau-)neutrino system -> massless missing-particle form.
alg_tau.signal_side do |s|
  s.charged(pip: 1)     # charged multiset: one pi+ and nothing else
  s.require_charge 1    # total signal-side charge = +1
  s.missing :nu_tau     # massless missing particle (D+ -> tau+ nu_tau, tau+ -> pi+ anti-nu_tau)
end

# Four-momentum-constrained (4C) kinematic fit: tag + pi+ + missing neutrino -> CMS.
alg_tau.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200        # loose chi2 cut in BOSS; nominal fit
end

alg_tau.note(:background_veto, "additional vetoes applied in the ROOT analysis after the nominal 4C fit: " \
              "EMC energy over momentum E/p < 0.95 for pion-like tracks (E_EMC > 300 MeV); " \
              "maximum extra-shower energy < 300 MeV; |cos(theta_missing)| < 0.95 (mu-like) / 0.75 (pi-like); " \
              "opening angle between the missing momentum and the charged track > 25 deg (mu-like) / 45 deg (pi-like)")
       .note(:peaking_background, "peaking backgrounds from D+ -> pi0 pi+ and D+ -> K_L0 pi+ are estimated from " \
              "data control samples with a kernel method and subtracted from the mu-like (E_EMC <= 300 MeV) and " \
              "pi-like (E_EMC > 300 MeV) m_miss^2 distributions, which are fitted simultaneously")

alg_tau.with_decay_card(decay_card_tau).apply   # tag analysis: apply takes no Selection argument
alg_tau.execute_on([psi3770_data, psi3770_incMC, exMC_tau])

# Normalization channel: D+ -> mu+ nu_mu — identical tag scheme and selection
# chain, differing only in the signal-side charged species and decay card.
alg_name_mu = "STDmMuNu"
alg_mu = TagAnalysis.new(alg_name_mu)
alg_mu.set_header(["#{alg_name_mu}Alg/#{alg_name_mu}.h"])
      .set_constant({"ECMS" => [:double, 3.773]})

alg_mu.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi
  t.charm -1
end

alg_mu.signal_side do |s|
  s.charged(mup: 1)     # one mu+ on the signal side
  s.require_charge 1
  s.missing :nu_mu      # massless missing neutrino
end

alg_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200        # same nominal 4C fit as the tau channel
end

alg_mu.note(:background_veto, "same E/p, extra-shower energy, |cos(theta_missing)| and opening-angle vetoes as " \
             "the signal channel, applied in ROOT (pi-like/mu-like EMC separation at 300 MeV)")
      .note(:peaking_background, "peaking backgrounds from D+ -> pi0 pi+ and D+ -> K_L0 pi+ estimated from data " \
             "control samples with a kernel method")

alg_mu.with_decay_card(decay_card_mu).apply
alg_mu.execute_on([psi3770_data, psi3770_incMC, exMC_mu])