# =============================================================================
# BOSS DSL — psi(3770) double-tag measurement of D+ -> eta l+ nu_l (l = e, mu)
# sqrt(s) = 3.773 GeV ; hadronic D- tag through six tag modes.
# Tag-based analysis -> TagAnalysis surface (no Selection object).
# =============================================================================

### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")      # psi(3770) real data (20.3 fb^-1)
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC sample

# --- shared decay-card fragments -------------------------------------------
# Six hadronic D- tag modes on the opposite side of the signal D+.
tag_block = <<~TAG
  Decay D-
  0.30 K+ pi- pi-        PHSP;
  0.15 K_S0 pi-          PHSP;
  0.20 K+ pi- pi- pi0    PHSP;
  0.12 K_S0 pi- pi0      PHSP;
  0.10 K_S0 pi+ pi- pi-  PHSP;
  0.13 K+ K- pi-         PHSP;
  Enddecay
TAG

# Neutral daughters used by both the signal (3pi mode) and the tag modes.
common_block = <<~COMMON
  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay
COMMON

# --- four signal decay cards (D+ signal, D- through the six tag modes) ------
decay_card_gg_e = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D-;
  Enddecay

  Decay D+
  1.0 eta e+ nu_e PHSP;
  Enddecay

  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay

  #{tag_block}
  #{common_block}
  End
DECAYCARD

decay_card_3pi_e = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D-;
  Enddecay

  Decay D+
  1.0 eta e+ nu_e PHSP;
  Enddecay

  Decay eta
  1.0 pi+ pi- pi0 PHSP;
  Enddecay

  #{tag_block}
  #{common_block}
  End
DECAYCARD

decay_card_gg_mu = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D-;
  Enddecay

  Decay D+
  1.0 eta mu+ nu_mu PHSP;
  Enddecay

  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay

  #{tag_block}
  #{common_block}
  End
DECAYCARD

decay_card_3pi_mu = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D-;
  Enddecay

  Decay D+
  1.0 eta mu+ nu_mu PHSP;
  Enddecay

  Decay eta
  1.0 pi+ pi- pi0 PHSP;
  Enddecay

  #{tag_block}
  #{common_block}
  End
DECAYCARD

# --- exclusive MC: 500k events for each of the four signal modes ------------
exMC_gg_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_etagg_enu"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_gg_e
  config.cross_section   = :default
end

exMC_3pi_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_eta3pi_enu"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_3pi_e
  config.cross_section   = :default
end

exMC_gg_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_etagg_munu"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_gg_mu
  config.cross_section   = :default
end

exMC_3pi_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_eta3pi_munu"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_3pi_mu
  config.cross_section   = :default
end

### Event selection (BOSS) — tag-based (TagAnalysis) ###

# ---------------------------------------------------------------------------
# Mode 1 : D+ -> eta(gamma gamma) e+ nu_e
# ---------------------------------------------------------------------------
alg_gg_e = TagAnalysis.new("DpEtaGGE")
alg_gg_e.set_header(["DpEtaGGEAlg/DpEtaGGE.h"])
        .set_constant({ "ECMS" => [:double, 3.773] })          # sqrt(s) = 3.773 GeV
        .set_alias({ "std::vector<double>" => "Vdouble" })
        .with_decay_card(decay_card_gg_e)

alg_gg_e.tag_side(:Dplus) do |t|
  # six hadronic D- tag modes
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1                                                   # tag the D-
  # tag mBC/deltaE are stored unconditionally (windowed in the ROOT M_BC fit)
end

alg_gg_e.signal_side do |s|
  s.photons 2                       # eta -> gamma gamma (no signal-photon energy threshold)
  s.charged(ep: 1)                  # exactly one e+; exact count vetoes extra good charged tracks
  s.missing :nu_e                   # missing (massless) neutrino
end

alg_gg_e.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)  # gamma gamma -> eta
  f.chi2_cut 200
end

alg_gg_e.note(:electron_pid, "e+ PID selection: L_e > 0.8*(L_e+L_pi+L_K), L_e > 0.001, E_EMC/p > 0.8")
alg_gg_e.note(:extra_photon_energy_veto, "additional (non-signal) photon showers required to have E < 0.25 GeV")
alg_gg_e.note(:efficiency_curve, "M(gamma gamma) acceptance window (0.50, 0.57) GeV/c^2 applied to the eta candidate")

alg_gg_e.apply
root_files_gg_e = alg_gg_e.execute_on([data_3773, incMC_3773, exMC_gg_e])

# ---------------------------------------------------------------------------
# Mode 2 : D+ -> eta(pi+ pi- pi0) e+ nu_e
# ---------------------------------------------------------------------------
alg_3pi_e = TagAnalysis.new("DpEta3PiE")
alg_3pi_e.set_header(["DpEta3PiEAlg/DpEta3PiE.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })
         .set_alias({ "std::vector<double>" => "Vdouble" })
         .with_decay_card(decay_card_3pi_e)

alg_3pi_e.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

alg_3pi_e.signal_side do |s|
  s.photons 2                       # two photons from pi0 -> gamma gamma
  s.charged(ep: 1, pip: 1, pim: 1)  # e+ plus the eta -> pi+ pi- pi0 pions; exact => extra tracks vetoed
  s.missing :nu_e
end

alg_3pi_e.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim, :gamma, :gamma).between(0.53, 0.57)   # M(pi+ pi- pi0) eta window
  f.chi2_cut 200
end

alg_3pi_e.note(:electron_pid, "e+ PID selection: L_e > 0.8*(L_e+L_pi+L_K), L_e > 0.001, E_EMC/p > 0.8")
alg_3pi_e.note(:extra_photon_energy_veto, "additional (non-signal) photon showers required to have E < 0.25 GeV")

alg_3pi_e.apply
root_files_3pi_e = alg_3pi_e.execute_on([data_3773, incMC_3773, exMC_3pi_e])

# ---------------------------------------------------------------------------
# Mode 3 : D+ -> eta(gamma gamma) mu+ nu_mu
# ---------------------------------------------------------------------------
alg_gg_mu = TagAnalysis.new("DpEtaGGMu")
alg_gg_mu.set_header(["DpEtaGGMuAlg/DpEtaGGMu.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })
         .set_alias({ "std::vector<double>" => "Vdouble" })
         .with_decay_card(decay_card_gg_mu)

alg_gg_mu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

alg_gg_mu.signal_side do |s|
  s.photons 2                       # eta -> gamma gamma
  s.charged(mup: 1)                 # exactly one mu+ (DSL fixed muon probability recipe applies)
  s.missing :nu_mu
end

alg_gg_mu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_gg_mu.note(:muon_pid, "mu+ PID selection: L_mu > L_K, L_mu > L_e, L_mu > 0.001, E_EMC in (0.101, 0.282) GeV; DSL default probability recipe requires prob_mu >= 0.001, prob_mu > prob_e, prob_mu > prob_K, prob_mu/(prob_mu+prob_pi+prob_K) >= 0.43")
alg_gg_mu.note(:extra_photon_energy_veto, "additional (non-signal) photon showers required to have E < 0.25 GeV")
alg_gg_mu.note(:extra_pi0_veto, "no extra pi0 allowed in the muon channel")
alg_gg_mu.note(:background_veto, "M(eta mu+) < 1.722 GeV/c^2 vetoed against the peaking D+ -> eta pi+ background")
alg_gg_mu.note(:efficiency_curve, "M(gamma gamma) acceptance window (0.50, 0.57) GeV/c^2 applied to the eta candidate")

alg_gg_mu.apply
root_files_gg_mu = alg_gg_mu.execute_on([data_3773, incMC_3773, exMC_gg_mu])

# ---------------------------------------------------------------------------
# Mode 4 : D+ -> eta(pi+ pi- pi0) mu+ nu_mu
# ---------------------------------------------------------------------------
alg_3pi_mu = TagAnalysis.new("DpEta3PiMu")
alg_3pi_mu.set_header(["DpEta3PiMuAlg/DpEta3PiMu.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
          .set_alias({ "std::vector<double>" => "Vdouble" })
          .with_decay_card(decay_card_3pi_mu)

alg_3pi_mu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

alg_3pi_mu.signal_side do |s|
  s.photons 2                       # two photons from pi0 -> gamma gamma
  s.charged(mup: 1, pip: 1, pim: 1) # mu+ plus the eta -> pi+ pi- pi0 pions; exact => extra tracks vetoed
  s.missing :nu_mu
end

alg_3pi_mu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim, :gamma, :gamma).between(0.53, 0.57)   # M(pi+ pi- pi0) eta window
  f.chi2_cut 200
end

alg_3pi_mu.note(:muon_pid, "mu+ PID selection: L_mu > L_K, L_mu > L_e, L_mu > 0.001, E_EMC in (0.101, 0.282) GeV; DSL default probability recipe requires prob_mu >= 0.001, prob_mu > prob_e, prob_mu > prob_K, prob_mu/(prob_mu+prob_pi+prob_K) >= 0.43")
alg_3pi_mu.note(:extra_photon_energy_veto, "additional (non-signal) photon showers required to have E < 0.25 GeV")
alg_3pi_mu.note(:extra_pi0_veto, "no extra pi0 allowed in the muon channel")
alg_3pi_mu.note(:background_veto, "M(eta mu+) < 1.710 GeV/c^2 vetoed against the peaking D+ -> eta pi+ background")

alg_3pi_mu.apply
root_files_3pi_mu = alg_3pi_mu.execute_on([data_3773, incMC_3773, exMC_3pi_mu])