### Dataset preparation ###
# Six BOSS 703 data points at sqrt(s) = 4.178 - 4.226 GeV
data_4180 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")
data_points = [data_4180, data_4190, data_4200, data_4210, data_4220, data_4230]

# Matching inclusive MC samples at the same six energy points
incmc_4180 = DatasetManager.inclusive_mc.find("703_4180")
incmc_4190 = DatasetManager.inclusive_mc.find("703_4190")
incmc_4200 = DatasetManager.inclusive_mc.find("703_4200")
incmc_4210 = DatasetManager.inclusive_mc.find("703_4210")
incmc_4220 = DatasetManager.inclusive_mc.find("703_4220")
incmc_4230 = DatasetManager.inclusive_mc.find("703_4230")
incmc_points = [incmc_4180, incmc_4190, incmc_4200, incmc_4210, incmc_4220, incmc_4230]

# Decay card for the D*+ -> e+ nu_e signal mode
# (top mother psi(4260) per KKMC convention; tag D*- -> anti-D0 pi-, anti-D0 -> K+ pi-)
decay_card_enu = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D*+ D*- PHSP;
    Enddecay

    Decay D*-
    1.0000 anti-D0 pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay D*+
    1.0000 e+ nu_e PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the D*+ -> mu+ nu_mu signal mode (identical tag chain)
decay_card_munu = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D*+ D*- PHSP;
    Enddecay

    Decay D*-
    1.0000 anti-D0 pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay D*+
    1.0000 mu+ nu_mu PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for each signal mode; one sample per energy point
exMC_enu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_dstartag_enu"
  config.events        = 100_000
  config.decay_card    = decay_card_enu
  config.cross_section = :default
end

exMC_munu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_dstartag_munu"
  config.events        = 100_000
  config.decay_card    = decay_card_munu
  config.cross_section = :default
end

### Event selection (BOSS) ###
# --- Channel I: single D*- tag + D*+ -> e+ nu_e (missing massless neutrino) ---
alg_e = TagAnalysis.new("DstarTagENu")
alg_e.set_header(["DstarTagENuAlg/DstarTagENu.h"])
     .set_constant({"ECMS" => [:double, 4.178]})
     .with_decay_card(decay_card_enu)

# Tag side: charm -1 D0 hadronic modes Kpi, Kpipi0, Kpipipi (D*+ tag modes not covered)
alg_e.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: exactly one charged positron of charge +1 and one missing massless neutrino
alg_e.signal_side do |s|
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

# 4C fit: tag D*-, signal lepton and neutrino constrained to the CMS four-momentum
alg_e.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_e.note(:photon_veto, "extra photons vetoed by requiring the maximum unused photon energy to be below 0.3 GeV; no DSL primitive exists for an energy ceiling on unused showers, so this cut must be enforced in the generated selection code")
     .note(:pid_correction_method, "signal positron identified with L'(e) > 0.001, L'(e)/(L'(e)+L'(pi)+L'(K)) > 0.8 and E/p > 0.8; the DSL ep key falls back to the default SimplePIDSvc recipe, so these custom electron-PID criteria must be re-implemented in the generated code")
     .note(:tag_mode_coverage, "only the three charm -1 D0 tag modes (Kpi, Kpipi0, Kpipipi) are included in the tag side; D+ tag modes proceeding through D- pi0 are deliberately not used")

alg_e.apply
alg_e.execute_on(data_points + incmc_points + exMC_enu)

# --- Channel II: single D*- tag + D*+ -> mu+ nu_mu (missing massless neutrino) ---
# Shares the identical tag-side selection chain with the electron channel.
alg_mu = TagAnalysis.new("DstarTagMuNu")
alg_mu.set_header(["DstarTagMuNuAlg/DstarTagMuNu.h"])
      .set_constant({"ECMS" => [:double, 4.178]})
      .with_decay_card(decay_card_munu)

alg_mu.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: exactly one charged muon of charge +1 and one missing massless neutrino
alg_mu.signal_side do |s|
  s.charged(mup: 1)
  s.require_charge 1
  s.missing :nu_mu
end

# Same 4C fit hypothesis: tag D*-, signal lepton and neutrino to the CMS four-momentum
alg_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mu.note(:photon_veto, "extra photons vetoed by requiring the maximum unused photon energy to be below 0.3 GeV; no DSL primitive exists for an energy ceiling on unused showers, so this cut must be enforced in the generated selection code")
      .note(:pid_correction_method, "signal muon identified with EMC energy in (0, 0.3) GeV together with MUC hit-depth requirements; the DSL mup key falls back to the fixed v1 probability recipe, so these custom muon-PID criteria must be re-implemented in the generated code")
      .note(:tag_mode_coverage, "only the three charm -1 D0 tag modes (Kpi, Kpipi0, Kpipipi) are included in the tag side; D+ tag modes proceeding through D- pi0 are deliberately not used")

alg_mu.apply
alg_mu.execute_on(data_points + incmc_points + exMC_munu)