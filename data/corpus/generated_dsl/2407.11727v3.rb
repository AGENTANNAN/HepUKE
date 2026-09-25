# =============================================================================
#  BOSS DSL — Leptonic decays D_s^+ -> mu^+ nu_mu and D_s^+ -> tau^+ nu_tau
#  (tau^+ -> e^+ nu nubar, mu^+ nu nubar, pi^+ nubar, rho^+ nubar with
#   rho^+ -> pi^+ pi^0) measured with the double-tag method in
#  e^+e^- -> D_s^{*+} D_s^{*-} in the psi(4260) energy region.
#  Tag side : D_s^{*-} -> gamma D_s^- (DTagTool tag candidates)
#  Signal   : the D_s^+ recoiling against the tag
# =============================================================================

### ---------------------------------------------------------------------------
### Dataset preparation
### ---------------------------------------------------------------------------

# 17 real-data samples spanning 4.237-4.699 GeV (10.64 fb^-1 in total)
data_points = [
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("703_4310"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("703_4390"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("703_4470"),
  DatasetManager.real_data.find("703_4530"),
  DatasetManager.real_data.find("703_4575"),
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("706_4610"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660")
]

# Matching inclusive MC at every energy point
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4237"),
  DatasetManager.inclusive_mc.find("703_4246"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4270"),
  DatasetManager.inclusive_mc.find("703_4280"),
  DatasetManager.inclusive_mc.find("703_4310"),
  DatasetManager.inclusive_mc.find("703_4360"),
  DatasetManager.inclusive_mc.find("703_4390"),
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("703_4470"),
  DatasetManager.inclusive_mc.find("703_4530"),
  DatasetManager.inclusive_mc.find("703_4575"),
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660")
]

# ---- Signal decay cards (one per signal mode) -------------------------------
# The tag D_s^- -> K^+ K^- pi^- is included as the reference tag mode;
# D_s^{*+-} -> gamma D_s^{+-} accounts for the D_s^* -> gamma D_s transition.

# Mode 1: D_s^+ -> mu^+ nu_mu
decay_card_munu = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s*-
  1.000 gamma D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 mu+ nu_mu PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode 2: D_s^+ -> tau^+ nu_tau, tau^+ -> e^+ nu_e anti-nu_tau
decay_card_taunu_e = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s*-
  1.000 gamma D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 tau+ nu_tau PHSP;
  Enddecay

  Decay tau+
  1.000 e+ nu_e anti-nu_tau PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode 3: D_s^+ -> tau^+ nu_tau, tau^+ -> mu^+ nu_mu anti-nu_tau
decay_card_taunu_mu = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s*-
  1.000 gamma D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 tau+ nu_tau PHSP;
  Enddecay

  Decay tau+
  1.000 mu+ nu_mu anti-nu_tau PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode 4: D_s^+ -> tau^+ nu_tau, tau^+ -> pi^+ anti-nu_tau
decay_card_taunu_pi = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s*-
  1.000 gamma D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 tau+ nu_tau PHSP;
  Enddecay

  Decay tau+
  1.000 pi+ anti-nu_tau PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode 5: D_s^+ -> tau^+ nu_tau, tau^+ -> rho^+ anti-nu_tau, rho^+ -> pi^+ pi^0
decay_card_taunu_rho = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s*- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s*-
  1.000 gamma D_s- PHSP;
  Enddecay

  Decay D_s+
  1.000 tau+ nu_tau PHSP;
  Enddecay

  Decay tau+
  1.000 rho+ anti-nu_tau PHSP;
  Enddecay

  Decay rho+
  1.000 pi+ pi0 VSS;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  End
DECAYCARD

# ---- Exclusive MC: 1M events at every energy point for each signal mode -----
exMC_munu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ds_munu"
  config.events        = 1_000_000
  config.decay_card    = decay_card_munu
  config.cross_section = :default
end

exMC_taunu_e = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ds_taunu_e"
  config.events        = 1_000_000
  config.decay_card    = decay_card_taunu_e
  config.cross_section = :default
end

exMC_taunu_mu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ds_taunu_mu"
  config.events        = 1_000_000
  config.decay_card    = decay_card_taunu_mu
  config.cross_section = :default
end

exMC_taunu_pi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ds_taunu_pi"
  config.events        = 1_000_000
  config.decay_card    = decay_card_taunu_pi
  config.cross_section = :default
end

exMC_taunu_rho = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ds_taunu_rho"
  config.events        = 1_000_000
  config.decay_card    = decay_card_taunu_rho
  config.cross_section = :default
end

# All samples the algorithms are run on
all_datasets = data_points + incMC_points +
               exMC_munu + exMC_taunu_e + exMC_taunu_mu + exMC_taunu_pi + exMC_taunu_rho

### ---------------------------------------------------------------------------
### Event selection (BOSS) — tag-based (TagAnalysis)
### ---------------------------------------------------------------------------

# 13 D_s^- tag modes
tag_modes = [
  :DstoKKPi,             # K^+ K^- pi^-
  :DstoKKPiPi0,          # K^+ K^- pi^- pi^0
  :DstoPiPiPi,           # pi^+ pi^- pi^-
  :DstoKsK,              # K_S^0 K^-
  :DstoKsKPi0,           # K_S^0 K^- pi^0
  :DstoKPiPi,            # K^+ pi^- pi^-
  :DstoKsKsPi,           # K_S^0 K_S^0 pi^-
  :DstoKsKPiPi,          # K_S^0 K^+ pi^- pi^-
  :DstoKsPiPiK,          # K_S^0 K^- pi^+ pi^-
  :DstoPiEta,            # pi eta
  :DstoEtaPrimePiPiPi,   # eta' pi pi pi
  :DstoPiEtaPrime,       # pi eta' (eta' -> rho^0 gamma)
  :DstoPiPi0Eta          # pi pi^0 eta
]

# ---------------------------------------------------------------------------
# Mode 1 — D_s^+ -> mu^+ nu_mu
# ---------------------------------------------------------------------------
alg_munu = TagAnalysis.new("DsToMuNu")
alg_munu.set_header(["DsToMuNuAlg/DsToMuNu.h"])
        .set_constant({ "ECMS" => [:double, 4.260] })
        .with_decay_card(decay_card_munu)

alg_munu.tag_side(:Ds) do |t|
  t.modes(*tag_modes)          # 13 D_s^- tag modes
  t.charm -1                   # tag the D_s^- side
  t.window :deltaE, abs: 0.05  # tag-side DeltaE window (best |DeltaE| candidate kept)
end

alg_munu.signal_side do |s|
  s.charged(mup: 1)            # exactly one extra charged track: mu^+
  s.require_charge(+1)         # net charge of the signal side
  s.missing :nu_mu             # the muon neutrino is undetected
end

alg_munu.fit do |f|
  f.constrain_four_momentum    # 4-momentum conservation
  f.chi2_cut 200               # chi^2 < 200
end

alg_munu
  .note(:background_veto, "K_S^0 peaking veto |M(pi+pi-) - m(K_S^0)| < 0.03 GeV/c^2
        applied for the pi+pi-pi- and K+pi+pi- tag modes to reject
        K_S^0 -> pi+pi- contamination of the tag side")
  .note(:pid_correction_method, "muon identification uses the MUC hit-depth requirement
        (mu_a) with a mu_b cross-check; the efficiency is not a flat per-track number and
        is taken from control samples")
  .note(:efficiency_curve, "per-mode DeltaE windows keep only the tag candidate with the
        smallest |DeltaE| in each mode and charge; the single-tag yields are extracted from
        a fit to the M_BC distribution, which is stored unconditionally")
  .apply

# ---------------------------------------------------------------------------
# Mode 2 — D_s^+ -> tau^+ nu_tau, tau^+ -> e^+ nu_e anti-nu_tau
# ---------------------------------------------------------------------------
alg_taunu_e = TagAnalysis.new("DsToTauNuE")
alg_taunu_e.set_header(["DsToTauNuEAlg/DsToTauNuE.h"])
           .set_constant({ "ECMS" => [:double, 4.260] })
           .with_decay_card(decay_card_taunu_e)

alg_taunu_e.tag_side(:Ds) do |t|
  t.modes(*tag_modes)
  t.charm -1
  t.window :deltaE, abs: 0.05
end

alg_taunu_e.signal_side do |s|
  s.charged(ep: 1)             # exactly one extra charged track: e^+
  s.require_charge(+1)
  s.missing :nu                # undetected neutrinos (nu_tau + nu_e + anti-nu_tau)
end

alg_taunu_e.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_taunu_e
  .note(:background_veto, "K_S^0 peaking veto |M(pi+pi-) - m(K_S^0)| < 0.03 GeV/c^2
        applied for the pi+pi-pi- and K+pi+pi- tag modes")
  .note(:efficiency_curve, "per-mode DeltaE windows keep only the tag candidate with the
        smallest |DeltaE| in each mode and charge; single-tag yields from a fit to M_BC")
  .apply

# ---------------------------------------------------------------------------
# Mode 3 — D_s^+ -> tau^+ nu_tau, tau^+ -> mu^+ nu_mu anti-nu_tau
# ---------------------------------------------------------------------------
alg_taunu_mu = TagAnalysis.new("DsToTauNuMu")
alg_taunu_mu.set_header(["DsToTauNuMuAlg/DsToTauNuMu.h"])
            .set_constant({ "ECMS" => [:double, 4.260] })
            .with_decay_card(decay_card_taunu_mu)

alg_taunu_mu.tag_side(:Ds) do |t|
  t.modes(*tag_modes)
  t.charm -1
  t.window :deltaE, abs: 0.05
end

alg_taunu_mu.signal_side do |s|
  s.charged(mup: 1)            # exactly one extra charged track: mu^+
  s.require_charge(+1)
  s.missing :nu                # undetected neutrinos (nu_tau + nu_mu + anti-nu_tau)
end

alg_taunu_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_taunu_mu
  .note(:background_veto, "K_S^0 peaking veto |M(pi+pi-) - m(K_S^0)| < 0.03 GeV/c^2
        applied for the pi+pi-pi- and K+pi+pi- tag modes")
  .note(:pid_correction_method, "muon identification uses the MUC hit-depth requirement
        (mu_a) with a mu_b cross-check; efficiency taken from control samples")
  .note(:efficiency_curve, "per-mode DeltaE windows keep only the tag candidate with the
        smallest |DeltaE| in each mode and charge; single-tag yields from a fit to M_BC")
  .apply

# ---------------------------------------------------------------------------
# Mode 4 — D_s^+ -> tau^+ nu_tau, tau^+ -> pi^+ anti-nu_tau
# ---------------------------------------------------------------------------
alg_taunu_pi = TagAnalysis.new("DsToTauNuPi")
alg_taunu_pi.set_header(["DsToTauNuPiAlg/DsToTauNuPi.h"])
            .set_constant({ "ECMS" => [:double, 4.260] })
            .with_decay_card(decay_card_taunu_pi)

alg_taunu_pi.tag_side(:Ds) do |t|
  t.modes(*tag_modes)
  t.charm -1
  t.window :deltaE, abs: 0.05
end

alg_taunu_pi.signal_side do |s|
  s.charged(pip: 1)            # exactly one extra charged track: pi^+
  s.require_charge(+1)
  s.missing :nu                # undetected neutrino (anti-nu_tau + nu_tau)
end

alg_taunu_pi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_taunu_pi
  .note(:background_veto, "K_S^0 peaking veto |M(pi+pi-) - m(K_S^0)| < 0.03 GeV/c^2
        applied for the pi+pi-pi- and K+pi+pi- tag modes")
  .note(:efficiency_curve, "per-mode DeltaE windows keep only the tag candidate with the
        smallest |DeltaE| in each mode and charge; single-tag yields from a fit to M_BC")
  .apply

# ---------------------------------------------------------------------------
# Mode 5 — D_s^+ -> tau^+ nu_tau, tau^+ -> rho^+ anti-nu_tau, rho^+ -> pi^+ pi^0
# ---------------------------------------------------------------------------
alg_taunu_rho = TagAnalysis.new("DsToTauNuRho")
alg_taunu_rho.set_header(["DsToTauNuRhoAlg/DsToTauNuRho.h"])
             .set_constant({ "ECMS" => [:double, 4.260] })
             .with_decay_card(decay_card_taunu_rho)

alg_taunu_rho.tag_side(:Ds) do |t|
  t.modes(*tag_modes)
  t.charm -1
  t.window :deltaE, abs: 0.05
end

alg_taunu_rho.signal_side do |s|
  s.charged(pip: 1)            # the pi^+ from rho^+ -> pi^+ pi^0
  s.photons 2                  # the two photons from pi^0 -> gamma gamma
  s.min_photon_angle 10.0      # minimum photon angle to charged tracks (degrees)
  s.require_charge(+1)
  s.missing :nu                # undetected neutrino (anti-nu_tau + nu_tau)
end

alg_taunu_rho.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # M(gamma gamma) -> m(pi^0)
  f.chi2_cut 200
end

alg_taunu_rho
  .note(:background_veto, "K_S^0 peaking veto |M(pi+pi-) - m(K_S^0)| < 0.03 GeV/c^2
        applied for the pi+pi-pi- and K+pi+pi- tag modes")
  .note(:efficiency_curve, "per-mode DeltaE windows keep only the tag candidate with the
        smallest |DeltaE| in each mode and charge; single-tag yields from a fit to M_BC")
  .note(:simultaneous_fit, "the four tau sub-modes (tau -> e nu nu, mu nu nu, pi nu, rho nu)
        are fitted simultaneously to extract the tau^+ nu_tau branching fraction; the
        simultaneous fit itself is performed in the ROOT stage")
  .apply

### ---------------------------------------------------------------------------
### Execution
### ---------------------------------------------------------------------------
root_files_munu      = alg_munu.execute_on(all_datasets)
root_files_taunu_e   = alg_taunu_e.execute_on(all_datasets)
root_files_taunu_mu  = alg_taunu_mu.execute_on(all_datasets)
root_files_taunu_pi  = alg_taunu_pi.execute_on(all_datasets)
root_files_taunu_rho = alg_taunu_rho.execute_on(all_datasets)