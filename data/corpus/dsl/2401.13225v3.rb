# =============================================================================
# Paper: 2401.13225v3
# Title: Observation of D+ → f0(500) μ+ ν_μ and study of
#        D+ → π+π- ℓ+ ν_ℓ decay dynamics
# Data:  2.93 fb^-1 e+e- collisions at √s = 3.773 GeV (ψ(3770))
# Method: TagAnalysis — single-tag D- mesons reconstructed in 6 hadronic
#         final states; semileptonic D+ → π+π- ℓ+ ν_ℓ on the signal side
# =============================================================================

# ---------------------------------------------------------------------------
# Dataset loads — real data, inclusive MC, and two exclusive MC samples
# ---------------------------------------------------------------------------
data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ---------------------------------------------------------------------------
# Decay cards for exclusive MC generation
# ---------------------------------------------------------------------------

# Muon channel: D+ → π+ π- μ+ ν_μ
decay_card_DpPiPiMuNu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-                            PHSP;
  Enddecay

  Decay D+
  1.0000 pi+ pi- mu+ nu_mu                PHOTOS  ISGW2;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi-                       PHSP;
  Enddecay

  End
DECAYCARD

# Electron channel: D+ → π+ π- e+ ν_e
decay_card_DpPiPiENu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-                            PHSP;
  Enddecay

  Decay D+
  1.0000 pi+ pi- e+ nu_e                  PHOTOS  ISGW2;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi-                       PHSP;
  Enddecay

  End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive MC samples
# ---------------------------------------------------------------------------
exMC_DpPiPiMuNu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Dp_pipi_mu_nu"
  c.related_dataset = data_3773
  c.events          = 500000
  c.decay_card      = decay_card_DpPiPiMuNu
  c.cross_section   = :default
end

exMC_DpPiPiENu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Dp_pipi_e_nu"
  c.related_dataset = data_3773
  c.events          = 500000
  c.decay_card      = decay_card_DpPiPiENu
  c.cross_section   = :default
end

# =============================================================================
# Algorithm 1 — Muon channel:  D+ → π+ π- μ+ ν_μ
# =============================================================================
alg_mu = TagAnalysis.new("DpToPiPiMuNu")
alg_mu.set_header(["DpToPiPiMuNuAlg/DpToPiPiMuNu.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .with_decay_card(decay_card_DpPiPiMuNu)

# Tag side: reconstruct D- in 6 hadronic modes (charm = -1 for D+ species = D-)
alg_mu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

# Signal side: D+ → π+ π- μ+ ν_μ
alg_mu.signal_side do |s|
  s.charged(pip: 1, pim: 1, mup: 1)
  s.require_charge 0
  s.missing :nu_mu                 # massless neutrino
end

# 4-momentum constrained kinematic fit
alg_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# ---------------------------------------------------------------------------
# Selection notes (paper-level cuts; most applied in ROOT analysis)
# ---------------------------------------------------------------------------
alg_mu.note(:U_miss_window,
            "U_miss selection window: events retained in (-0.025, 0.025) GeV; " \
            "applied at ROOT analysis level, not in BOSS selection.")
       .note(:E_extra_max_cut,
             "Maximum energy of unused (extra) photons < 0.25 GeV (E_extra_max cut); " \
             "suppresses background with additional neutral energy.")
       .note(:N_extra_char_zero,
             "No extra charged tracks allowed beyond the three signal-side tracks: " \
             "N_extra_char = 0 (applied in ROOT).")
       .note(:muon_PID,
             "Muon PID thresholds (applied in ROOT): CL_μ > 0.001, CL_μ > CL_e, " \
             "CL_μ > CL_π, CL_μ > CL_K; E_EMC in (0.09, 0.31) GeV " \
             "(muon MIP energy deposition window in EMC).")
       .note(:resonance_analysis,
             "Resonance content (f0(500) and ρ0) extracted from fits to the " \
             "M(π+π-) invariant mass spectrum in bins of q^2; performed at ROOT level.")
       .note(:q2_binning,
             "Form-factor analysis uses q^2 binning of the π+π- invariant mass " \
             "for simultaneous extraction of f0(500) and ρ0 contributions.")
       .note(:U_miss_improvement,
             "U_miss resolution improved using tag-side D- momentum (constrained " \
             "from beam energy and D- mass); applied at ROOT level.")

alg_mu.apply
alg_mu.execute_on([data_3773, incMC_3773, exMC_DpPiPiMuNu])

# =============================================================================
# Algorithm 2 — Electron channel:  D+ → π+ π- e+ ν_e
# =============================================================================
alg_e = TagAnalysis.new("DpToPiPiENu")
alg_e.set_header(["DpToPiPiENuAlg/DpToPiPiENu.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .with_decay_card(decay_card_DpPiPiENu)

# Tag side: same 6 D- hadronic modes as muon channel
alg_e.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

# Signal side: D+ → π+ π- e+ ν_e
alg_e.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)
  s.require_charge 0
  s.missing :nu_e
end

# 4-momentum constrained kinematic fit
alg_e.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# ---------------------------------------------------------------------------
# Selection notes (paper-level cuts; most applied in ROOT analysis)
# ---------------------------------------------------------------------------
alg_e.note(:U_miss_window,
           "U_miss selection window: events retained in (-0.025, 0.025) GeV; " \
           "applied at ROOT analysis level, not in BOSS selection.")
      .note(:E_extra_max_cut,
            "Maximum energy of unused (extra) photons < 0.25 GeV (E_extra_max cut); " \
            "suppresses background with additional neutral energy.")
      .note(:N_extra_char_zero,
            "No extra charged tracks allowed beyond the three signal-side tracks: " \
            "N_extra_char = 0 (applied in ROOT).")
      .note(:resonance_analysis,
            "Resonance content (f0(500) and ρ0) extracted from fits to the " \
            "M(π+π-) invariant mass spectrum in bins of q^2; performed at ROOT level.")
      .note(:q2_binning,
            "Form-factor analysis uses q^2 binning of the π+π- invariant mass " \
            "for simultaneous extraction of f0(500) and ρ0 contributions.")
      .note(:U_miss_improvement,
            "U_miss resolution improved using tag-side D- momentum (constrained " \
            "from beam energy and D- mass); applied at ROOT level.")
      .note(:combined_fit,
            "Combined fit of electron and muon channels for the f0(500) form-factor " \
            "extraction (common signal model, channel-specific efficiencies); " \
            "performed at ROOT level.")

alg_e.apply
alg_e.execute_on([data_3773, incMC_3773, exMC_DpPiPiENu])