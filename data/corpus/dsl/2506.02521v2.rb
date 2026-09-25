# =====================================================================
# BESIII: Measurement of BFs of D+ → η e+ ν_e and D+ → η μ+ ν_μ
#   [arXiv:2506.02521v2]
#
# ψ(3770) → D Dbar at 3.773 GeV, 20.3 fb⁻¹. Double-tag method.
# ST: 6 hadronic D- tag modes.
# DT signal: D+ → η ℓ+ ν_ℓ (ℓ = e, μ), η → γγ, π+π-π0.
# TagAnalysis (D-tag, double-tag semileptonic).
# =====================================================================

psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

sub_pi0_gg = <<~DECAYCARD
  Decay pi0
  1.0  gamma gamma  PHSP;
  Enddecay
DECAYCARD

sub_eta_gg = <<~DECAYCARD
  Decay eta
  1.0  gamma gamma  PHSP;
  Enddecay
DECAYCARD

sub_eta_3pi = <<~DECAYCARD
  Decay eta
  1.0  pi+ pi- pi0  PHSP;
  Enddecay
DECAYCARD

sub_Ks_pipi = <<~DECAYCARD
  Decay K_S0
  1.0  pi+ pi-  PHSP;
  Enddecay
DECAYCARD

# =====================================================================
# Algorithm 1: D+ → η(→γγ) e+ ν_e
# =====================================================================

decay_card_e_gg = <<~DECAYCARD
  Decay psi(3770)
  1.0  D+ D-  PHSP;
  Enddecay

  Decay D-
  0.1667  K+ pi- pi-          PHSP;
  0.1667  K_S0 pi-            PHSP;
  0.1667  K+ pi- pi- pi0      PHSP;
  0.1667  K_S0 pi- pi0        PHSP;
  0.1667  K_S0 pi+ pi- pi-    PHSP;
  0.1667  K+ K- pi-           PHSP;
  Enddecay

  Decay D+
  1.0  eta e+ nu_e            PHSP;
  Enddecay

  #{sub_eta_gg}
  #{sub_Ks_pipi}
  #{sub_pi0_gg}
  End
DECAYCARD

exMC_e_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dp_to_eta_gg_e_nu"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_e_gg
  config.cross_section   = :default
end

alg_e_gg = TagAnalysis.new("DpTagEtaGGElectronNu")
alg_e_gg.set_header(["DpTagEtaGGElectronNuAlg/DpTagEtaGGElectronNu.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .with_decay_card(decay_card_e_gg)

alg_e_gg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,
          :DptoKsPi,
          :DptoKPiPiPi0,
          :DptoKsPiPi0,
          :DptoKsPiPiPi,
          :DptoKKPi
  t.charm -1
end

alg_e_gg.signal_side do |s|
  s.photons 2                                    # 2 photons from η→γγ
  s.charged(ep: 1)                               # positron (electron PID)
  s.missing :nu_e
end

alg_e_gg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_e_gg
  .note(:eta_selection,
        "η → γγ: M(γγ) ∈ (0.50, 0.57) GeV/c². 1C kinematic fit to constrain γγ mass " \
        "to nominal η mass, χ² < 50. Best candidate by minimum χ².")
  .note(:electron_pid,
        "Positron: L_e > 0.8 × (L_e + L_π + L_K), L_e > 0.001. " \
        "E_EMC / p_track > 0.8 for D+ → η e+ ν_e.")
  .note(:background_veto,
        "No additional good charged tracks (N_extra^char = 0) on signal side. " \
        "Max extra photon energy E_extra^max(γ) < 0.25 GeV. " \
        "No extra π0 allowed for D+ → η μ+ ν_μ.")
  .note(:signal_extraction,
        "Signal yield from unbinned ML fit to U_miss ≡ E_miss - |p_miss|c. " \
        "Signal: MC shape ⊗ Gaussian. Combinatorial background: inclusive MC shape. " \
        "U_miss peaks near 0 for signal. p_D+ = -p_D- × sqrt(E_beam²/c² - m_D-²c²).")
  .note(:st_yields,
        "ST yields from unbinned ML fit to M_BC. MC signal shape ⊗ double-Gaussian. " \
        "Background: ARGUS function. M_BC ∈ (1.863, 1.877) GeV/c² for D- tag window. " \
        "ΔE requirements: mode-dependent, from Table 2 of the paper.")

alg_e_gg.apply
alg_e_gg.execute_on([psi3770_data, psi3770_incMC, exMC_e_gg])

# =====================================================================
# Algorithm 2: D+ → η(→π+π-π0) e+ ν_e
# =====================================================================

decay_card_e_3pi = <<~DECAYCARD
  Decay psi(3770)
  1.0  D+ D-  PHSP;
  Enddecay

  Decay D-
  0.1667  K+ pi- pi-          PHSP;
  0.1667  K_S0 pi-            PHSP;
  0.1667  K+ pi- pi- pi0      PHSP;
  0.1667  K_S0 pi- pi0        PHSP;
  0.1667  K_S0 pi+ pi- pi-    PHSP;
  0.1667  K+ K- pi-           PHSP;
  Enddecay

  Decay D+
  1.0  eta e+ nu_e            PHSP;
  Enddecay

  #{sub_eta_3pi}
  #{sub_Ks_pipi}
  #{sub_pi0_gg}
  End
DECAYCARD

exMC_e_3pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dp_to_eta_3pi_e_nu"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_e_3pi
  config.cross_section   = :default
end

alg_e_3pi = TagAnalysis.new("DpTagEta3PiElectronNu")
alg_e_3pi.set_header(["DpTagEta3PiElectronNuAlg/DpTagEta3PiElectronNu.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_e_3pi)

alg_e_3pi.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,
          :DptoKsPi,
          :DptoKPiPiPi0,
          :DptoKsPiPi0,
          :DptoKsPiPiPi,
          :DptoKKPi
  t.charm -1
end

alg_e_3pi.signal_side do |s|
  s.photons 2                                    # 2 photons from π0
  s.charged(pip: 1, pim: 1, ep: 1)               # π+π- from η, plus e+
  s.missing :nu_e
end

alg_e_3pi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_e_3pi
  .note(:eta_3pi_selection,
        "η → π+π-π0: M(π+π-π0) ∈ (0.53, 0.57) GeV/c². " \
        "η invariant mass constrained via 1C fit. " \
        "The η mass constraint on the π+π-π0 system is not directly expressible in the " \
        "Tag DSL's fit block (which operates on participant-level symbols only).")
  .note(:electron_pid,
        "Positron PID: L_e > 0.8 × (L_e + L_π + L_K), L_e > 0.001. E/p > 0.8.")
  .note(:background_veto,
        "No extra charged tracks. E_extra^max(γ) < 0.25 GeV. " \
        "No extra π0 for μ+ channel. " \
        "For D+ → η μ+ ν_μ: M(ημ+) < 1.722 GeV/c² (γγ mode) or < 1.710 GeV/c² (3π mode) " \
        "to suppress D+ → η π+ peaking background.")

alg_e_3pi.apply
alg_e_3pi.execute_on([psi3770_data, psi3770_incMC, exMC_e_3pi])

# =====================================================================
# Algorithm 3: D+ → η(→γγ) μ+ ν_μ
# =====================================================================

decay_card_mu_gg = <<~DECAYCARD
  Decay psi(3770)
  1.0  D+ D-  PHSP;
  Enddecay

  Decay D-
  0.1667  K+ pi- pi-          PHSP;
  0.1667  K_S0 pi-            PHSP;
  0.1667  K+ pi- pi- pi0      PHSP;
  0.1667  K_S0 pi- pi0        PHSP;
  0.1667  K_S0 pi+ pi- pi-    PHSP;
  0.1667  K+ K- pi-           PHSP;
  Enddecay

  Decay D+
  1.0  eta mu+ nu_mu          PHSP;
  Enddecay

  #{sub_eta_gg}
  #{sub_Ks_pipi}
  #{sub_pi0_gg}
  End
DECAYCARD

exMC_mu_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dp_to_eta_gg_mu_nu"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_mu_gg
  config.cross_section   = :default
end

alg_mu_gg = TagAnalysis.new("DpTagEtaGGMuonNu")
alg_mu_gg.set_header(["DpTagEtaGGMuonNuAlg/DpTagEtaGGMuonNu.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_mu_gg)

alg_mu_gg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,
          :DptoKsPi,
          :DptoKPiPiPi0,
          :DptoKsPiPi0,
          :DptoKsPiPiPi,
          :DptoKKPi
  t.charm -1
end

alg_mu_gg.signal_side do |s|
  s.photons 2                                    # 2 photons from η→γγ
  s.charged(mup: 1)                              # muon (μ+)
  s.missing :nu_mu
end

alg_mu_gg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_mu_gg
  .note(:muon_pid,
        "Muon: L_μ > L_K, L_μ > L_e, L_μ > 0.001. " \
        "E_EMC ∈ (0.101, 0.282) GeV (optimized by S/√(S+B)). " \
        "Fixed v1 muon PID thresholds in generated code: probMuon ≥ 0.001, " \
        "probMuon > probElectron, probMuon > probKaon, " \
        "probMuon/(probMuon+probPion+probKaon) ≥ 0.43.")
  .note(:eta_selection,
        "η → γγ: M(γγ) ∈ (0.50, 0.57) GeV/c². 1C fit to η mass, χ² < 50.")
  .note(:background,
        "D+ → η π+ π0 peaking background in muon channel: yields fixed to 389 (γγ) " \
        "and 151 (3π) from mis-ID rate and PDG BFs in U_miss fits.")

alg_mu_gg.apply
alg_mu_gg.execute_on([psi3770_data, psi3770_incMC, exMC_mu_gg])

# =====================================================================
# Algorithm 4: D+ → η(→π+π-π0) μ+ ν_μ
# =====================================================================

decay_card_mu_3pi = <<~DECAYCARD
  Decay psi(3770)
  1.0  D+ D-  PHSP;
  Enddecay

  Decay D-
  0.1667  K+ pi- pi-          PHSP;
  0.1667  K_S0 pi-            PHSP;
  0.1667  K+ pi- pi- pi0      PHSP;
  0.1667  K_S0 pi- pi0        PHSP;
  0.1667  K_S0 pi+ pi- pi-    PHSP;
  0.1667  K+ K- pi-           PHSP;
  Enddecay

  Decay D+
  1.0  eta mu+ nu_mu          PHSP;
  Enddecay

  #{sub_eta_3pi}
  #{sub_Ks_pipi}
  #{sub_pi0_gg}
  End
DECAYCARD

exMC_mu_3pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Dp_to_eta_3pi_mu_nu"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_mu_3pi
  config.cross_section   = :default
end

alg_mu_3pi = TagAnalysis.new("DpTagEta3PiMuonNu")
alg_mu_3pi.set_header(["DpTagEta3PiMuonNuAlg/DpTagEta3PiMuonNu.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .with_decay_card(decay_card_mu_3pi)

alg_mu_3pi.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,
          :DptoKsPi,
          :DptoKPiPiPi0,
          :DptoKsPiPi0,
          :DptoKsPiPiPi,
          :DptoKKPi
  t.charm -1
end

alg_mu_3pi.signal_side do |s|
  s.photons 2                                    # 2 photons from π0
  s.charged(pip: 1, pim: 1, mup: 1)              # π+π- from η, plus μ+
  s.missing :nu_mu
end

alg_mu_3pi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mu_3pi
  .note(:eta_3pi_selection,
        "η → π+π-π0: M(π+π-π0) ∈ (0.53, 0.57) GeV/c². 1C fit constraining mass to η. " \
        "Sub-composition mass constraint (η on π+π-π0) not directly expressible in Tag fit block.")
  .note(:muon_pid,
        "Muon PID thresholds as described above. M(ημ+) < 1.722 (γγ) / < 1.710 (3π) GeV/c².")
  .note(:peaking_background,
        "D+ → η π+ π0 peaking background in muon channel fixed in U_miss fits.")

alg_mu_3pi.apply
alg_mu_3pi.execute_on([psi3770_data, psi3770_incMC, exMC_mu_3pi])