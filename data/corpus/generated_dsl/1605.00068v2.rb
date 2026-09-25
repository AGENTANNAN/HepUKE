### Dataset preparation ###
# ψ(3770) real data at 3.773 GeV (2.93 fb⁻¹) + matching inclusive MC
psipp_data  = DatasetManager.real_data.find("712_3773")
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card 1 — signal MC: D+ → K̄0 μ+ν with K̄0 → K_S0 → π+π−
decay_card_ks_pipi = <<~DECAYCARD
  Decay psi(3770)
  1.000 D+ D- PHSP;
  Enddecay

  Decay D-
  1.000 K+ pi- pi- PHSP;
  Enddecay

  Decay D+
  1.000 K_S0 mu+ nu_mu PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Decay card 2 — signal MC: D+ → K̄0 μ+ν with K̄0 → K_S0 → π0π0
decay_card_ks_pi0pi0 = <<~DECAYCARD
  Decay psi(3770)
  1.000 D+ D- PHSP;
  Enddecay

  Decay D-
  1.000 K+ pi- pi- PHSP;
  Enddecay

  Decay D+
  1.000 K_S0 mu+ nu_mu PHSP;
  Enddecay

  Decay K_S0
  1.000 pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Two 500k-event exclusive signal samples
exMC_ks_pipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_Dp_KS0munu_KS0topipim"
  config.related_dataset = psipp_data
  config.events          = 500_000
  config.decay_card      = decay_card_ks_pipi
  config.cross_section   = :default
end

exMC_ks_pi0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_Dp_KS0munu_KS0topi0pi0"
  config.related_dataset = psipp_data
  config.events          = 500_000
  config.decay_card      = decay_card_ks_pi0pi0
  config.cross_section   = :default
end

### Event selection (BOSS) — D− tag + D+ → K̄0 μ+ν double tag ###
# Tag side  : D− in six hadronic single-tag modes (pre-stored tag candidates)
# Signal side: the tracks/showers the tag did not use + one missing ν_μ

# ------------------------------------------------------------------
# Chain A — signal D+ → K̄0 μ+ν, K̄0 → K_S0 → π+π−
# ------------------------------------------------------------------
alg_kspipi = TagAnalysis.new("KsPiPiMuNuTag")
alg_kspipi.set_header(["KsPiPiMuNuTagAlg/KsPiPiMuNuTag.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          # tag-side charged-track quality (not exposed by the tag DSL)
          .note(:tag_track_selection,
                "tag-side charged tracks require |cosθ| < 0.93 and, unless originating from a K_S0, V_xy < 1.0 cm and |V_z| < 10 cm")
          # tag-side K_S0 reconstruction
          .note(:tag_ks0_reconstruction,
                "tag-side K_S0 formed from π+π− without PID; daughter |V_z| < 20 cm, |M(π+π−) − m_K_S0| < 12 MeV, decay length > 2σ")
          # tag-side K/π separation
          .note(:pid_correction_method,
                "tag-side K/π separation uses the combined dE/dx + TOF confidence levels: CL_K > CL_π for kaons, CL_π > CL_K for pions")
          # tag-side photon and π0 treatment
          .note(:tag_pi0_reconstruction,
                "tag-side photons require EMC shower time within 700 ns, E > 25 MeV (barrel) / 50 MeV (endcap) and > 10° to any charged track; π0 candidates required in (0.115, 0.150) GeV with a mass-constrained fit")
          # per-mode tag combination choice
          .note(:tag_best_candidate,
                "for each tag mode only the minimum-|ΔE| combination is retained; ΔE windows ±25 MeV (K+π−π−, K_S0π−, K_S0π+π−π−, K+K−π−) and −55/+40 MeV (K+π−π−π0, K_S0π−π0)")
          # signal-side muon identification
          .note(:signal_muon_pid,
                "signal-side muon identified by CL_μ > CL_K, CL_μ > CL_e, CL_μ > 0.001 and 0.1 < E_EMC < 0.3 GeV")
          # signal-side K_S0 selection
          .note(:signal_ks0_selection,
                "signal-side K_S0 → π+π− required to satisfy |M(π+π−) − m_K_S0| < 12 MeV")
          # competing-process veto (applied with the ROOT-level U_miss fit)
          .note(:background_veto,
                "D+ → K̄0 π+(π0) background rejected via M(K̄0 μ+) < 1.6 GeV and E_max^extra γ < 0.15 GeV; the U_miss fit is performed afterwards in ROOT")
          .with_decay_card(decay_card_ks_pipi)

# one tag_side call ⇒ single tag; the D− is tagged in six hadronic modes
alg_kspipi.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :deltaE, min: -0.055, max: 0.040   # union of the requested per-mode ΔE windows
  t.window :mBC,    min: 1.863,  max: 1.877   # single-tag signal region
end

# signal side: μ+ plus the two K_S0 → π+π− daughters, and the missing ν_μ
alg_kspipi.signal_side do |s|
  s.charged(mup: 1, pip: 1, pim: 1)
  s.missing :nu_mu
  s.require_charge 1
end

# 4C kinematic fit: tag (D−) + π+π− + μ+ + massless ν = √s
alg_kspipi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_kspipi.apply
alg_kspipi.execute_on([psipp_data, psipp_incMC, exMC_ks_pipi])

# ------------------------------------------------------------------
# Chain B — signal D+ → K̄0 μ+ν, K̄0 → K_S0 → π0π0
# ------------------------------------------------------------------
alg_pi0pi0 = TagAnalysis.new("Pi0Pi0MuNuTag")
alg_pi0pi0.set_header(["Pi0Pi0MuNuTagAlg/Pi0Pi0MuNuTag.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .note(:tag_track_selection,
                "tag-side charged tracks require |cosθ| < 0.93 and, unless originating from a K_S0, V_xy < 1.0 cm and |V_z| < 10 cm")
          .note(:tag_ks0_reconstruction,
                "tag-side K_S0 formed from π+π− without PID; daughter |V_z| < 20 cm, |M(π+π−) − m_K_S0| < 12 MeV, decay length > 2σ")
          .note(:pid_correction_method,
                "tag-side K/π separation uses the combined dE/dx + TOF confidence levels: CL_K > CL_π for kaons, CL_π > CL_K for pions")
          .note(:tag_pi0_reconstruction,
                "tag-side photons require EMC shower time within 700 ns, E > 25 MeV (barrel) / 50 MeV (endcap) and > 10° to any charged track; π0 candidates required in (0.115, 0.150) GeV with a mass-constrained fit")
          .note(:tag_best_candidate,
                "for each tag mode only the minimum-|ΔE| combination is retained; ΔE windows ±25 MeV (K+π−π−, K_S0π−, K_S0π+π−π−, K+K−π−) and −55/+40 MeV (K+π−π−π0, K_S0π−π0)")
          .note(:signal_muon_pid,
                "signal-side muon identified by CL_μ > CL_K, CL_μ > CL_e, CL_μ > 0.001 and 0.1 < E_EMC < 0.3 GeV")
          # π0π0 reconstruction of the signal K̄0
          .note(:pi0pi0_selection,
                "the four signal photons must give M(π0π0) in (0.45, 0.51) GeV; the combination with the smallest summed π0 mass-constrained χ² is kept")
          .note(:background_veto,
                "D+ → K̄0 π+(π0) background rejected via M(K̄0 μ+) < 1.6 GeV and E_max^extra γ < 0.15 GeV; the U_miss fit is performed afterwards in ROOT")
          .with_decay_card(decay_card_ks_pi0pi0)

alg_pi0pi0.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :deltaE, min: -0.055, max: 0.040
  t.window :mBC,    min: 1.863,  max: 1.877
end

# signal side: μ+ plus the four photons from K_S0 → π0π0 → 4γ, and the missing ν_μ
alg_pi0pi0.signal_side do |s|
  s.photons 4
  s.charged(mup: 1)
  s.missing :nu_mu
  s.require_charge 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# 4C kinematic fit + the two π0 mass constraints
alg_pi0pi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_pi0pi0.apply
alg_pi0pi0.execute_on([psipp_data, psipp_incMC, exMC_ks_pi0pi0])