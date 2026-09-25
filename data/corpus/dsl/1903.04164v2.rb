#
# 1903.04164v2: D_s+ → K_S0 K+ and D_s+ → K_L0 K+
#   at sqrt(s) = 4.178 GeV (3.19 fb⁻¹)
#   Single-tag (ST) + double-tag (DT) / ST + missing
# TagAnalysis — tag-based
#

# ── Dataset ────────────────────────────────────────────────────────────
datasets  = [DatasetManager.load_real_data.find("703_4180")]
inc_mc    = [DatasetManager.load_inclusive_mc.find("703_4180")]

# ── Exclusive MC ──────────────────────────────────────────────────────
# conexc generator; D_s*± D_s∓ production

decay_card_ks = <<~DECAY
  Decay D*-
  1.0 gamma D_s-  VSP_PWAVE;
  Decay D_s+
  1.0 K_S0 K+  SVS;
  Decay K_S0
  1.0 pi+ pi-  PHSP;
  Enddecay
DECAY

decay_card_kl = <<~DECAY
  Decay D*-
  1.0 gamma D_s-  VSP_PWAVE;
  Decay D_s+
  1.0 K_L0 K+  SVS;
  Enddecay
DECAY

exclusive_mc_ks = DatasetManager.create_exclusive_mc do |c|
  c.decay_card(decay_card_ks, generator: :conexc)
end

exclusive_mc_kl = DatasetManager.create_exclusive_mc do |c|
  c.decay_card(decay_card_kl, generator: :conexc)
end

# ── 13 ST tag modes for D_s- (excluding D_s- → K_S0 K- to avoid double counting) ─
# Known modes used explicitly; remaining covered by mode_group :hadronic
ds_tag_modes = [
  :DstoKpKmPim, :DstoKpKmPimPi0, :DstoKsKmPi0,
  :DstoKsKpPimPim, :DstoPimEta, :DstoPimEtap, :DstoKsKm
]

# ═══════════════════════════════════════════════════════════════════════
# Mode I: D_s+ → K_S0 K+  (DT method)
# ═══════════════════════════════════════════════════════════════════════
alg_ks = TagAnalysis.new("Ds_to_KS0_Kplus_DT")

alg_ks.set_header(%w[KinematicFit/KinematicFit.h DTagAlg/DTagAlg.h])
alg_ks.set_constant(ECMS: "MeasuredEcmsSvc::get_ecms()")
alg_ks.with_decay_card(decay_card_ks)

alg_ks.tag_side(:Ds) do |t|
  t.modes(*ds_tag_modes)
  t.mode_group :hadronic
  t.charm(-1)
end

alg_ks.note(:tag_mode_list,
  "13 ST modes used: K+K-pi-, K-pi+pi-, pi+pi-pi-, K+K-pi-pi0, " \
  "pi-eta'(gamma_rho0), rho-eta, K_S0K-pi+pi-, K_S0K+pi-pi-, " \
  "K_S0K-pi0, K_S0K_S0pi-, pi-eta, pi-eta'(pi+pi-eta), pi-eta(pi+pi-pi0). " \
  "D_s- → K_S0 K- excluded to avoid double counting in K_S0 K+ measurement.")

alg_ks.signal_side do |s|
  s.charged(kp: 1)
  s.require_charge 1
end

alg_ks.note(:KS0_reconstruction,
  "K_S0 → pi+pi- reconstructed on signal side from unused tracks. " \
  "Secondary vertex fit applied; K_S0 mass window (0.487, 0.511) GeV/c^2.")
alg_ks.note(:extra_track_veto,
  "Events with additional charged tracks satisfying |cosθ|<0.93 and " \
  "|Vz|<20 cm are rejected to suppress combinatorial background.")
alg_ks.note(:dt_fit,
  "2D unbinned maximum likelihood fit to M(K_S0 K+) vs M_tag. " \
  "Signal yield: 1782 ± 47.  Branching fraction: (1.425 ± 0.038_stat)%.")

# No kinematic fit block — the DT signal yield is extracted from a 2D fit
# to M_KS0K+ vs M_tag, not from a constrained kinematic fit.
alg_ks.note(:no_kinematic_fit,
  "No constrained kinematic fit for the K_S0 K+ DT mode. " \
  "Signal yield extracted via 2D unbinned ML fit to invariant mass distributions.")

alg_ks.apply

# ═══════════════════════════════════════════════════════════════════════
# Mode II: D_s+ → K_L0 K+  (ST + missing)
# ═══════════════════════════════════════════════════════════════════════
alg_kl = TagAnalysis.new("Ds_to_KL0_Kplus_ST")

alg_kl.set_header(%w[KinematicFit/KinematicFit.h DTagAlg/DTagAlg.h])
alg_kl.set_constant(ECMS: "MeasuredEcmsSvc::get_ecms()")
alg_kl.with_decay_card(decay_card_kl)

alg_kl.tag_side(:Ds) do |t|
  t.modes(*ds_tag_modes)
  t.mode_group :hadronic
  t.charm(-1)
end

alg_kl.note(:tag_mode_list,
  "Same 13 ST modes as Mode I.  D_s- → K_S0 K- excluded.")

# K_L0 treated as missing particle; K+ on signal side
alg_kl.signal_side do |s|
  s.charged(kp: 1)
  s.require_charge 1
  s.missing :K_L0
end

# 4C kinematic fit: 4-momentum conservation
# gamma_dir from D_s* → gamma D_s selected in the fit
alg_kl.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 40
end

alg_kl.note(:gamma_direct_selection,
  "Photon from D_s* → gamma D_s selected by trying both tag-side and " \
  "signal-side hypotheses; the one with minimum chi2 retained. " \
  "Extra photons with E > 250 MeV and angle to missing particle > 15° rejected.")
alg_kl.note(:missing_mass_squared,
  "MM^2 = (P_e+e- - P_Ds- - P_gamma - P_K+)^2 used to extract K_L0 signal. " \
  "Fit includes peaking backgrounds from D_s+ → K_S0 K+ and D_s+ → eta K+.")
alg_kl.note(:DTagAlg_mass_constraints,
  "Kinematic fit constrains masses of ST D_s-, signal D_s+, and D_s* " \
  "intermediate state, plus initial 4-momenta.")

alg_kl.apply

# ── Execute ───────────────────────────────────────────────────────────
alg_ks.execute_on(datasets + inc_mc + [exclusive_mc_ks])
alg_kl.execute_on(datasets + inc_mc + [exclusive_mc_kl])