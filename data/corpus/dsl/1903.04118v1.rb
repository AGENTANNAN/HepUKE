#
# 1903.04118v1: Amplitude analysis of D_s+ → pi+ pi0 eta
#   at sqrt(s) = 4.178 GeV (3.19 fb⁻¹)
#   Double-tag (DT) analysis with 7 ST modes for D_s-
# TagAnalysis — tag-based
#

# ── Dataset ────────────────────────────────────────────────────────────
datasets  = [DatasetManager.load_real_data.find("703_4180")]
inc_mc    = [DatasetManager.load_inclusive_mc.find("703_4180")]

# ── Exclusive MC ──────────────────────────────────────────────────────
# conexc generator for open-charm production at 4.178 GeV

decay_card = <<~DECAY
  Decay D*-
  1.0 gamma D_s-  VSP_PWAVE;
  Decay D_s+
  1.0 pi+ pi0 eta  PHSP;
  Decay pi0
  1.0 gamma gamma  PHSP;
  Decay eta
  1.0 gamma gamma  PHSP;
  Enddecay
DECAY

exclusive_mc = DatasetManager.create_exclusive_mc do |c|
  c.decay_card(decay_card, generator: :conexc)
end

# ── TagAnalysis ────────────────────────────────────────────────────────
alg = TagAnalysis.new("Ds_to_pi_pi0_eta_DT")

alg.set_header(%w[KinematicFit/KinematicFit.h DTagAlg/DTagAlg.h])
alg.set_constant(ECMS: "MeasuredEcmsSvc::get_ecms()")
alg.with_decay_card(decay_card)

# ── Tag side (D_s-) — 7 ST modes ──────────────────────────────────────
alg.tag_side(:Ds) do |t|
  t.modes(:DstoKsKm, :DstoKpKmPim, :DstoKsKmPi0, :DstoKpKmPimPi0,
          :DstoKsKpPimPim, :DstoPimEta, :DstoPimEtap)
  t.charm(-1)
end

# ── Signal side (D_s+) ────────────────────────────────────────────────
alg.signal_side do |s|
  s.photons 4               # pi0 → yy, eta → yy (4 photons)
  s.charged(pip: 1)         # one pi+
  s.require_charge 1
end

# ── 7C kinematic fit: 4-momentum + pi0 mass + eta mass + D_s mass ─────
# The tag side particles enter automatically.
# Signal side: pi+, 2 gamma (pi0), 2 gamma (eta)
# The direct photon (gamma_dir) from D_s* → gamma D_s is also fitted.
# NOTE: the 7C fit with gamma_dir and D_s mass constraint is complex;
# the DSL captures the core constraints.
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 1000
  f.store_fitted_momenta
end

alg.note(:gamma_direct,
  "Direct photon from D_s* → gamma D_s selected in 7C/7CA fit. " \
  "The photon yielding smallest chi2_7C is selected as gamma_direct.")
alg.note(:D_s_mass_constraint,
  "Signal D_s+ mass constrained to nominal D_s mass in 7C fit. " \
  "7CA fit constrains gamma_dir + D_s(+/-) to D_s* mass. " \
  "These multi-constraint fits require custom C++ code beyond v1 DSL.")
alg.note(:BDT,
  "Boosted Decision Tree classifier applied in ROOT stage using " \
  "M(gamma_gamma)_eta, p(low-E_gamma)_eta, p(gamma_dir) as inputs. " \
  "Retains 77.8% signal, rejects 73.4% background.")
alg.note(:cos_theta_eta_veto, "Events with cos(theta_eta) < 0.998 vetoed.")
alg.note(:amplitude_analysis, "Amplitude analysis performed in ROOT stage on accepted DT candidates.")

alg.apply

# ── Execute ───────────────────────────────────────────────────────────
alg.execute_on(datasets + inc_mc + [exclusive_mc])