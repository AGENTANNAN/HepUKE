# ============================================================
# Paper: 2005.11211v2
# Measurement of absolute branching fraction of inclusive decay
#   Lambda_c+ -> K_S^0 X
# Double-tag technique at sqrt(s) = 4.6 GeV
# ============================================================

data_4600 = DatasetManager.real_data.find("703_4600")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.0000 anti-Lambda_c- Lambda_c+  PHSP;
    Enddecay
    Decay anti-Lambda_c-
    1.0000 anti-p- K+ pi-  PHSP;
    Enddecay
    Decay Lambda_c+
    1.0000 K_S0 X  PHSP;
    Enddecay
    Decay K_S0
    1.0000 pi+ pi-  PHSP;
    Enddecay
    End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Lc_KsX"
  config.related_dataset = data_4600
  config.events          = 500000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

alg = TagAnalysis.new("LcKsX")
alg.set_header(["LcKsXAlg/LcKsX.h"])
    .set_constant({ "ECMS" => [:double, 4.600] })
    .with_decay_card(decay_card)

# Single tag (ST): anti-Lambda_c- (charm -1) via hadronic decay modes
alg.tag_side(:Lambdac) do |t|
  # Available modes from the authoritative tag-mode table
  t.modes :LambdacPtoKsP,        # pbar K_S^0
          :LambdacPtoKPiP,       # pbar K+ pi-
          :LambdacPtoKsPi0P,     # pbar K_S^0 pi0
          :LambdacPtoKsPiPiP,    # pbar K_S^0 pi+ pi-
          :LambdacPtoKPiPi0P,    # pbar K+ pi- pi0
          :LambdacPtoLambdaPi,   # Lambda pi-
          :LambdacPtoLambdaPiPi0,    # Lambda pi- pi0
          :LambdacPtoLambdaPiPiPi    # Lambda pi- pi+ pi-
  t.charm -1
end

# Signal side: K_S^0 -> pi+ pi- in the recoil system (inclusive Lambda_c+ -> K_S^0 X)
alg.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.require_charge 0
  # Inclusive: the rest of Lambda_c+ decay products (X) are undetected.
  # Use massless missing particle to account for the unseen remainder.
  s.missing :X0, mass: nil
end

alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 999999   # effectively no cut; K_S^0 mass selection via ROOT 2D fit
end

alg
  .note(:tag_mode_unavailable, "Three ST modes from the paper are not in the
    authoritative tag-mode table: Sigma_bar0 pi-, Sigma- pi0, Sigma- pi+ pi-.
    These 3 modes correspond to ~20% of the total ST BF and are excluded from
    the DSL spec. The paper uses 11 modes total (8 available here).")
  .note(:ks0_signal_side, "K_S^0 candidates are selected among tracks remaining
    after ST Lambda_c- tag with the same vertex-fit criteria as the tag side:
    pi+pi- vertex fit chi2 < 100, decay vertex separated from IP by >2*sigma,
    |cos(theta)| < 0.93, Vz < 20 cm, no Vr constraint. If multiple K_S^0
    candidates, the one with minimum vertex chi2 is kept.")
  .note(:two_dim_fit, "Signal yield extracted from 2D unbinned ML fit to
    M_BC (beam-constrained mass from ST tag) vs M(pi+pi-). Signal function is
    product of Lambda_c- signal (MC shape convolved with Gaussian) and K_S^0
    signal (Gaussian). Three background components modeled. The absolute BF:
    B(Lc+ -> K_S^0 X) = N_sig / [B(K_S0->pi+pi-) * Sum_i(N_tag_i * eps_DT_i / eps_ST_i)].")
  .note(:background_veto, "Vetoes applied in ROOT:
    - Lambda veto for pbar K_S^0 pi0, pbar K_S^0 pi+ pi-, Sigma- pi+ pi- modes:
      M(pbar pi+) in [1.110, 1.120] rejected
    - K_S^0 veto for Lambda pi- pi+ pi-, Sigma- pi0, Sigma- pi+ pi- modes:
      M(pi+pi-) or M(pi0pi0) in [0.480, 0.520] rejected
    - Sigma- veto for pbar K_S^0 pi0 mode:
      M(pbar pi0) in [1.170, 1.200] rejected")
  .note(:systematic_uncertainties, "Systematic sources: ST-related 1.2%,
    K_S^0 reconstruction 1.5%, B(K_S0->pi+pi-) 0.1%, signal yield 3.4%.
    Total 3.9%.")

alg.apply
alg.execute_on([data_4600, incMC_4600, exMC])