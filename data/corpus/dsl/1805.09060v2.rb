#==============================================================================
# Paper: 1805.09060v2
# Title: Measurement of the absolute branching fraction of the inclusive
#        semileptonic Lambda_c+ decay with the BESIII detector
# Journal: Phys. Rev. Lett. 121, 062003 (2018)
# Data: sqrt(s) = 4.6 GeV, 567 pb^-1, BOSS 703
# Method: Double-tag (ST + missing nu) with Lambda_c- hadronic tags
# Signal: Lambda_c+ -> X e+ nu_e (inclusive)
# Tag modes: Lambda_c- -> pbar K_S^0, Lambda_c- -> pbar K^+ pi^-
#==============================================================================

# --- Decay card: e+e- -> Lambda_c+ anti-Lambda_c- at 4.6 GeV ---
decay_card = <<~DECAY
Decay vpho
1.0 Lambda_c+ anti-Lambda_c- PHSP;
Enddecay
DECAY

# --- Algorithm ---
alg = TagAnalysis.new("LambdacIncSL")

alg.set_header(["LambdacIncSLAlg/LambdacIncSL.h"])
    .set_constant({ "ECMS" => [:double, 4.600] })
    .with_decay_card(decay_card)

# --- Tag side: Lambda_c- reconstructed via two hadronic modes ---
# Both modes are charge-conjugated by charm(-1) to select anti-Lambda_c
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP, :LambdacPtoKPiP
  t.charm -1
end

# --- Signal side: inclusive semileptonic ---
# Only the positron is explicitly identified; the hadronic system X
# produces extra tracks that must be tolerated via at_least: true.
# The neutrino is treated as a missing massless particle (2-arg p4 overload).
# Auto-stored observables: m_P4_miss_fit, m_Umiss, m_Umiss2, m_q2
alg.signal_side do |s|
  s.charged(ep: 1, at_least: true)
  s.missing :nu_e
end

# --- Kinematic fit: 4C (tag + e+ + nu = ecms_lab) ---
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.note(:pid_unfolding,
  "Positron PID uses CL(e) > 0.001 and CL(e)/(CL(e)+CL(pi)+CL(K)+CL(p)) > 0.8 " \
  "plus E_e/p_e > 0.8. Hadron misidentification unfolded via PID efficiency " \
  "matrix in ROOT stage. Not expressible in BOSS DSL.")

alg.note(:tracking_efficiency,
  "Tracking efficiency corrected via unfolding matrix N_true = sum_j T(i|j) N_prod(j) " \
  "in ROOT stage.")

alg.note(:momentum_extrapolation,
  "Positron momentum spectrum below 200 MeV/c extrapolated via fit to exclusive " \
  "semileptonic mode templates. Fraction f(p_e < 200 MeV/c) = (5.6 +/- 1.5)% applied " \
  "in ROOT stage.")

alg.note(:ws_subtraction,
  "Secondary positrons (gamma conversions, pi0 Dalitz) subtracted via wrong-sign " \
  "sample. WS subtraction validated in MC. M_BC sideband subtraction with scale " \
  "factor 0.78 applied in ROOT stage.")

alg.apply

# --- Datasets ---
data_4600  = DatasetManager.real_data.find("703_4600")
inc_mc_4600 = DatasetManager.inclusive_mc.find("703_4600")

sig_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_lambdac_inc_sl"
  config.related_dataset = data_4600
  config.events          = 500_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# --- Execute ---
alg.execute_on([data_4600, inc_mc_4600, sig_mc])