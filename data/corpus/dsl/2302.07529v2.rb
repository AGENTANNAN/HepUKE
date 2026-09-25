# Search for semi-leptonic decays Lambda_c+ -> Lambda pi+ pi- e+ nu_e and
#   Lambda_c+ -> p KS0 pi- e+ nu_e
#   [arXiv:2302.07529]
#
# Tag-based analysis: Lambda_c ST + DT at 4.600-4.699 GeV
# 12 ST tag modes for anti-Lambda_c- (3 unavailable in authoritative tag-mode list)
# Two independent signal modes -> two TagAnalysis objects
# 4.5 fb-1 total integrated luminosity

### Dataset preparation ###
data_4600 = DatasetManager.real_data.find("703_4600")  # 4.600 GeV
data_4612 = DatasetManager.real_data.find("706_4610")  # 4.612 GeV
data_4628 = DatasetManager.real_data.find("706_4620")  # 4.628 GeV
data_4641 = DatasetManager.real_data.find("706_4640")  # 4.641 GeV
data_4661 = DatasetManager.real_data.find("706_4660")  # 4.661 GeV
data_4682 = DatasetManager.real_data.find("706_4680")  # 4.682 GeV
data_4699 = DatasetManager.real_data.find("706_4700")  # 4.699 GeV

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
# incMC for 4.612 not in standard table; use closest

all_data = [data_4600, data_4612, data_4628, data_4641, data_4661, data_4682, data_4699]
all_incMC = [incMC_4600, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

# Tag-side ST modes for anti-Lambda_c- (charge conjugate of Lambda_c+):
# From 12 original modes, 9 are available in the authoritative tag-mode list.
# 3 modes dropped: Sigmabar0 pi-, Sigmabar- pi0, Sigmabar- pi- pi+
# (no corresponding DTagAlg channel symbols exist)

# ---- Decay card: Signal mode 1 (Lambda_c+ -> Lambda pi+ pi- e+ nu_e) ----
decay_card_mode1 = <<~DECAYCARD
  Decay vpho
  1.0000 Lambda_c+ anti-Lambda_c-   PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Lambda pi+ pi- e+ nu_e   PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K+ pi-   PHSP;
  Enddecay

  Decay Lambda
  1.0000 p+ pi-   PHSP;
  Enddecay

  End
DECAYCARD

# ---- Decay card: Signal mode 2 (Lambda_c+ -> p KS0 pi- e+ nu_e) ----
decay_card_mode2 = <<~DECAYCARD
  Decay vpho
  1.0000 Lambda_c+ anti-Lambda_c-   PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 p+ K_S0 pi- e+ nu_e   PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K+ pi-   PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-   PHSP;
  Enddecay

  End
DECAYCARD

exMC_mode1 = DatasetManager.create_exclusive_mc_for(all_data) do |c|
  c.sample_name     = "Lc_to_Lambda_pipi_enu"
  c.events          = 200_000
  c.decay_card      = decay_card_mode1
  c.cross_section   = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc_for(all_data) do |c|
  c.sample_name     = "Lc_to_pKSpi_enu"
  c.events          = 200_000
  c.decay_card      = decay_card_mode2
  c.cross_section   = :default
end

# =============================================================================
# TAG ANALYSIS 1: Lambda_c+ -> Lambda pi+ pi- e+ nu_e
# =============================================================================
alg1 = TagAnalysis.new("LcToLambdaPiPiENu")

# Multi-energy: skip ECMS
alg1.set_header(["LcToLambdaPiPiENuAlg/LcToLambdaPiPiENu.h"])
    .with_decay_card(decay_card_mode1)

# Tag side: anti-Lambda_c- reconstructed from 9 ST modes
# charm(-1) = tag the anti-Lambda_c
alg1.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP, :LambdacPtoKPiP, :LambdacPtoKsPi0P,
          :LambdacPtoKsPiPiP, :LambdacPtoKPiPi0P,
          :LambdacPtoLambdaPi, :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi, :LambdacPtoPiPiP
  t.charm(-1)
end

# Signal side: Lambda_c+ -> Lambda pi+ pi- e+ nu_e
# Lambda -> p pi-
# Signal particles: p + pi+(bachelor) + pi-(bachelor) + e+ + nu_e
#   + p(from Lambda) + pi-(from Lambda) = 1 proton, 1 pi+, 2 pi-, 1 e+
# Total 5 charged tracks, neutrino missing
alg1.signal_side do |s|
  s.photons 0
  s.charged(prp: 1, pip: 1, pim: 2, ep: 1)
  s.missing :nu_e           # massless neutrino
  s.require_charge 1        # Lambda_c+ charge = +1
end

# Kinematic fit
alg1.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg1.apply

alg1.note(:tag_mode_unavailable,
        "Three ST modes from the original paper are not available in the DTagAlg " \
        "channel list: Sigmabar0 pi-, Sigmabar- pi0, and Sigmabar- pi- pi+. " \
        "Only 9 of 12 ST modes are used in this spec. The original paper reports " \
        "a total ST yield of 123,509 +/- 461 across all 12 modes.")
  .note(:st_selection,
        "ST anti-Lambda_c- candidates selected using M_BC and Delta_E. " \
        "Delta_E requirements are mode-dependent (see Ref. [33]). " \
        "If multiple ST candidates, the one with minimum |Delta_E| is retained. " \
        "M_BC signal region around Lambda_c mass; sideband [2.25, 2.27] GeV/c^2 " \
        "used for background estimation.")
  .note(:signal_selection_mode1,
        "Signal: Lambda_c+ -> Lambda pi+ pi- e+ nu_e. " \
        "Five charged tracks on signal side: p + pi+(bachelor) + pi-(bachelor) + " \
        "p(from Lambda) + pi-(from Lambda) + e+. " \
        "Lambda candidates: p pi- from secondary vertex fit (chi2 < 100), " \
        "M(p pi-) in [1.09, 1.14] GeV/c^2, positive decay length. " \
        "Loose tracks for Lambda daughters: |Vz| < 20 cm, no Vr requirement. " \
        "Tight tracks for others: |Vz| < 10 cm, Vr < 1 cm, |cos(theta)| < 0.93.")
  .note(:pid_requirements_mode1,
        "Proton PID: L(p) > L(K), L(p) > L(pi), L(p) > 0. " \
        "Pion PID: L(pi) > L(K), L(pi) > 0. " \
        "Positron PID: L(e) > 0.001 and L(e)/(L(e)+L(pi)+L(K)) > 0.99. " \
        "FSR recovery: showers within 5 degree cone of positron track added back.")
  .note(:background_suppression_mode1,
        "cos_theta(e,pi) < 0.88 to suppress gamma-conversion. " \
        "M(Lambda pi+ pi- e(pi)+) < 2.20 GeV/c^2 to suppress " \
        "Lambda_c+ -> Lambda pi+ pi- pi+. " \
        "cos_theta(P_miss, gamma_most_energetic) < 0.82 to suppress " \
        "Lambda_c+ -> Lambda pi+ omega/eta backgrounds. " \
        "U_miss = E_miss - |p_miss| in [-0.08, 0.08] GeV.")
  .note(:signal_extraction_mode1,
        "3 observed events in signal region, 9 in M_BC sideband. " \
        "No significant signal observed. Upper limit at 90% CL: " \
        "BF < 3.9 x 10^-4. DT efficiency = 9.69 +/- 0.03%.")

alg1.execute_on(all_data + all_incMC + exMC_mode1)

# =============================================================================
# TAG ANALYSIS 2: Lambda_c+ -> p KS0 pi- e+ nu_e
# =============================================================================
alg2 = TagAnalysis.new("LcToPKS0PiENu")

alg2.set_header(["LcToPKS0PiENuAlg/LcToPKS0PiENu.h"])
    .with_decay_card(decay_card_mode2)

# Same tag side as mode 1
alg2.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP, :LambdacPtoKPiP, :LambdacPtoKsPi0P,
          :LambdacPtoKsPiPiP, :LambdacPtoKPiPi0P,
          :LambdacPtoLambdaPi, :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi, :LambdacPtoPiPiP
  t.charm(-1)
end

# Signal side: Lambda_c+ -> p KS0 pi- e+ nu_e
# KS0 -> pi+ pi-
# Signal particles: p + pi+(from KS0) + pi-(from KS0) + pi-(bachelor) + e+ + nu_e
# Total 5 charged tracks + missing neutrino
alg2.signal_side do |s|
  s.photons 0
  s.charged(prp: 1, pip: 1, pim: 2, ep: 1)
  s.missing :nu_e
  s.require_charge 1
end

alg2.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg2.apply

alg2.note(:tag_mode_unavailable,
        "Same as alg1: 3 ST modes (Sigmabar0 pi-, Sigmabar- pi0, Sigmabar- pi- pi+) " \
        "not available in DTagAlg channel list.")
  .note(:st_selection,
        "Same ST selection as mode 1. Tag-side windows from Ref. [33]. " \
        "ST yield: 123,509 +/- 461 total across 12 modes at all energies.")
  .note(:signal_selection_mode2,
        "Signal: Lambda_c+ -> p KS0 pi- e+ nu_e. " \
        "Five charged tracks on signal side: p + pi+(from KS0) + pi-(from KS0) + " \
        "pi-(bachelor) + e+. " \
        "KS0 candidates: pi+ pi- from secondary vertex fit (chi2 < 100), " \
        "M(pi+ pi-) in [0.490, 0.504] GeV/c^2, positive decay length. " \
        "Loose tracks for KS0 daughters: |Vz| < 20 cm, no Vr requirement. " \
        "Tight tracks for others: |Vz| < 10 cm, Vr < 1 cm, |cos(theta)| < 0.93.")
  .note(:pid_requirements_mode2,
        "Same PID as mode 1 except: e+ PID uses L(e)/(L(e)+L(pi)+L(K)) > 0.98 " \
        "(slightly looser than mode 1's 0.99).")
  .note(:background_suppression_mode2,
        "cos_theta(e,pi) < 0.92 to suppress gamma-conversion. " \
        "M(p KS0 pi- e(pi)+) < 2.28 GeV/c^2 to suppress " \
        "Lambda_c+ -> p KS0 pi+ pi-. " \
        "cos_theta(P_miss, gamma_most_energetic) < 0.90 to suppress " \
        "Lambda_c+ -> p KS0 eta backgrounds. " \
        "U_miss = E_miss - |p_miss| in [-0.08, 0.08] GeV.")
  .note(:signal_extraction_mode2,
        "2 observed events in signal region, 0 in M_BC sideband. " \
        "No significant signal observed. Upper limit at 90% CL: " \
        "BF < 3.3 x 10^-4. DT efficiency = 13.58 +/- 0.02%.")
  .note(:systematics,
        "Systematic uncertainties (Table 3): tracking (pi 2.6%, e 0.5%), " \
        "PID (pi 0.7%, e 2.8%), Lambda reconstruction 2.2% (mode 1), " \
        "KS0 reconstruction 3.2% (mode 2), signal model 2.2-5.6%. " \
        "Total 5.2% (mode 1), 7.5% (mode 2).")

alg2.execute_on(all_data + all_incMC + exMC_mode2)