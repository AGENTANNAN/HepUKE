# Evidence for the semileptonic decays Lambda_c+ -> Sigma+/- pi-/+ e+ nu_e
#   [arXiv:2512.05178]
#
# Double-tag (DT) analysis at sqrt(s) = 4.600-4.699 GeV with 4.5 fb^-1.
# ST Lambda_c- reconstructed from 12 hadronic tag modes; signal Lambda_c+ -> Sigma+/- pi-/+ e+ nu_e
# searched for in the tracks recoiling against the ST.
# This is a Lambda_c tag analysis: one tag side (ST Lambda_c-bar) + signal side + missing neutrino.

### Dataset description ###
# Lambda_c+ scan data at sqrt(s) = 4.600, 4.612, 4.628, 4.641, 4.661, 4.682, 4.699 GeV
data_4600 = DatasetManager.real_data.find("703_4600")
data_4612 = DatasetManager.real_data.find("706_4610")  # sample 4612
data_4628 = DatasetManager.real_data.find("706_4620")  # sample 4628
data_4641 = DatasetManager.real_data.find("706_4640")  # sample 4640
data_4661 = DatasetManager.real_data.find("706_4660")
data_4682 = DatasetManager.real_data.find("706_4680")
data_4699 = DatasetManager.real_data.find("706_4700")

lambda_c_scan_data = [data_4600, data_4612, data_4628, data_4641, data_4661, data_4682, data_4699]

# Inclusive MC for each energy point
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4612 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4628 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4641 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4661 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4682 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4699 = DatasetManager.inclusive_mc.find("706_4700")

lambda_c_scan_incMC = [incMC_4600, incMC_4612, incMC_4628, incMC_4641, incMC_4661, incMC_4682, incMC_4699]

# ---------------------------------------------------------------- decay cards
# Signal decay card: e+e- -> Lambda_c+ Lambda_c-,
#   Lambda_c- -> tag modes, Lambda_c+ -> Sigma+ pi- e+ nu_e (and Sigma- pi+ e+ nu_e)
# Sub-decay blocks
sub_pi0_gg = <<~DECAYCARD
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
DECAYCARD

sub_ks_pipi = <<~DECAYCARD
  Decay K_S0
  1.0000 pi+ pi-  PHSP;
  Enddecay
DECAYCARD

sub_lambda_ppi = <<~DECAYCARD
  Decay Lambda0
  1.0000 p+ pi-  HypWK;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+  HypWK;
  Enddecay
DECAYCARD

sub_sigma_decays = <<~DECAYCARD
  Decay Sigma+
  0.5160 p+ pi0  PHSP;
  0.4840 n0 pi+  PHSP;
  Enddecay

  Decay anti-Sigma-
  1.0000 anti-p- pi0  PHSP;
  Enddecay

  Decay Sigma-
  1.0000 n0 pi-  PHSP;
  Enddecay

  Decay anti-Sigma+
  1.0000 anti-n0 pi+  PHSP;
  Enddecay
DECAYCARD

sub_sigma0_lamgam = <<~DECAYCARD
  Decay Sigma0
  1.0000 Lambda0 gamma  PHSP;
  Enddecay

  Decay anti-Sigma0
  1.0000 anti-Lambda0 gamma  PHSP;
  Enddecay
DECAYCARD

# Combined decay card for the signal process:
# Lambda_c+ -> Sigma+ pi- e+ nu_e with Sigma+ -> p pi0 (signal mode 1)
# Lambda_c- -> tag side (inclusive, modelled by all hadronic modes)
decay_card_signal_mode1 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c-  PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Sigma+ pi- e+ nu_e  PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K_S0  PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0  PHSP;
  Enddecay

  #{sub_ks_pipi}
  #{sub_pi0_gg}

  End
DECAYCARD

# Signal MC for mode 2+3 (neutron channels): Lambda_c+ -> Sigma+/- pi-/+ e+ nu_e with Sigma -> n pi
decay_card_signal_mode23 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c-  PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Sigma+ pi- e+ nu_e  PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K_S0  PHSP;
  Enddecay

  Decay Sigma+
  1.0000 n0 pi+  PHSP;
  Enddecay

  #{sub_ks_pipi}

  End
DECAYCARD

# Generate exclusive MC for Lambda_c scan energy points
# Mode 1: Lambda_c+ -> Sigma+ pi- e+ nu_e, Sigma+ -> p pi0
exMC_signal_mode1 = DatasetManager.create_exclusive_mc_for(lambda_c_scan_data) do |config|
  config.sample_name   = "Lc_Sigma_pepi_nu"
  config.events        = 200_000
  config.decay_card    = decay_card_signal_mode1
  config.cross_section = :default
end

# Mode 2+3: Lambda_c+ -> Sigma+/- pi-/+ e+ nu_e, Sigma+/- -> n pi+/-
exMC_signal_mode23 = DatasetManager.create_exclusive_mc_for(lambda_c_scan_data) do |config|
  config.sample_name   = "Lc_Sigma_npi_epi_nu"
  config.events        = 200_000
  config.decay_card    = decay_card_signal_mode23
  config.cross_section = :default
end

### Event selection (BOSS) — TagAnalysis (Lambda_c tag) ###
alg_name = "LcSemiLeptonicSigmaPiENu"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.600]})   # nominal; per-run measured Ecms from MeasuredEcmsSvc
   .with_decay_card(decay_card_signal_mode1)

# Tag side — ST Lambda_c-bar reconstructed via 12 hadronic decay modes.
# Both charges are scanned (charm-less findSTag overload) to cover Lambda_c+ and Lambda_c-.
alg.tag_side(:Lambdac) do |t|
  t.mode_group :hadronic
end

# Signal side — the particles recoiling against the tag:
# Mode 1: Sigma+ -> p pi0, final state p, pi-, e+ + missing nu_e (pi0 partially reconstructed)
# Modes 2+3: Sigma+/- -> n pi+/-, final state n, pi+/-, e+ + missing nu_e (neutron narrow shower)
# The signal_side declares the minimum reconstruction: pion, positron + missing neutrino.
# The proton/neutron and pi0 reconstruction are handled per-submode at ROOT level.
alg.signal_side do |s|
  s.photons 0                           # pi0 reconstructed at ROOT level
  s.charged(pim: 1, ep: 1, at_least: true)  # common: one pion, one positron; at_least allows extra proton/pion
  s.require_charge 0                    # -1 (pim) + 1 (ep) = 0 net charge for SL decay
  s.missing :nu_e                       # massless neutrino
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# Kinematic fit: 4C (tag + pi + e + nu = ecms_lab)
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.note(:tag_modes,
         "Twelve hadronic ST modes for Lambda_c-bar used: " \
         "pbar K_S0, pbar K+ pi-, pbar K_S0 pi0, pbar K_S0 pi- pi+, " \
         "pbar K+ pi- pi0, Lambdabar pi-, Lambdabar pi- pi0, " \
         "Lambdabar pi- pi+ pi-, Sigma0bar pi-, Sigma-bar pi0, " \
         "Sigma-bar pi- pi+, pbar pi- pi+. " \
         "If multiple ST candidates per tag mode per charge per event, " \
         "the candidate with minimum |DeltaE| is retained. " \
         "ST yields and efficiencies obtained by fits to M_BC distributions.")
  .note(:signal_submodes,
         "Three signal reconstruction modes are analyzed together in a simultaneous fit: " \
         "(1) Lambda_c+ -> Sigma+ pi- e+ nu_e, Sigma+ -> p pi0 — requires 3 charged tracks (p, pi-, e+), " \
         "pi0 -> gamma gamma, and M(p pi0) in [1.176, 1.200] GeV/c^2 with kinematic mass constraint; " \
         "(2) Lambda_c+ -> Sigma+ pi- e+ nu_e, Sigma+ -> n pi+ — uses EMC neutron shower with U_miss=0 constraint; " \
         "(3) Lambda_c+ -> Sigma- pi+ e+ nu_e, Sigma- -> n pi- — same neutron method. " \
         "Mode 1 uses M_miss^2 for fitting; modes 2+3 use U_miss. The BFs of the two SL modes " \
         "are assumed equal under isospin symmetry in the simultaneous fit.")
  .note(:background_suppression,
         "Background suppression requirements applied at ROOT level: " \
         "(a) M(p pi-) > 1.13 GeV/c^2 to suppress Lambda_c+ -> Lambda e+ nu_e, Lambda -> p pi-; " \
         "(b) M(Sigma+ pi- pi(e)+) < 2.27 GeV/c^2 to suppress Lambda_c+ -> Sigma+ pi+ pi-; " \
         "(c) cos(theta_pi,e) < 0.95 to suppress gamma-conversion backgrounds; " \
         "(d) pi0 veto for neutron modes: showers with E3x3/E5x5 > 0.9 and second moment < 20 cm^2, " \
         "M(gamma gamma) in [0.115, 0.150], chi2_1C < 20; " \
         "(e) Mass windows: M(n pi+) in [1.15, 1.23], M(n pi-) in [1.16, 1.24] GeV/c^2.")
  .note(:simultaneous_fit,
         "Simultaneous unbinned maximum likelihood fit to M_miss^2 (mode 1) and U_miss (modes 2+3) " \
         "distributions. Signal shapes from signal MC; peaking backgrounds (Lambda_c+ -> Sigma pi pi, " \
         "Lambda_c+ -> Sigma omega) from dedicated MC with yields fixed by known BFs. " \
         "Cross-feed between neutron modes estimated from signal MC (~3%). " \
         "Additive systematic uncertainties from fit variations: 5.6% for simultaneous fit.")
  .note(:mc_model,
         "Signal MC generated uniformly in phase space. Systematic uncertainty from MC model " \
         "estimated by reweighting to Lambda(1405) + Lambda(1520) resonance model from quark model " \
         "predictions (Ref.[16]). Differences in efficiencies: 16.7% (proton mode), 1.3% (Sigma+ n pi+), " \
         "3.9% (Sigma- n pi-).")
  .note(:neutron_reconstruction,
         "Neutron-induced showers in EMC identified by largest deposited energy among candidates " \
         "satisfying U_miss = 0 constraint for momentum magnitude. Systematic uncertainty 8.0% " \
         "from control sample Lambda_c+ -> Sigma+ pi+ pi-, Sigma+ -> n pi+. " \
         "pi0 veto systematic uncertainty 11.7% from same control sample.")
  .note(:energy_scan,
         "Data taken at seven energy points (4.600, 4.612, 4.628, 4.641, 4.661, 4.682, 4.699 GeV). " \
         "ST yields are summed over all energy points; efficiencies averaged by yield-weighted mean. " \
         "Total ST yield N_ST = 120350 +/- 464 (statistical only). The Lambda_c+ momentum is " \
         "determined from the ST Lambda_c-bar direction and beam energy.")
  .note(:bf_calculation,
         "BF = N_DT / (B_inter * N_ST * epsilon_s). Averaged signal efficiency epsilon_s " \
         "= (1/N_ST) * sum_{i,j} (N_ST_ij / epsilon_ST_ij * epsilon_DT_ij), " \
         "where i = ST mode, j = energy point. Missing mass defined as " \
         "E_miss = E_beam - E_SL, p_miss = p_Lambda_c+ - p_SL; " \
         "U_miss = E_miss - c|p_miss|, M_miss^2 = E_miss^2/c^4 - |p_miss|^2/c^2. " \
         "Lambda_c+ momentum: p_Lambda_c+ = -p_hat_tag * sqrt(E_beam^2/c^2 - m_Lambdac^2 c^2).")

alg.apply
alg.execute_on(lambda_c_scan_data + lambda_c_scan_incMC + exMC_signal_mode1 + exMC_signal_mode23)