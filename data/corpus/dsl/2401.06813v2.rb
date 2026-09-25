# BOSS DSL for 2401.06813v2: Lambda_c+ → n K_S0 pi+ pi0 via double tag (DT)
# Data: e+e- collision at 7 c.m. energies using double tag method
# Tag side: Lambda_c- with 11 ST tag modes
# Signal side: Lambda_c+ → n K_S0 pi+ pi0 (neutron via missing mass)
# Multi-energy: 4599.53, 4611.86, 4628.00, 4640.91, 4661.24, 4681.92, 4698.82 MeV
#
# NOTE: Lambda_c DT is NOT supported by DTagAlg/DTagTool (decayMode() >= 1000 rejected).
# The DSL TagAnalysis framework raises [DSL:tag_dt_lambdac_unsupported] for Lambda_c DT.
# This file approximates the analysis as a TagAnalysis with ST + missing pattern
# and documents all known expressibility gaps.

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_4600 = DatasetManager.real_data.find("703_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")

all_data = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]
inc_mc_4600 = DatasetManager.inclusive_mc.find("703_4600")
inc_mc_4610 = DatasetManager.inclusive_mc.find("706_4610")
inc_mc_4620 = DatasetManager.inclusive_mc.find("706_4620")
inc_mc_4640 = DatasetManager.inclusive_mc.find("706_4640")
inc_mc_4660 = DatasetManager.inclusive_mc.find("706_4660")
inc_mc_4680 = DatasetManager.inclusive_mc.find("706_4680")
inc_mc_4700 = DatasetManager.inclusive_mc.find("706_4700")

all_incMC = [inc_mc_4600, inc_mc_4610, inc_mc_4620, inc_mc_4640, inc_mc_4660, inc_mc_4680, inc_mc_4700]

decay_card = <<~DECAYCARD
  Decay anti-Lambda_c-
  1.0000 anti-p- K_S0 pi+ pi0 PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMCs = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "sig_lambdac_n_ks_pi_pi0"
  config.events        = 1_000_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

alg = TagAnalysis.new("LambdacToNKsPiPi0")

alg.set_header(["LambdacToNKsPiPi0/LambdacToNKsPiPi0.h"])

# Tag side: anti-Lambda_c- with 11 ST modes
# 11 ST tag modes from paper, mapped to DTagAlg symbols:
# 1. p K_S0                → LambdacPtoKsP
# 2. p K+ pi-              → LambdacPtoKPiP
# 3. p K+ pi- pi0          → (unavailable — no known DTagAlg symbol)
# 4. p K_S0 pi- pi+        → (unavailable — no known DTagAlg symbol)
# 5. p K+ pi- pi0          → (same as 3)
# 6. anti-Lambda pi-       → LambdacPtoLambdaPi
# 7. anti-Lambda pi- pi0   → LambdacPtoLambdaPiPi0
# 8. anti-Lambda pi- pi+ pi- → LambdacPtoLambdaPiPiPi
# 9. anti-Sigma- pi0       → LambdacPtoSigmaPiPi0
# 10. anti-Sigma0 pi-      → LambdacPtoSigma0Pi
# 11. anti-Sigma- pi- pi+  → LambdacPtoSigmaPiPi
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP, :LambdacPtoKPiP,
          :LambdacPtoLambdaPi, :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi,
          :LambdacPtoSigmaPiPi0, :LambdacPtoSigma0Pi,
          :LambdacPtoSigmaPiPi
  t.charm -1
end

# Signal side: Lambda_c+ → n K_S0 pi+ pi0
# K_S0 → pi+ pi- (2 charged tracks)
# pi+ (1 charged track from Lambda_c+ decay)
# pi0 → gamma gamma (2 photons)
# n → missing (neutron not detected)
alg.signal_side do |s|
  s.photons 2
  s.missing :n, mass: 0.939565   # neutron mass
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.note(:tag_dt_lambdac_unsupported,
  "DTagTool rejects Lambda_c DT (decayMode() >= 1000). " \
  "This analysis uses the double tag (DT) method for Lambda_c which is " \
  "NOT supported by the DSL TagAnalysis framework. The DSL approximates " \
  "as ST + missing pattern but the actual BOSS implementation uses " \
  "a custom tag-based approach, not standard DTagAlg/DTagTool.")

alg.note(:tag_mode_unavailable,
  "Modes not expressible in DSL: " \
  "pbar K+ pi- pi0 (modes 3 and 5 in paper's Table II — 'pKpi' and 'pK+pi-pi'), " \
  "pbar K_S0 pi- pi+ (mode 4 — 'pKpi-pi+'). " \
  "These have no known DTagAlg channel-name symbols in the frozen BOSS 7.0.6/7.1.2 releases. " \
  "Only 8 of the paper's 11 tag modes are expressible.")

alg.note(:signal_side_ks_reconstruction,
  "Signal-side K_S0 reconstruction (pi+ pi- secondary vertex with L/sigma_L > 2, " \
  "|M(pi+pi-)-M_KS0| mass window) cannot be expressed in the TagAnalysis DSL. " \
  "The DSL signal_side block declares only photon/charged/missing content, " \
  "not intermediate particle reconstruction. The actual BOSS analysis reconstructs " \
  "K_S0 from signal-side pi+ pi- pairs with secondary vertex fit " \
  "before the kinematic fit.")

alg.note(:signal_side_pi0_reconstruction,
  "Signal-side pi0 reconstruction (gamma gamma with 0.115 < M_gg < 0.150 GeV/c^2 " \
  "and 1C mass-constrained fit chi2_1C < 200) cannot be expressed. " \
  "Best pi0 candidate chosen by minimum chi2_1C. " \
  "The DSL Kalman fit is not available inside TagAnalysis signal_side.")

alg.note(:signal_side_charged_track,
  "Signal side has 3 charged tracks: pi+ (from Lambda_c+), pi+ and pi- (from K_S0). " \
  "Additionally, there is a pi+ from Lambda_c+ decay that needs PID: L(pi) > L(K). " \
  "The K_S0 pions have no PID requirement. " \
  "The DSL signal_side charged() multiset cannot distinguish K_S0 tracks from direct tracks. " \
  "Vz < 20 cm for K_S0 daughter tracks, |Vz| < 10 cm for direct pi+.")

alg.note(:peaking_background_vetos,
  "Post-fit peaking background vetoes (applied in ROOT): " \
  "M(n pi-) - M(n) outside (0.22, 0.27) GeV/c^2; " \
  "M(n pi0) - M(n) > 0.20 GeV/c^2; " \
  "M(n pi+) - M(n) outside (0.23, 0.28) GeV/c^2 for both pi+ combinations.")

alg.note(:multi_energy,
  "Analysis uses 7 c.m. energy points: 4599.53, 4611.86, 4628.00, 4640.91, " \
  "4661.24, 4681.92, 4698.82 MeV. " \
  "No ECMS constant is set (multi-energy scan). " \
  "ST yields and DT/ST efficiencies are per-mode per-energy.")

alg.note(:dtag_method,
  "Signal yield obtained by 2D unbinned maximum likelihood fit to M(n) and M(pi+pi-) " \
  "spectra combined from 7 energies. BF computed from DT yields normalized by ST yields. " \
  "Full analysis in ROOT stage, not expressible at BOSS DSL level.")

alg.apply
alg.execute_on(all_data + all_incMC + exMCs)