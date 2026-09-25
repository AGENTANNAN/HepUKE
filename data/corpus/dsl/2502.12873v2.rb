# BESIII DSL: Measurement of strong-phase difference between D0 and D0bar -> K+K-pi+pi-
# in bins of phase space
# Paper: 2502.12873v2
# CMS energy: 3.773 GeV (psi(3770)), luminosity 20.3 fb^-1
# Method: Double-tag with quantum-correlated DD pairs

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================
# Decay cards for exclusive signal MC
# Signal: D0 -> K+ K- pi+ pi-
# Various tag modes listed below
# ============================================================

decay_card_signal = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.000 K+ K- pi+ pi- PHSP;
  Enddecay
  Decay anti-D0
  1.000 K+ pi- PHSP;
  Enddecay
  End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_D0_KKpipi"
  config.related_dataset = data_3773
  config.events          = 2_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# ============================================================
# TagAnalysis: D0 -> K+K-pi+pi-
# Double-tag with quantum-correlated D0 anti-D0 pairs
#
# Tag categories (4 types):
#   Flavor: K-pi+, K-pi+pi0, K-pi+pi-pi+, K-e+nu_e (semileptonic)
#   CP-even: K+K-, pi+pi-, KS0pi0pi0, pi+pi-pi0, KL0pi0
#   CP-odd:  KS0pi0, KS0eta, KS0eta'(pipieta), KS0eta'(rho0gamma), KS0pi+pi-pi0
#   Mixed-CP: KS0pi+pi-, KL0pi+pi-
#
# Many tag modes are custom-reconstructed; standard DTagAlg modes
# used where available.
# ============================================================

alg = TagAnalysis.new("D0toKKpipiStrongPhase")
alg.set_header(["D0toKKpipiStrongPhase/D0toKKpipiStrongPhase.h"])
alg.set_constant({ "ECMS" => [:double, 3.773] })
alg.with_decay_card(decay_card_signal)

# Tag side: D0 with all available tag modes
# Standard DTagAlg modes for flavor tags
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi,
          :D0toKPiPi0,
          :D0toKPiPiPi
  t.charm 0
end

# Signal side: D0 -> K+ K- pi+ pi- (fully reconstructed)
alg.signal_side do |s|
  s.charged(km: 1, kp: 1, pim: 1, pip: 1)
  s.require_charge 0
end

# Standard 4C kinematic fit (no missing particle)
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.apply

alg.note(:quantum_correlation, "Quantum-correlated DD pairs at psi(3770): C=-1 initial state")
alg.note(:binning_scheme, "2x4 bins in phase space (i = +/-1,2,3,4) optimised for gamma measurement")
alg.note(:tag_categories, "Four categories: flavor (standard DTagAlg modes), CP-even, CP-odd, mixed-CP")
alg.note(:additional_tag_modes, <<~NOTE
  Many tag modes are custom-reconstructed outside DTagAlg:
  - CP-even: D0->K+K-, D0->pi+pi-, D0->KS0pi0pi0, D0->pi+pi-pi0, D0->KL0pi0
  - CP-odd:  D0->KS0pi0, D0->KS0eta, D0->KS0eta'(pipieta), D0->KS0eta'(rho0gamma), D0->KS0pi+pi-pi0
  - Mixed-CP: D0->KS0pi+pi-, D0->KL0pi+pi-
  - Semileptonic: D0->K-e+nu_e (flavor)
  These are reconstructed in separate BOSS modules and combined at analysis level.
NOTE
)
alg.note(:partial_reco_tags, <<~NOTE
  Partially reconstructed tag modes with missing particle:
  - D0->KL0pi0 (miss KL0): M_miss^2 peak at KL0 mass, no extra tracks or pi0
  - D0->K-e+nu_e (miss nu_e): U_miss = E_miss - |p_miss| peaks at 0
  - Partially reconstructed signal (miss charged kaon from K+K-pi+pi-):
    tagged by D0->K+K-, D0->KS0pi0, D0->KS0pi+pi-; M_miss^2 peaks at K+/- mass
NOTE
)
alg.note(:ks_veto, "Remove events with pi+pi- invariant mass in [477, 507] MeV/c^2 (KS0 veto)")
alg.note(:dt_selection, "DT events: both D decays reconstructed; best combination by average M_BC closest to D0 mass")
alg.note(:mbc_fit, "ST and DT yields from 1D M_BC fits; signal shape from MC convolved with Gaussian resolution")
alg.note(:strong_phase_fit, <<~NOTE
  Maximum-likelihood fit to DT yields in 8 bins (2x4) using Eqs. (8)-(10):
  Free parameters: c_i, s_i (i=1..4), R_i (recursive K_i fractions),
  B(D0->K+K-pi+pi-), r_D^Kpi cos(delta), r_D^Kpi sin(delta) [Gaussian constrained]
  External inputs: K_S^0,K_L^0 pi+pi- strong phases from BESIII/CLEO combination,
  F_+ fractions for pi+pi-pi0 and KS0pi+pi-pi0, r_D parameters for flavor tags
NOTE
)
alg.note(:cp_even_fraction, "F_+ = 0.754 +/- 0.010(stat) +/- 0.008(syst)")
alg.note(:bf_result, "B(D0->K+K-pi+pi+) = (2.863 +/- 0.028 +/- 0.045) x 10^-3")

alg.execute_on([data_3773, incMC_3773, sig_mc])