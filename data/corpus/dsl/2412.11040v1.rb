# DSL for arXiv:2412.11040v1
# Amplitude analysis and BF measurement of D+ → K- pi+ pi+ pi0
# Using 7.93 fb^-1 at sqrt(s)=3.773 GeV, BESIII
# TagAnalysis: ST (tag D- → K+ pi- pi-) + fully reconstructed signal D+

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# =============================================================================
# Datasets — psi(3770) at 3.773 GeV
# =============================================================================

psipp_data = DatasetManager.real_data.find("712_3773")
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")

# =============================================================================
# Decay card — KKMC + psi(3770) top mother
# D+ → K- pi+ pi+ pi0; pi0 → gamma gamma
# Tag: D- → K+ pi- pi-
# =============================================================================

decay_Dp_Kpipipi0 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-   PHSP;
  Enddecay
  Decay D+
  1.0000 K- pi+ pi+ pi0   PHSP;
  Enddecay
  Decay D-
  1.0000 K+ pi- pi-   PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay
  End
DECAYCARD

# =============================================================================
# Exclusive MC samples
# =============================================================================

# Signal MC — generated with phase space for normalization integral,
# then reweighted by amplitude analysis results for efficiency
exMC_Dp_Kpipipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Dp_Kpipipi0"
  config.related_dataset = psipp_data
  config.events          = 1_000_000
  config.decay_card      = decay_Dp_Kpipipi0
  config.cross_section   = :default
end

# =============================================================================
# TagAnalysis — ST D- tag + fully reconstructed D+ → K- pi+ pi+ pi0
# Tag: D- → K+ pi- pi- (single hadronic tag mode)
# Signal: K-, pi+, pi+, pi0 (pi0 → gamma gamma)
# No missing particle — signal fully reconstructed
# =============================================================================

alg = TagAnalysis.new("DpTagKpipipi0")
alg.set_header(["DpTagKpipipi0Alg/DpTagKpipipi0.h"])
alg.set_constant({ "ECMS" => [:double, 3.773] })
alg.note(:amplitude_analysis, "unbinned maximum-likelihood amplitude analysis using
  isobar model with covariant tensor formalism; 15 intermediate processes identified
  with significances > 5 sigma; dominant component D+[S] → K*(892)0 rho(770)+
  (FF = 66.6%); fit performed at ROOT level with MC integration for normalization")
alg.note(:six_constraint_fit, "6C kinematic fit: 4C (four-momentum conservation) +
  1C (pi0 → gamma gamma mass constraint) + 1C (D+ mass constraint);
  chi2 < 100 retained; fitted four-momenta used in amplitude analysis")
alg.note(:st_yield, "ST D- yield from 1D binned M_BC fit: MC shape convolved with
  double-Gaussian + ARGUS background; N_ST(D-) = 2215326 ± 1589;
  DeltaE window [-0.025, 0.024] GeV applied")
alg.note(:dt_yield, "DT yield from 2D M_BC(sig) vs M_BC(tag) unbinned ML fit;
  N_DT = 35481 ± 220; purity = (98.4±0.1)% in signal region [1.863, 1.879] GeV/c^2;
  26709 signal candidates for amplitude analysis")
alg.note(:D0D0bar_veto, "D0 D0bar background suppressed by rejecting events with
  1.863 < M_BC^W < 1.867 GeV/c^2 under D0 D0bar mispartition hypothesis")
alg.note(:deltaE_cut, "amplitude analysis applies further DeltaE cuts:
  [-0.062, 0.034] GeV for D+, [-0.025, 0.025] GeV for D-")
alg.note(:bf, "B(D+ → K- pi+ pi+ pi0) = (6.06±0.04±0.07)%;
  B(D+ → K*(892)0 rho(770)+) = (4.15±0.07±0.17)% (dominant intermediate process);
  B(D+ → Kbar*(892)0 rho(770)+) = (6.23±0.11±0.25)% (isospin-corrected)")
alg.note(:pi0_reco, "pi0 reconstructed from photon pairs with invariant mass
  [0.115, 0.150] GeV/c^2; 1C kinematic fit constraining to nominal pi0 mass with
  chi2 < 30; at least one photon in barrel EMC required")

# Tag side: D- → K+ pi- pi-
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi      # D- → K+ pi- pi- (charge conjugate of D+ → K- pi+ pi+)
  t.charm -1              # tag D- (anti-D+)
end

# Signal side: K-, pi+, pi+, pi0 (→ gamma gamma)
alg.signal_side do |s|
  s.photons 2                    # two photons for pi0 → gamma gamma
  s.min_photon_angle 10.0        # minimum opening angle to charged tracks > 10 deg
  s.min_photon_energy 0.025      # E > 25 MeV barrel / 50 MeV endcap
  s.charged(km: 1, pip: 2)       # K- and two pi+ from D+ decay
  s.require_charge 1             # D+ total charge: -1 + 1 + 1 = +1
end

# Kinematic fit: 4C + pi0 mass constraint (6C with optional D+ tag constraint)
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Dplus)
  f.chi2_cut 100          # chi2_6C < 100 as stated in paper
end

alg.with_decay_card(decay_Dp_Kpipipi0).apply
alg.execute_on([psipp_data, psipp_incMC, exMC_Dp_Kpipipi0])