# Paper: 2309.05760v2 — Amplitude analysis of D+ -> KS0 pi+ eta
# Method: Double-tag (DT) with 6 hadronic tag modes for D-
# Data: 2.93 fb-1 at psi(3770) (3.773 GeV)
# Result: First observation of D+ -> KS0 a0(980)+; BF(D+ -> KS0 pi+ eta) = (1.27 +/- 0.04 +/- 0.03)%

### Dataset preparation ###
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card: psi(3770) -> D+ D-
# Signal: D+ -> KS0 pi+ eta, eta -> gamma gamma
decay_card = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K_S0 pi+ eta PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Also decay card for D- tag modes (for exclusive MC of the tag side)
# Tag modes: K+pi-pi-, KS0pi-, K+pi-pi-pi0, KS0pi-pi0, KS0pi-pi-pi+, K+K-pi-
decay_card_tag = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Dp_KS0_pi_eta"
  config.related_dataset = psi3770_data
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

# 6 tag modes for D- (from paper)
# Mapped to authoritative DTagAlg channel names:
#   K+pi-pi-     -> DptoKPiPi
#   KS0pi-       -> DptoKsPi
#   K+pi-pi-pi0  -> DptoKPiPiPi0
#   KS0pi-pi0    -> DptoKsPiPi0
#   KS0pi-pi-pi+ -> DptoKsPiPiPi
#   K+K-pi-      -> DptoKKPi
dp_tag_modes = [:DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
                :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi]

### Algorithm: D+ -> KS0 pi+ eta (KS0 -> pi+ pi-, eta -> gamma gamma) ###
alg = TagAnalysis.new("DpKS0PiEta")
alg.set_header(["DpKS0PiEtaAlg/DpKS0PiEta.h"])
alg.set_constant({ "ECMS" => [:double, 3.773] })

# Tag side: D- reconstructed via 6 hadronic modes
alg.tag_side(:D) do |t|
  t.modes(*dp_tag_modes)
  t.charm -1       # tag the D- (charm = -1)
end

# Signal side: D+ -> KS0 pi+ eta
# KS0 -> pi+ pi-, eta -> gamma gamma
# Charged: pi+ (from D+ decay), pi+ and pi- (from KS0) = 2 pi+, 1 pi-
alg.signal_side do |s|
  s.photons 2          # eta -> gamma gamma
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.charged(pip: 2, pim: 1)
  s.require_charge 1   # 2*(+1) + 1*(-1) = +1
end

# 4C kinematic fit with D mass and intermediate resonance constraints
# Invariant masses of (gamma gamma)_eta, (pi+ pi-)_KS0, and D+/- constrained
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D)
  f.chi2_cut 200
end

alg
  .note(:amplitude_analysis,
    "first amplitude analysis of D+ -> KS0 pi+ eta; unbinned maximum-likelihood fit performed in ROOT; " \
    "two dominant amplitudes: D+ -> KS0 a0(980)+ (FF ~ 105%) and D+ -> K0*(1430)0 pi+ (FF ~ 10.8%); " \
    "destructive interference of 15.83 +/- 1.53(stat) +/- 1.65(syst)% observed between the two amplitudes")
  .note(:mBC_window,
    "D+/- candidates required 1.865 < M_BC < 1.875 GeV")
  .note(:deltaE_window,
    "DeltaE windows: tag modes with pi0: -0.055 < DeltaE < 0.040 GeV; " \
    "other tag modes: -0.025 < DeltaE < 0.025 GeV; " \
    "signal side: -0.020 < DeltaE < 0.020 GeV")
  .note(:best_candidate_selection,
    "if multiple tag candidates, the one with minimum |DeltaE| is chosen; " \
    "best signal candidate with minimum |DeltaE| selected at recoiling side")
  .note(:extra_shower_veto,
    "energy of largest unused photons required < 0.23 GeV")
  .note(:dt_yield,
    "1113 DT events obtained with signal purity (98.2 +/- 0.4)% from 2D M_BC fit")
  .note(:a0_line_shape,
    "a0(980)+ propagator parametrized with Flatte formula; parameters fixed to BESIII Ds+ -> pi+ pi0 eta measurement")
  .note(:kstar0_line_shape,
    "K0*(1430)0 propagator parametrized as relativistic Breit-Wigner; parameters from CLEO measurement")
  .note(:other_resonances_tested,
    "also tested but found insignificant: D+ -> K1(1270)0 pi+, K*(1410)+ eta, K2*(1430)+ eta, " \
    "K2*(1580)0 pi+, KS0 pi1(1400)+, KS0 a0(1450)+, KS0 a2(1320)+, K*(892)+ eta, " \
    "(KS0 pi+) S-wave eta (LASS parametrization)")
  .note(:bf_update,
    "BF(D+ -> KS0 pi+ eta) updated to (1.27 +/- 0.04 +/- 0.03)% using amplitude analysis model for efficiency; " \
    "previous measurement was (1.31 +/- 0.05)% with BODY3 generator; " \
    "BF(D+ -> KS0 a0(980)+) = (1.33 +/- 0.05 +/- 0.04)%, " \
    "BF(D+ -> K0*(1430)0 pi+) = (0.14 +/- 0.02 +/- 0.02)%")
  .note(:detection_efficiency,
    "detection efficiency estimated with signal MC generated according to amplitude analysis model; " \
    "uncertainty from amplitude model 0.7%")
  .note(:background_subtraction,
    "background subtracted via sPlot technique on M_BC variable; " \
    "inclusive MC used for background estimation with weight w = (1-purity)*N_D/N_B")

alg.with_decay_card(decay_card).apply
alg.execute_on([psi3770_data, psi3770_incMC, exMC_signal])