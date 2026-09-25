# BOSS DSL for: Amplitude analysis and branching fraction measurement of D0 -> K- pi+ pi0 pi0
# Paper: arXiv:1903.06316v1 (BESIII, 2.93 fb^-1 at sqrt(s)=3.773 GeV at psi(3770))
# DTag technique: tag anti-D0 -> K+ pi- (DTagAlg mode D0toKPi), signal D0 -> K- pi+ pi0 pi0

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card: psi(3770) -> D0 anti-D0, D0 -> K- pi+ pi0 pi0 (signal), anti-D0 -> K+ pi- (tag)
decay_card = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ pi0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the full tag+signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "D0toKPiPi0Pi0_DTag_signal"
  config.related_dataset = data_3773
  config.events = 1_000_000
  config.decay_card = decay_card
  config.cross_section = :default
end

# Tag analysis: tag side anti-D0 -> K+ pi- (DTagAlg), signal side D0 -> K- pi+ pi0 pi0 from remaining tracks/showers
alg = TagAnalysis.new("D0toKPiPi0Pi0_PWA")
alg.set_header(["D0toKPiPi0Pi0_PWAAlg/D0toKPiPi0Pi0_PWA.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })

# Tag side: anti-D0 -> K+ pi- through DTagAlg
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi
  t.charm -1
  # Tag side windows as specified in the paper
  t.window :deltaE, min: -0.03, max: 0.02
  t.window :mBC, min: 1.8575, max: 1.8775
end

# Signal side: D0 -> K- pi+ pi0 pi0 from remaining tracks/showers unused by the tag
alg.signal_side do |s|
  s.photons 4                  # 4 photons for two pi0 -> gamma gamma
  s.charged(km: 1, pip: 1)     # K- and pi+ from the signal D0
  s.min_photon_angle 10.0      # > 10 degrees from charged tracks
  s.min_photon_energy 0.025    # photon energy threshold
  s.require_charge -1          # total signal charge = -1 (K-)
end

# Kinematic fit over the combined tag + signal system
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # 1st pi0 -> gamma gamma
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # 2nd pi0 -> gamma gamma
  f.chi2_cut 80                  # paper's 3C kinematic fit chi2 cut on the signal side
end

# Notes for inexpressible BOSS-side procedures
alg.note(:vertex_fit_signal, "Vertex fit performed on signal side K- and pi+ tracks; chi2 < 100 required before kinematic fit. The tag DSL fit does not support per-track vertex fits on signal charged tracks.")
   .note(:signal_side_windows, "Signal side mBC [1.8600, 1.8730] GeV/c^2 and DeltaE [-0.04, 0.02] GeV applied. In the tag DSL these are stored unconditionally and windowed in ROOT, not applied at BOSS level for the signal side.")
   .note(:pi0_mass_window, "Diphoton invariant mass required in [0.115, 0.150] GeV/c^2 for pi0 candidates; at least one barrel photon required. The tag DSL does not express diphoton mass windows or barrel/endcap photon distinctions.")
   .note(:ks_veto, "K_S0 mass veto: M(pi0 pi0) not in [0.458, 0.520] GeV/c^2 applied on the signal side to suppress D0 -> K- K_S0 pi+ peaking background.")
   .note(:photon_selection, "Photon timing < 700 ns from event start; barrel (|cos(theta)|<0.80) energy > 25 MeV; endcap (0.86<|cos(theta)|<0.92) energy > 50 MeV. The tag DSL's min_photon_energy applies uniformly; barrel/endcap distinction not expressed.")
   .note(:kinematic_fit_strategy, "Paper uses a 3C kinematic fit on the signal side alone (constraining D mass + 2 pi0 masses to PDG values), without 4C constraint to CMS energy. The tag DSL always combines tag+signal with constrain_four_momentum, yielding a different fit configuration (4C + 2 pi0 mass constraints = 6C). The chi2_cut of 80 from the paper's 3C fit is carried over to the tag DSL fit.")
   .note(:amplitude_analysis, "The amplitude analysis (PWA) with 26 intermediate amplitudes, unbinned maximum likelihood fit, covariant tensor formalism, Blatt-Weisskopf barrier factors, and relativistic Breit-Wigner / Gounaris-Sakurai propagators is performed at the ROOT analysis level. The BOSS DSL covers only dataset preparation and event selection up to the kinematic fit.")
   .note(:branching_fraction, "BF measurement: B(D0->K-pi+pi0pi0) = (8.86 +/- 0.13(stat) +/- 0.19(syst))%. ST yield 534,581 +/- 769 from M_BC(K+pi-) fit (CB+Gaussian signal, ARGUS background). DT yield 6,101 +/- 83 from 2D M_BC fit. ST efficiency (66.01+/-0.03)%, corrected DT efficiency (8.50+/-0.04)%. These fits and BF calculation are ROOT-level procedures.")
   .note(:efficiency_corrections, "Data-MC efficiency corrections applied: weighted-average differences of -0.69% (pi0 reconstruction), 1.83% (kaon tracking), 0.22% (pion tracking). PID efficiency difference negligible. These corrections are applied at the ROOT analysis level.")

alg.with_decay_card(decay_card)

alg.apply

alg.execute_on([data_3773, incMC_3773, exMC_signal])