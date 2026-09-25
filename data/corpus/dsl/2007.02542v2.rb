# ============================================================
# Paper: 2007.02542v2
# Measurement of D -> omega pi pi and D -> eta pi pi
# Double-tag technique at psi(3770), 2.93 fb^-1
# omega/eta -> pi+ pi- pi0; D0 -> omega/eta pi+ pi- (pi0 pi0)
#                           D+ -> omega/eta pi+ pi0
# ============================================================

data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Representative decay card: psi(3770) -> D0 anti-D0,
# signal D0 -> omega pi+ pi-, omega -> pi+ pi- pi0
decay_card = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0  PHSP;
    Enddecay
    Decay D0
    1.0000 omega pi+ pi-   PHSP;
    Enddecay
    Decay anti-D0
    1.0000 K- pi+   PHSP;
    Enddecay
    Decay omega
    1.0000 pi+ pi- pi0   PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma   PHSP;
    Enddecay
    End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_D0_omega_pipi"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# ============================================================
# Process 1: D0 -> omega/eta pi+ pi-  &  D0 -> omega/eta pi0 pi0
# Tag side: anti-D0 via hadronic modes
# Signal side: D0 -> (omega/eta) + (pi+ pi- / pi0 pi0)
# ============================================================

alg_D0 = TagAnalysis.new("D0OmegaEtaPiPi")
alg_D0.set_header(["D0OmegaEtaPiPiAlg/D0OmegaEtaPiPi.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .with_decay_card(decay_card)

alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi,        # anti-D0 -> K+ pi-
          :D0toKPiPi0,     # anti-D0 -> K+ pi- pi0
          :D0toKPiPiPi     # anti-D0 -> K+ pi- pi- pi+
  t.charm -1               # tag anti-D0
end

alg_D0.signal_side do |s|
  s.charged(pip: 1, pim: 1, at_least: true)
  s.photons 2..6
  s.require_charge 0
end

alg_D0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_D0
  .note(:signal_modes, "Four signal modes measured from remaining tracks:
    (1) D0 -> omega pi+ pi-  (pip:2, pim:2, photons:2)
    (2) D0 -> omega pi0 pi0  (pip:1, pim:1, photons:6)
    (3) D0 -> eta pi+ pi-    (pip:2, pim:2, photons:2)
    (4) D0 -> eta pi0 pi0    (pip:1, pim:1, photons:6)
    Track/photon multiplicities differ across modes. The DSL signal_side
    uses minimal charged tracks (pip:1,pim:1) and photon range (2..6)
    with at_least:true to accommodate all modes. Mode-specific selection
    (exact track count, pi0 count) is applied in ROOT.")
  .note(:tag_modes, "Three anti-D0 tag modes: K+ pi-, K+ pi- pi0, K+ pi- pi- pi+.
    The K+ pi- pi+ pi- mode includes a K_S^0 veto:
    |M(pi+ pi-) - M(K_S^0)| > 30 MeV/c^2 to suppress D0 -> K_S^0 K+ pi- background.
    This veto is applied via DTagAlg/DTagTool internal selection.")
  .note(:omega_eta_reconstruction, "omega/eta candidates: pi+ pi- pi0 combinations
    with M(pi+ pi- pi0) < 0.9 GeV/c^2 retained. For D0->omega/eta pi+ pi-:
    four possible pi+ pi- pi0 combinations; for D0->omega/eta pi0 pi0:
    three combinations. omega signal region: 0.74 < M(3pi) < 0.82 GeV/c^2.
    eta signal region: 0.52 < M(3pi) < 0.57 GeV/c^2.
    Signal extraction via fits to M(pi+ pi- pi0) distribution.")
  .note(:ks0_veto, "K_S^0 veto applied to suppress CF backgrounds:
    0.475 < M(pi+ pi-) < 0.520 GeV/c^2 or 0.448 < M(pi0 pi0) < 0.548 GeV/c^2
    candidates rejected. Applied at ROOT level after event reconstruction.")
  .note(:two_dim_fit, "DT signal yield from 2D unbinned ML fit to
    M_BC^tag vs M_BC^sig. PDF includes signal + three background types:
    BKGI (one D correct, one incorrect), BKGII (q qbar continuum),
    BKGIII (neither D correct). Sideband subtraction for background
    with same final states but no omega/eta (BKGIV). Peaking background
    from D0 -> K_S^0 omega/eta (BKGV) estimated from known BFs.")
  .note(:pi0_reconstruction, "pi0 candidates: gamma-gamma pair with
    0.115 < M_gg < 0.150 GeV/c^2. At least one photon in EMC barrel.
    Kinematic fit constrains M_gg to pi0 nominal mass. For multi-pi0
    modes (D0->omega/eta pi0 pi0 with 3 pi0), each pi0 constrained
    separately before D candidate formation (not expressible in single
    TagAnalysis fit which constrains only one pi0 pair).")
  .note(:quantum_correlation, "At psi(3770), D0 D0bar produced in C=-1
    entangled state. Strong-phase correction factors c_f^i applied for
    flavor tags. B(D0->omega/eta pi pi) corrected: B_CP± = B_sig /
    [1 - c_f^i (2 f_CP+ - 1)]. f_CP+ uncertainty 7.3% for omega modes,
    0.8% for eta pi+ pi-. Applied in ROOT after BF calculation.")
  .note(:deltaE_requirement, "DT signal: deltaE within 3.0 (3.5) times
    resolution for D0 -> omega pi+ pi- (D0 -> omega pi0 pi0). If multiple
    combinations, select minimum |deltaE|. Mode-dependent deltaE cut
    applied in ROOT analysis stage.")
  .note(:systematic_uncertainties, "Additive: signal PDFs, fit bias,
    background PDF, BKGIV/BKGV contributions. Multiplicative: tracking 2%,
    PID 2%, pi0 reco 2-6%, deltaE requirement, K_S^0 veto, omega/eta
    signal region, MC generator 2-3%, ST yield 1.2%, strong-phase 7.3%
    (omega)/0.8% (eta), intermediate BFs. Total additive 4.3-9.3 events,
    multiplicative 5.3-10.0% depending on mode.")

alg_D0.apply

# ============================================================
# Process 2: D+ -> omega/eta pi+ pi0
# Tag side: D- via hadronic modes
# Signal side: D+ -> (omega/eta) + pi+ pi0
# ============================================================

alg_Dp = TagAnalysis.new("DpOmegaEtaPiPi0")
alg_Dp.set_header(["DpOmegaEtaPiPi0Alg/DpOmegaEtaPiPi0.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .with_decay_card(decay_card)

alg_Dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,       # D- -> K+ pi- pi-
          :DptoKPiPiPi0,    # D- -> K+ pi- pi- pi0
          :DptoKsPi,        # D- -> K_S^0 pi-
          :DptoKsPiPi0,     # D- -> K_S^0 pi- pi0
          :DptoKsPiPiPi,    # D- -> K_S^0 pi- pi- pi+
          :DptoKKPi         # D- -> K+ K- pi-
  t.charm -1                # tag D-
end

alg_Dp.signal_side do |s|
  s.charged(pip: 2, pim: 1, at_least: true)
  s.photons 4
  s.require_charge 1
end

alg_Dp.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_Dp
  .note(:signal_modes, "Two signal modes from remaining tracks:
    (1) D+ -> omega pi+ pi0  (pip:2, pim:1, photons:4 from 2 pi0)
    (2) D+ -> eta pi+ pi0    (same final state topology)
    The DSL signal_side declares pip:2,pim:1 with at_least:true
    and photons:4 minimum. Mode distinction via M(pi+ pi- pi0)
    mass window in ROOT (omega: 0.74-0.82, eta: 0.52-0.57 GeV/c^2).")
  .note(:tag_modes, "Six D- tag modes comprising ~28% of total D- decays:
    K+ pi- pi-, K+ pi- pi- pi0, K_S^0 pi-, K_S^0 pi- pi0,
    K_S^0 pi- pi- pi+, K+ K- pi-.")
  .note(:tag_side_vetoes, "Lambda veto for K_S^0 pi- pi0, K_S^0 pi- pi- pi+
    modes: M(pbar pi+) in [1.110, 1.120] rejected. K_S^0 veto for
    K_S^0 pi- pi- pi+ mode: M(pi+ pi-) in [0.480, 0.520] rejected.
    Sigma- veto for K_S^0 pi- pi0 mode: M(pbar pi0) in [1.170, 1.200]
    rejected. Applied via DTagAlg internal selection.")
  .note(:omega_eta_reconstruction, "For D+ -> omega/eta pi+ pi0:
    four possible pi+ pi- pi0 combinations. omega signal region:
    0.74 < M(3pi) < 0.82 GeV/c^2; eta: 0.52 < M(3pi) < 0.57 GeV/c^2.
    Signal extraction via fit to M(pi+ pi- pi0) with Crystal Ball
    signal + reversed ARGUS background.")
  .note(:ks0_veto, "Same K_S^0 veto as D0 analysis: 0.475 < M(pi+ pi-) < 0.520
    or 0.448 < M(pi0 pi0) < 0.548 GeV/c^2 candidates rejected.")
  .note(:two_dim_fit, "DT signal yield from 2D unbinned ML fit to
    M_BC^tag vs M_BC^sig. Same background model (BKGI/BKGII/BKGIII)
    with sideband subtraction for BKGIV and known-BF peaking background
    subtraction for BKGV. Projection plots used for signal extraction.")
  .note(:pi0_reconstruction, "Two pi0 candidates reconstructed per event
    (pi0 from omega/eta + bachelor pi0). Each pi0: 0.115 < M_gg < 0.150
    GeV/c^2. Kinematic fit constrains each M_gg to pi0 nominal mass.
    TagAnalysis fit constrains only one pi0 pair; second pi0 constrained
    via kalman_kinematic_fit at ROOT level or additional pre-fit step.")
  .note(:deltaE_requirement, "DT signal: deltaE within 3.5 times resolution
    for D+ -> omega/eta pi+ pi0. Best combination by minimum |deltaE|.
    Mode-dependent deltaE cut applied in ROOT.")
  .note(:systematic_uncertainties, "Additive from fit: signal PDFs, fit bias,
    background shapes, BKGIV/BKGV subtraction. Multiplicative: tracking 1.5%,
    PID 1.5%, pi0 reco 4%, deltaE 0.3%, K_S^0 veto 0.8%, omega/eta region
    0.2%, MC generator 3.5%, ST yield 0.4%, intermediate BFs.
    Total systematic 5.9-6.0% depending on mode.")

alg_Dp.apply

# ============================================================
# Execute
# ============================================================

alg_D0.execute_on([data_3773, incMC_3773, exMC])
alg_Dp.execute_on([data_3773, incMC_3773, exMC])