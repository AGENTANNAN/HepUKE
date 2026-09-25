# ============================================================
# Paper: 2006.02800v1
# Analysis of D0 -> K_S^0 K+ K- at psi(3770)
# Includes untagged (BF measurement) and flavor-tagged (Dalitz analysis)
# ============================================================

data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card: psi(3770) -> D0 anti-D0, signal D0 -> K_S^0 K+ K-
decay_card = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0  PHSP;
    Enddecay
    Decay D0
    1.0000 K_S0 K+ K-   PHSP;
    Enddecay
    Decay anti-D0
    1.0000 K+ pi-       PHSP;
    Enddecay
    Decay K_S0
    1.0000 pi+ pi-      PHSP;
    Enddecay
    End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_D0_KSKK"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# ============================================================
# Untagged analysis: reconstruct D0 -> K_S^0 K+ K-
# ============================================================

alg_untagged = Algorithm.new("D0KSKKUntagged")
alg_untagged.set_header(["D0KSKKUntaggedAlg/D0KSKKUntagged.h"])
              .set_constant({ "ECMS" => [:double, 3.773] })

sel_untagged = Selection.new
sel_untagged.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"       # K+ + pi+ from K_S^0
    nChrn ">=2"       # K- + pi- from K_S^0
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp ">=1"
    nkm ">=1"
    npip ">=1"
    npim ">=1"
  end
  # K_S^0 -> pi+ pi- (secondary vertex only; mass window: ROOT-level)
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # 1C kinematic fit: constrain (K_S^0 K+ K-) to D0 nominal mass; chi2 < 20
  .kinematic_fit([:K_S0, :kp, :km]) do
    nominal
    invariant_mass_of(:K_S0, :kp, :km).constrain_to_nominal_mass_of(:D0)
    chi2_cut 20
  end

alg_untagged
  .note(:ks0_track_requirements, "K_S^0 daughter tracks (pi+pi-):
    Vz < 20 cm, no Vr constraint, |cos(theta)| < 0.93. No PID applied to
    K_S^0 daughters. The DSL select_track block applies common Vz=10cm and
    Vr=1cm to all tracks; the per-track-type relaxation for K_S^0 pions
    cannot be expressed separately.")
  .note(:ks0_selection, "K_S^0 selection: secondary vertex chi2 < 100,
    flight distance / sigma(flight distance) > 2 for untagged (>0 for tagged).
    Invariant mass window m_ks in [0.487, 0.511] GeV/c^2 (or m_bc vs m_ks
    2D box of +/-4 sigma around peak). Mass window applied at ROOT level.")
  .note(:pid_relative, "Kaon candidates: P(K) > P(pi); pion candidates:
    P(pi) > P(K). The DSL prob_cut + identify with :against is the closest
    approximation to relative-probability comparison PID.")
  .note(:best_candidate_untagged, "If multiple D0 candidates, select the one
    with smallest |E_D0 - E_beam|. Applied at ROOT level by analyzing
    all combinations stored by the kinematic fit.")
  .note(:two_dim_fit, "Signal yield from 2D unbinned ML fit to m_bc vs m_ks.
    Signal: product of Crystal Ball in m_bc and Crystal Ball in m_ks.
    Three background components: (1) non-K_S^0: CB(m_bc) * poly1(m_ks);
    (2) combinatorial: [f*ARGUS+(1-f)*Gauss](m_bc) * [g*poly1+(1-g)*CB](m_ks);
    (3) additional non-K_S^0 from tag side (untagged only).
    Untagged yield: 11660+/-118; tagged yield: 1856+/-45 (96.37% purity).")
  .note(:tagged_sample, "Flavor-tagged sample: both D0 and anti-D0 reconstructed.
    Tag modes: K-pi+, K-pi+pi0, K-pi+pi+pi-pi0, K-pi+pi+pi-, K-pi+pi0pi0,
    K-pi+eta. Tag D0 reconstructed from scratch (not DTagAlg).
    Best combination: average m_bc closest to D0 nominal mass.")
  .note(:tag_mode_unavailable, "Flavor tag channels with pi0/eta require
    kalman_kinematic_fit for pi0/eta reconstruction before D0 formation.
    Also, the D0 mass constraint (1C) is a different fit topology than
    4C at CMS. The combined tag+signal reconstruction requires a separate
    algorithm instance per tag mode, similar to Rule T1.
    The tag D0 kinematic fit also constrains mass to D0 nominal (1C).")
  .note(:quantum_entanglement, "At psi(3770), D0-D0bar produced in C=-1
    entangled state. The Dalitz amplitude accounts for interference:
    |M|^2 approx |A_bar_tag * A_3K * lambda_tag - A_3K|^2.
    For self-conjugate K_S^0 K+ K-: A_bar_3K(m^2_KSK+, m^2_KSK-) =
    A_3K(m^2_KSK-, m^2_KSK+). Hadronic parameters from global fit:
    r_D^2(Kpi) = 0.344%, delta_D(Kpi) = 9.8 deg.")
  .note(:dalitz_resonances, "Dalitz model: a0(980)^0, a0(980)^+,
    phi(1020), a2(1320)^+, a2(1320)^-, a0(1450)^-. S-wave KK modeled
    by Flatte formula with g_KK = 3.77+/-0.24+/-0.35 GeV.
    Dalitz analysis uses ComPWA framework; not BOSS/DSL expressible.")
  .note(:systematic_uncertainties, "BF: tracking/PID ~1%/track, K_S^0 reco
    1.5%, 2D fit, signal MC model, resolution modeling. Total 3.6%.
    Dalitz: resonance parameterization, background model, fit bias,
    efficiency variations, quantum correlation parameters.")
  .with_decay_card(decay_card)
  .apply(sel_untagged)

alg_untagged.execute_on([data_3773, incMC_3773, exMC])