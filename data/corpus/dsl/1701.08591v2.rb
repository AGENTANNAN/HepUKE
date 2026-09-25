# =============================================================================
# 1701.08591v2
# Amplitude analysis of D0 -> K- pi+ pi+ pi-  (BESIII)
#
# Data: 2.93 fb^-1 collected at sqrt(s) = 3.773 GeV (psi(3770) resonance).
#
# The analysis reconstructs the D0 D0bar pair with the double-tag technique:
# one D0 is tagged through the hadronic mode D0 -> K+ pi- and the signal mode
# D0 -> K- pi+ pi+ pi- is reconstructed from the tracks not used by the tag.
# Because the tag is taken from the pre-stored EvtRecDTag collection, the spec
# is a TagAnalysis (DTagAlg / DTagTool), not an Algorithm + Selection.
#
# BOSS-side scope: decay cards, exclusive MC, tag declaration, signal-side
# declaration and the 5C kinematic fit. The amplitude model, the unbinned
# maximum-likelihood fit, the fit fractions and the phases are ROOT-level and
# are not part of this spec.
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets
### ---------------------------------------------------------------------------
psipp_data  = DatasetManager.real_data.find("712_3773")     # psi(3770), 2.93 fb^-1
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")   # generic inclusive MC

### ---------------------------------------------------------------------------
### Decay cards
### ---------------------------------------------------------------------------
# Signal: psi(3770) -> D0 D0bar with the tagged anti-D0 -> K+ pi- and the
# signal D0 -> K- pi+ pi+ pi-. Phase space is used for the four-body signal
# decay; the nominal amplitude model is applied when generating the SIGNAL MC
# sample used for the efficiency (see the notes).
decay_card_signal = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0000 K- pi+ pi+ pi- PHSP;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Dedicated background sample: psi(3770) -> D0 D0bar with the tag
# anti-D0 -> K+ pi- and the peaking background D0 -> K_S0 K- pi+, which has
# the same final state as the signal.
decay_card_kskpi_bkg = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0000 K_S0 K- pi+ PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi- PHSP;
  Enddecay

  End
DECAYCARD

### ---------------------------------------------------------------------------
### Exclusive MC
### ---------------------------------------------------------------------------
# PHSP signal sample — used for the MC integration in the amplitude fit.
exMC_signal_phsp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "D0toKPiPiPi_phsp_mc"
  config.related_dataset = psipp_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Dedicated K_S0 K- pi+ background MC sample.
exMC_kskpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "D0toKsKPi_bkg_mc"
  config.related_dataset = psipp_data
  config.events          = 500_000
  config.decay_card      = decay_card_kskpi_bkg
  config.cross_section   = :default
end

### ---------------------------------------------------------------------------
### Tag analysis (BOSS)
### ---------------------------------------------------------------------------
alg_name = "D0ToKPiPiPiDTag"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })
   .with_decay_card(decay_card_signal)

# Tag side: the anti-D0 -> K+ pi- tag (charm -1). The tag mBC / DeltaE are
# normally stored unconditionally; the paper applies explicit windows, which
# are declared here as the sanctioned opt-in tag-side windows.
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi
  t.charm -1
  t.window :mBC,    min: 1.8575, max: 1.8775   # 1.8575 < M_BC < 1.8775 GeV/c^2
  t.window :deltaE, min: -0.03,  max: 0.03     # -0.03 < DeltaE < 0.03 GeV (K+ pi- tag)
end

# Signal side: the four tracks recoiling against the tag, K- pi+ pi+ pi-.
# They have distinct momenta from the tag daughters, so the multiset is
# unambiguous; the total charge -1 + 2 - 1 = 0 is checked against the tag.
alg.signal_side do |s|
  s.charged(km: 1, pip: 2, pim: 1)
  s.require_charge 0
end

# 5C kinematic fit: the total four-momentum of all final-state particles is
# constrained to the initial e+e- four-momentum and the invariant mass of the
# signal-side K- pi+ pi+ pi- is constrained to the PDG D0 mass.
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:km, :pip, :pip, :pim).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200     # loose BOSS-pass cut; the paper's chi2 < 40 is applied in ROOT
end

# BOSS-side procedures without a DSL construct.
alg
  .note(:kinematic_fit_chi2, "the paper requires the chi2 of the 5C kinematic fit to be less than 40; BOSS applies only the loose default chi2_cut 200 and the chi2 < 40 requirement is applied in the ROOT analysis")
  .note(:track_selection, "good charged tracks must have a point of closest approach to the IP within 10 cm along the beam axis and within 1 cm in the plane perpendicular to the beam, with |cos(theta)| < 0.93; the tag side inherits the DTagAlg/DtagTool track quality and the signal-side tracks are the tag's leftover tracks")
  .note(:pid_correction_method, "PID combines the MDC dE/dx and the TOF information; charged kaons require P(K) > P(pi) and pions require P(pi) > P(K), and tracks without PID information are rejected; the signal-side kaon/pion classification is driven by the SimplePIDSvc recipes inside the tag fit")
  .note(:tag_side_selection, "the D0 D0bar pair must have opposite charm and no tracks in common; this is enforced by the DTagTool pairing and the signal-side leftover-track requirement")
  .note(:vertex_fit, "a vertex fit with the hypothesis that all tracks originate from the IP is performed on the D0 D0bar pair and the fit chi2 is required to be less than 200; this pre-fit has no DSL counterpart and is applied in ROOT")
  .note(:delta_e_window, "the signal-side DeltaE window -0.033 < DeltaE < 0.033 GeV for the K- pi+ pi+ pi- final state (and the shared M_BC window 1.8575 < M_BC < 1.8775 GeV/c^2) are applied to the signal-side observables in ROOT; only the tag-side windows are expressible in the tag DSL")
  .note(:ks0_veto, "for any pi+ pi- pair on the signal side whose invariant mass satisfies |m(pi+ pi-) - m_K_S0| < 0.03 GeV/c^2, a vertex-constrained fit is performed and the event is rejected if the decay-length significance exceeds 2 sigma; this removes about 80% of the D0 -> K_S0 K- pi+ background while retaining about 99% of the signal; the pi+ pi- combination scan and the decay-length significance are not expressible in the DSL")
  .note(:background_estimation, "the peaking background D0 -> K_S0 K- pi+ is described by the dedicated K_S0 K pi MC sample and its normalisation (96.8 +/- 14.5 events) is fixed in the amplitude fit; it is assigned a negative weight in the likelihood")
  .note(:efficiency_curve, "the detection efficiency entering the normalisation integral depends on the final four-momenta; the MC integration weights each generated event by the ratio of data to MC tracking and PID efficiencies, gamma_eps(p) = prod_j eps_j,data(p_j)/eps_j,MC(p_j)")
  .note(:amplitude_model, "the SIGNAL MC sample used to estimate the efficiency and the goodness of fit is generated according to the amplitude model obtained in this analysis (23 amplitudes in seven components, with the D0 -> K- a1+(1260), a1+(1260)[S] -> rho0 pi+ amplitude fixed to unit magnitude and zero phase); the amplitude fit itself is a ROOT-level procedure")
  .note(:bose_symmetry, "the two pi+ mesons of the K- pi+ pi+ pi- final state impose an explicit Bose symmetrisation of the total amplitude; this is part of the ROOT amplitude model")

alg.apply

### ---------------------------------------------------------------------------
### Execution
### ---------------------------------------------------------------------------
root_files = alg.execute_on([psipp_data, psipp_incMC, exMC_signal_phsp, exMC_kskpi])
