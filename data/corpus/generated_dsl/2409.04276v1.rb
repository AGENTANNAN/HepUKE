### Dataset preparation ###
# ψ(3770) = 3.773 GeV  ->  BOSS sample name "712_3773"
data_3773  = DatasetManager.real_data.find("712_3773")        # real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")     # corresponding inclusive MC

# Signal decay card (EvtGen): ψ(3770) → D0 D0bar with
#   signal  D0  → π⁻ π⁰ e⁺ ν_e  (semileptonic)
#   tag     D0bar → K⁺ π⁻        (one of the hadronic tag modes)
#   π⁰ → γ γ
decay_card_signal = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0000 pi- pi0 e+ nu_e PHSP;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi- PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC: 100k events of ψ(3770) → D Dbar (D0 semileptonic / D0bar → K⁺π⁻)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0_semilep_pip0enu"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection — double tag (TagAnalysis) ###
# One tag_side (hadronic D) + signal_side with a missing neutrino = DT pattern.
alg = TagAnalysis.new("D0SemilepDTag")
alg.set_header(["D0SemilepDTagAlg/D0SemilepDTag.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(decay_card_signal)

# --- Tag side: the other D reconstructed in the listed hadronic modes ---
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKsPiPi, :D0toKPiPi0Pi0, :D0toKPiPiPiPi0
end

# --- Signal side: D0 → π⁻ π⁰ e⁺ ν_e  (one π⁻, one e⁺, two photons, missing ν_e) ---
alg.signal_side do |s|
  s.charged(pim: 1, ep: 1)      # exactly one π⁻ and one e⁺
  s.photons 2                   # π⁰ → γγ
  s.min_photon_angle  10.0      # ≥ 10° away from any charged track
  s.min_photon_energy 0.025     # barrel energy floor (endcap 50 MeV captured in note)
  s.require_charge 0            # net charge zero
  s.missing :nu_e               # massless missing neutrino
end

# --- Kinematic fit: 4-momentum conservation + γγ mass constrained to nominal π⁰ ---
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)
  f.chi2_cut 200
end

# --- BOSS-side procedures that cannot be expressed in the tag DSL ---
alg.note(:tag_reconstruction_cuts, "tag-side D taken from the pre-stored EvtRecDTag
   collection (DTagAlg); charged tracks not originating from a K_S0 must satisfy
   |Vz|<10 cm and |Vxy|<1 cm, and K_S0 candidates require a secondary-vertex fit
   chi2<100, decay length >2 sigma and |M(pi+pi-)-m_K_S0|<12 MeV")
   .note(:signal_pid, "signal-side positron identified by L_e/(L_e+L_pi+L_K)>0.8,
   L_e>0.001 and E/(cp)>0.8; pi/K separated with dE/dx and TOF (L_pi>L_K for pions,
   L_K>L_pi for kaons). These likelihood cuts are not tunable through the signal_side
   ep key (fixed SimplePIDSvc recipe)")
   .note(:photon_selection, "signal photons require E>25 MeV in the barrel
   (|cos theta|<0.80) and E>50 MeV in the endcap (0.86<|cos theta|<0.92), at least
   10 degrees from any charged track, with EMC time within [0,700] ns")
   .note(:preliminary_1c_fit, "a preliminary 1C kinematic fit constraining the two
   photons to the nominal pi0 mass is required to give chi2<50 before the main 4C fit;
   the best pi0 candidate is the gamma-gamma pair whose invariant mass is closest to
   the PDG pi0 mass")
   .note(:background_veto, "background vetoes: E_max(extra gamma)<0.25 GeV, no extra
   charged tracks and no pi0 built from unused photons, M(pi- pi0 e+)<1.70 GeV,
   |M(pi- pi0)-m_K|<70 MeV for the BF measurement, |U_miss|<0.03 GeV for the amplitude
   analysis, and rejection of pi0 candidates with both photons in the EMC endcap")

alg.apply
alg.execute_on([data_3773, incMC_3773, exMC_signal])