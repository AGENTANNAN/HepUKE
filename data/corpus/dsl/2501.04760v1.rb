# Paper 2501.04760v1 — Search for D^+ → e^+ ν_e
# with 20.3 fb⁻¹ at √s = 3.773 GeV (ψ(3770))
# Double-tag method
# Tag side: D^- → K^+π^-π^-, K_S^0π^-, K^+π^-π^-π^0, K_S^0π^-π^0, K_S^0π^+π^-π^-, K^+K^-π^-
# Signal side: D^+ → e^+ ν_e (semileptonic, missing neutrino)

### Dataset preparation ###
psipp_data = DatasetManager.real_data.find("712_3773")
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")

decay_card = <<~DECAYCARD
  Decay D+
  1.0000 e+ nu_e PHSP;
  Enddecay
  End
DECAYCARD

# Signal exclusive MC — D^+ → e^+ ν_e
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Dp_enu"
  config.related_dataset = psipp_data
  config.events = 500_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (TagAnalysis — ST + missing neutrino) ###
alg = TagAnalysis.new("DpToENu")

alg.set_header(["DpToENuAlg/DpToENu.h"])
    .set_constant("ECMS" => [:double, 3.773])

alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

alg.signal_side do |s|
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"D+")
  f.chi2_cut 50
end

alg.note(:tag_selection, "ST D^- tagged via 6 hadronic modes: K^+π^-π^- (DptoKPiPi), K_S^0π^- (DptoKsPi), K^+π^-π^-π^0 (DptoKPiPiPi0), K_S^0π^-π^0 (DptoKsPiPi0), K_S^0π^+π^-π^- (DptoKsPiPiPi), K^+K^-π^- (DptoKKPi)")
  .note(:st_selection, "ST D^- candidates selected with mode-dependent ΔE requirements (Table 1: [-25,24] to [-62,49] MeV) and M_BC ∈ (1.863, 1.877) GeV/c². Multiple candidates: keep the one with minimum |ΔE|. ST yields from fits to M_BC distributions (ARGUS background, double-Gaussian ⊗ MC signal shape)")
  .note(:track_selection, "Charged tracks: |cos θ| < 0.93, V_r < 1 cm, |V_z| < 10 cm (except K_S^0 daughters: |V_z| < 20 cm). PID: L(K) > L(π) for kaons; L(π) > L(K) for pions")
  .note(:Ks_reconstruction, "K_S^0 → π^+π^-: secondary vertex fit with χ² < 100; decay length L/σ_L > 2; M(π^+π^-) ∈ (0.487, 0.511) GeV/c². Fitted K_S^0 four-vectors used in kinematic calculations")
  .note(:photon_selection, "EMC showers: E > 25 MeV (barrel, |cos θ| < 0.80), E > 50 MeV (endcap, 0.86 < |cos θ| < 0.92). Shower-track angle > 10°. EMC time ∈ [0, 700] ns")
  .note(:pi0_reconstruction, "π^0 → γγ: M(γγ) ∈ (0.115, 0.150) GeV/c²; 1C mass-constrained fit (→ π^0 nominal mass) with χ² < 50. Fitted π^0 four-momentum used")
  .note(:e_pid, "Positron PID: L(e) > 0.001, L(e)/(L(e)+L(K)+L(π)) > 0.8. E/p > 0.8 required to suppress misidentified hadrons. Bremsstrahlung/FSR recovery: photons within 5° cone of e^+ direction added to e^+ four-momentum")
  .note(:kinematic_fit, "Kinematic fit constraining total 4-momentum to initial state + tag D^- and signal D^+ invariant masses to known D^± mass. Neutrino four-momentum determined by fit; χ² < 50 required")
  .note(:extra_particle_veto, "Suppress backgrounds: N_π0^extra = 0 (no extra π^0), N_char^extra = 0 (no extra good charged tracks), E_max,γ^extra < 0.2 GeV (max extra photon energy). Optimized by maximizing ε/(1.5+√B)")
  .note(:signal_extraction, "Signal yield N_DT = 0.3^(+2.9)_(-3.4) (stat) from fit to M_miss² = (E_beam − E_e+)^2 − (−p_D− − p_e+)². Signal shape from signal MC; background shape from inclusive MC smoothed with RooKeysPDF")
  .note(:upper_limit, "No significant signal observed. Upper limit B(D^+ → e^+ν_e) < 9.7×10⁻⁷ at 90% CL (Bayesian approach with smeared likelihood incorporating multiplicative and additive systematics)")
  .with_decay_card(decay_card)
  .apply

root_files = alg.execute_on([psipp_data, psipp_incMC, exMC_signal])