# =============================================================================
# Λc+ → Σ−π+π+ (π0)  at √s = 4.600 GeV  —  single-tag (DTag) analysis
#
# The anti-Λc− is taken from the pre-stored DTag candidates (tag side, eleven
# hadronic modes); the recoiling Λc+ is the signal side: three pions + an
# optional π0 → γγ, with the neutron from Σ− → nπ− undetected and inferred from
# four-momentum conservation.  No Selection / select_track / select_photon /
# pid is emitted — the tag already carries its own selected, PID'd tracks and
# showers, and the signal side comes from the tag's unused tracks/showers.
# =============================================================================

### Datasets ###
data_4600  = DatasetManager.real_data.find("703_4600")      # 4.600 GeV real data (567 pb^-1)
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")   # 4.600 GeV inclusive MC

### Signal decay cards (EvtGen); the tag side decays generically and is left
### without a decay block, so only the signal Λc+ is forced. ###
decay_card_sigma_pipi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c-   PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 Sigma- pi+ pi+   PHSP;
    Enddecay

    Decay Sigma-
    1.0000 n0 pi-   PHSP;
    Enddecay

    End
DECAYCARD

decay_card_sigma_pipi_pi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c-   PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 Sigma- pi+ pi+ pi0   PHSP;
    Enddecay

    Decay Sigma-
    1.0000 n0 pi-   PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma   PHSP;
    Enddecay

    End
DECAYCARD

### 200k-event exclusive phase-space MC for each signal mode ###
exMC_sigma_pipi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4600_lambdac_to_sigmapippi"
  config.related_dataset = data_4600          # associated real dataset
  config.events          = 200_000
  config.decay_card      = decay_card_sigma_pipi
  config.cross_section   = :default
end

exMC_sigma_pipi_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4600_lambdac_to_sigmapipppi0"
  config.related_dataset = data_4600
  config.events          = 200_000
  config.decay_card      = decay_card_sigma_pipi_pi0
  config.cross_section   = :default
end

### The eleven hadronic tag modes of the anti-Λc− ###
TAG_MODES = [
  :LambdacPtoKsP,           # p K_S0
  :LambdacPtoKPiP,          # p K− π+
  :LambdacPtoKPiPPi0,       # p K− π+ π0
  :LambdacPtoKsPPi0,        # p K_S0 π0
  :LambdacPtoKsPPiPi,       # p K_S0 π+ π−
  :LambdacPtoLambdaPi,      # Λ π+
  :LambdacPtoLambdaPiPi0,   # Λ π+ π0
  :LambdacPtoLambdaPiPiPi,  # Λ π+ π− π+
  :LambdacPtoSigmaPi0,      # Σ+ π0
  :LambdacPtoSigma0Pi,      # Σ0 π+
  :LambdacPtoSigmaPiPi      # Σ+ π+ π−
].freeze

# =============================================================================
# Channel I :  Λc+ → Σ−π+π+     (no photon in the final state)
# =============================================================================
alg_3pi = TagAnalysis.new("LcTagSigmaPiPiPi")
alg_3pi.set_header(["LcTagSigmaPiPiPiAlg/LcTagSigmaPiPiPi.h"])
       .set_constant({"ECMS" => [:double, 4.600]})
       .with_decay_card(decay_card_sigma_pipi)

# Tag side: the anti-Λc−, reconstructed in the eleven hadronic modes.
alg_3pi.tag_side(:Lambdac) do |t|
  t.modes(*TAG_MODES)
  t.charm(-1)                       # pin the tagged anti-Λc−
  t.window :deltaE, abs: 0.02       # opt-in tag-side ΔE window (loose; per-mode ±3σ in ROOT)
end

# Signal side: the recoiling Λc+ → Σ−π+π+; the neutron of Σ− → nπ− is missing.
alg_3pi.signal_side do |s|
  s.charged(pip: 2, pim: 1)         # exactly two π+ and one π− recoil against the tag
  s.require_charge 1                # net charge of the signal side is +1
  s.missing :n                      # undetected neutron (massive), inferred from 4-momentum
end

# 4C kinematic fit: tag + π+π+π− + n to the beam energy, tag Λc+ at nominal mass.
alg_3pi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

# Inexpressible tag-side / signal-side procedures (BOSS side).
alg_3pi
  .note(:tag_track_selection, "tag-side charged tracks are required to satisfy |cosθ| < 0.93, |Vz| < 10 cm and Vr < 1 cm, with p/K/π separation from the combined dE/dx and TOF likelihoods; the K_S0 → π+π−, Λ → pπ−, Σ0 → γΛ, Σ+ → pπ0 and π0 → γγ candidates are built inside the DTagAlg tag reconstruction rather than re-implemented here")
  .note(:tag_candidate_selection, "for each tag mode the candidate with the minimum |ΔE| is retained, and the mode-dependent (≈±3σ) ΔE windows are applied per mode; the ST yield is extracted from a fit to M_BC in 2.280–2.296 GeV/c², both handled at ROOT level")
  .note(:signal_track_cuts, "the two signal-side π+ are required to satisfy |Vz| < 10 cm and Vr < 1 cm, while the π− track is exempt from the vertex requirement — the signal-side declaration cannot express this charge-dependent asymmetry")
  .note(:signal_pion_pid, "all three signal-side pions are required to satisfy L(π) > L(K) from the combined dE/dx and TOF likelihoods")
  .note(:signal_photon_veto, "no photon is allowed on the signal side of Λc+ → Σ−π+π+; a zero-photon veto on the recoil is not expressible in the signal-side declaration")

alg_3pi.apply                          # no Selection argument for a tag analysis
root_files_3pi = alg_3pi.execute_on([data_4600, incMC_4600, exMC_sigma_pipi])

# =============================================================================
# Channel II : Λc+ → Σ−π+π+π0, π0 → γγ   (two photons in the final state)
# =============================================================================
alg_3pi_pi0 = TagAnalysis.new("LcTagSigmaPiPiPiPi0")
alg_3pi_pi0.set_header(["LcTagSigmaPiPiPiPi0Alg/LcTagSigmaPiPiPiPi0.h"])
           .set_constant({"ECMS" => [:double, 4.600]})
           .with_decay_card(decay_card_sigma_pipi_pi0)

# Tag side: identical anti-Λc− tag.
alg_3pi_pi0.tag_side(:Lambdac) do |t|
  t.modes(*TAG_MODES)
  t.charm(-1)
  t.window :deltaE, abs: 0.02
end

# Signal side:  π+π+π− recoil against the tag, two photons from π0 → γγ, missing n.
alg_3pi_pi0.signal_side do |s|
  s.charged(pip: 2, pim: 1)
  s.require_charge 1
  s.photons 2                       # two signal photons consumed by the fit (π0 → γγ)
  s.min_photon_angle 10.0           # >10° from the nearest charged track
  s.min_photon_energy 0.025         # E > 25 MeV (barrel floor)
  s.missing :n                      # undetected neutron from Σ− → nπ−
end

# 4C kinematic fit + tag Λc+ mass constraint + γγ mass constraint to the nominal π0.
alg_3pi_pi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.invariant_mass_of(:gamma, :gamma).between(0.110, 0.155)         # π0 mass window
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_3pi_pi0
  .note(:tag_track_selection, "tag-side charged tracks are required to satisfy |cosθ| < 0.93, |Vz| < 10 cm and Vr < 1 cm, with p/K/π separation from the combined dE/dx and TOF likelihoods; the K_S0 → π+π−, Λ → pπ−, Σ0 → γΛ, Σ+ → pπ0 and π0 → γγ candidates are built inside the DTagAlg tag reconstruction rather than re-implemented here")
  .note(:tag_candidate_selection, "for each tag mode the candidate with the minimum |ΔE| is retained, and the mode-dependent (≈±3σ) ΔE windows are applied per mode; the ST yield is extracted from a fit to M_BC in 2.280–2.296 GeV/c², both handled at ROOT level")
  .note(:signal_track_cuts, "the two signal-side π+ are required to satisfy |Vz| < 10 cm and Vr < 1 cm, while the π− track is exempt from the vertex requirement")
  .note(:signal_pion_pid, "all three signal-side pions are required to satisfy L(π) > L(K) from the combined dE/dx and TOF likelihoods")
  .note(:signal_photon_selection, "signal-side photons must be isolated EMC clusters with |cosθ| ≤ 0.80 (barrel) or 0.86–0.92 (endcap) and EMC time in (0,700) ns, with E > 25 MeV in the barrel and E > 50 MeV in the endcap; only the 10° isolation and the 25 MeV floor are expressible in the signal-side declaration")

alg_3pi_pi0.apply
root_files_3pi_pi0 = alg_3pi_pi0.execute_on([data_4600, incMC_4600, exMC_sigma_pipi_pi0])