# =============================================================================
# BESIII: observation of Lambda_c+ -> Sigma- pi+ pi+ pi0
# arXiv:1705.11109v1
#
# 567 pb^-1 of e+e- collisions at sqrt(s) = 4.6 GeV, just above the
# Lambda_c+ anti-Lambda_c- threshold. Single-tag (ST) / double-tag (DT) method:
#   * the anti-Lambda_c- is reconstructed in eleven hadronic decay modes taken
#     from the pre-stored BOSS tag collection -> single tag (ST) sample;
#   * the recoiling Lambda_c+ is reconstructed as
#         Lambda_c+ -> Sigma- pi+ pi+ (pi0),   Sigma- -> n pi-,
#     the neutron being undetected and inferred from four-momentum conservation
#     -> double tag (DT) content.
# The absolute branching fraction follows from the DT/ST probability.
#
# DTagTool::findDTag rejects Lambda_c modes (decayMode() >= 1000), so the DT is
# not declared as two tag_side blocks: it is one Lambda_c tag side (ST) plus a
# signal_side carrying the recoiling Lambda_c+ final state and the missing
# neutron (ST + missing pattern).
# =============================================================================

### Dataset description ###
data_4600  = DatasetManager.real_data.find("703_4600")     # sqrt(s) = 4.599 GeV
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

### Decay cards — signal Lambda_c+ -> Sigma- pi+ pi+ (pi0), Sigma- -> n pi- ###
# The tag side (anti-Lambda_c-) is generated generically; the signal side is
# generated in phase space, matching the paper's phase-space signal samples.

decay_card_modes = {}

decay_card_modes[:SigmaPiPi] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-   PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  Sigma-  pi+  pi+            PHSP;
    Enddecay

    Decay Sigma-
    1.0000  n0  pi-                     PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                    PHSP;
    Enddecay

    End
DECAYCARD

decay_card_modes[:SigmaPiPiPi0] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-   PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  Sigma-  pi+  pi+  pi0       PHSP;
    Enddecay

    Decay Sigma-
    1.0000  n0  pi-                     PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                    PHSP;
    Enddecay

    End
DECAYCARD

exMC_modes = decay_card_modes.map do |key, card|
  DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "LcToSigmaPiPi#{key}_exclusive_mc"
    config.related_dataset = data_4600
    config.events          = 200_000
    config.decay_card      = card
    config.cross_section   = :default
  end
end

### Tag-side mode list — the eleven ST modes (DTagAlg channel names) ###
# Charge-conjugate modes are implied throughout (the tagged baryon may be
# Lambda_c- or Lambda_c+); the tag side is pinned to the charm -1 (anti-Lambda_c-)
# side, the recoiling baryon being the Lambda_c+.
tag_modes = [:LambdacPtoKsP,             # p K_S0
             :LambdacPtoKPiP,            # p K- pi+
             :LambdacPtoKPiPi0P,         # p K- pi+ pi0
             :LambdacPtoKsPi0P,          # p K_S0 pi0
             :LambdacPtoKsPiPiP,         # p K_S0 pi+ pi-
             :LambdacPtoLambdaPi,        # Lambda pi+
             :LambdacPtoLambdaPiPi0,     # Lambda pi+ pi0
             :LambdacPtoLambdaPiPiPi,    # Lambda pi+ pi- pi+
             :LambdacPtoPi0SIGMAPi0P,    # Sigma+ pi0            (Sigma+ -> p pi0)
             :LambdacPtoPiSIGMA0LambdaGam, # Sigma0 pi+          (Sigma0 -> Lambda gamma)
             :LambdacPtoPiPiSIGMAPi0P]   # Sigma+ pi+ pi-        (Sigma+ -> p pi0)

### Notes shared by the two signal channels ###
common_notes = [
  [:single_double_tag_method,
   "The anti-Lambda_c- is reconstructed (ST) in the eleven hadronic modes from the " \
   "pre-stored tag collection; the kinematic variable is the beam-energy-constrained " \
   "mass M_BC c^2 = sqrt(E_beam^2 - p^2 c^2) with E_beam the beam energy and p the " \
   "anti-Lambda_c- momentum in the e+e- c.m. system. For each mode only the candidate " \
   "with the minimum |Delta E| is kept, Delta E = E_beam - E(anti-Lambda_c-), and the " \
   "mode-dependent Delta E windows of Table 1 (about +/-3 sigma, ranging from " \
   "(-0.025,0.028) to (-0.049,0.054) GeV) are applied on the stored Delta E. The ST " \
   "yield is obtained from a fit to M_BC in the signal region " \
   "2.280 < M_BC < 2.296 GeV/c^2; the total ST yield is N_tot = 14415 +/- 159."],
  [:tag_side_selection,
   "The tag-side tracks carry the standard BESIII quality requirements: " \
   "|cos(theta)| < 0.93 and a distance of closest approach to the interaction " \
   "point below 10 cm along the beam axis and below 1 cm in the perpendicular " \
   "plane (pions from K_S0 decays are exempt). Protons, kaons and pions are " \
   "separated with combined dE/dx (MDC) and TOF likelihoods. Intermediate " \
   "candidates: K_S0 -> pi+ pi-, Lambda -> p pi-, Sigma0 -> gamma Lambda with " \
   "Lambda -> p pi-, Sigma+ -> p pi0, pi0 -> gamma gamma, using the criteria of " \
   "Ref. [14]."],
  [:signal_side_selection,
   "The three charged tracks recoiling against the tag must satisfy " \
   "|cos(theta)| < 0.93; the two pi+ from the Lambda_c+ must have distances of " \
   "closest approach to the interaction point within +/-10 cm along the beam " \
   "direction and within 1 cm in the perpendicular plane, while the pi- from " \
   "Sigma- is exempt from this requirement. The three charged pions are required " \
   "to satisfy L(pi) > L(K). Photons are isolated EMC clusters with " \
   "|cos(theta)| <= 0.80 (barrel) or 0.86 <= |cos(theta)| <= 0.92 (end cap), " \
   "energy above 25 (50) MeV in the barrel (end cap), angle to the nearest " \
   "charged track above 10 degrees, and EMC time within (0,700) ns of the event " \
   "start time. pi0 candidates from photon pairs with " \
   "0.110 < M(gamma gamma) < 0.155 GeV/c^2 are mass-constrained to the nominal " \
   "pi0 mass."],
  [:missing_neutron,
   "The neutron from Sigma- -> n pi- is not reconstructed; its kinematic " \
   "properties are deduced from four-momentum conservation (declared here as the " \
   "missing particle of the fit). The reconstructed neutron mass " \
   "M_n = sqrt[(E_beam - E_pipipim(pi0))^2 - |p_Lambda_c+ - p_pipipim(pi0)|^2] and " \
   "the Sigma- mass M_npi- replace the undetected neutron; the expected " \
   "Lambda_c+ momentum is p_Lambda_c+ = -p_hat_tag sqrt(E_beam^2 - m_Lambda_c+^2). " \
   "The signal yield is extracted from an unbinned maximum-likelihood fit to " \
   "M_npi- - M_n: the signal shape is a non-parametric function from signal MC " \
   "convoluted with a Gaussian (resolution data vs MC), the background a " \
   "second-order polynomial. Found yields: N_SigmaPiPi = 161 +/- 15 and " \
   "N_SigmaPiPiPi0 = 88 +/- 14."],
  [:background_veto,
   "Backgrounds from non-Lambda_c+ decays are estimated from the M_BC sideband " \
   "(2.252, 2.272) GeV/c^2 in data and found to be negligible; the inclusive MC " \
   "study shows no peaking background for either channel."],
  [:signal_modelling,
   "The Lambda_c+ -> Sigma- pi+ pi+ (pi0) decay in the MC is generated in phase " \
   "space and reweighted to approximate the observed kinematic distributions in " \
   "data; the reweighting factors are varied within their statistical " \
   "uncertainties for the systematic error on the efficiency."],
  [:charge_conjugate_modes,
   "Throughout the analysis charge-conjugate modes are implicitly assumed: the " \
   "tagged baryon may be Lambda_c- or Lambda_c+, the recoiling baryon " \
   "correspondingly Lambda_c+ or Lambda_c-. Only the anti-Lambda_c- tag (charm " \
   "pinned to -1) with the Lambda_c+ recoil is declared; the charge-conjugate " \
   "selection differs only by the sign of every charge."],
]

### Event selection (BOSS) — TagAnalysis, one algorithm per signal channel ###

# ---------------------------------------------------------------------------
# Channel 1: Lambda_c+ -> Sigma- pi+ pi+   (Sigma- -> n pi-, neutron missing)
# ---------------------------------------------------------------------------
alg_sigpipi = TagAnalysis.new("LcToSigmaPiPi")
alg_sigpipi.set_header(["LcToSigmaPiPiAlg/LcToSigmaPiPi.h"])
           .set_constant({ "ECMS" => [:double, 4.600] })
           .with_decay_card(decay_card_modes[:SigmaPiPi])

alg_sigpipi.tag_side(:Lambdac) do |t|     # tagged anti-Lambda_c- (charm pinned)
  t.modes(*tag_modes)
  t.charm -1
end

alg_sigpipi.signal_side do |s|            # recoiling Lambda_c+ -> Sigma- pi+ pi+
  s.charged(pip: 2, pim: 1)              # two pi+ and the pi- from Sigma- -> n pi-
  s.require_charge 1                     # +2 - 1 = +1
  s.photons 0
  s.missing :n0                          # neutron (massive form), inferred by 4-momentum conservation
end

alg_sigpipi.fit do |f|
  f.constrain_four_momentum              # 4C: tag + pi+ pi+ pi- + n = ecms_lab
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_sigpipi.note(key, text) }
alg_sigpipi.apply
alg_sigpipi.execute_on([data_4600, incMC_4600] + exMC_modes)

# ---------------------------------------------------------------------------
# Channel 2: Lambda_c+ -> Sigma- pi+ pi+ pi0   (Sigma- -> n pi-, neutron missing)
# ---------------------------------------------------------------------------
alg_sigpipipi0 = TagAnalysis.new("LcToSigmaPiPiPi0")
alg_sigpipipi0.set_header(["LcToSigmaPiPiPi0Alg/LcToSigmaPiPiPi0.h"])
              .set_constant({ "ECMS" => [:double, 4.600] })
              .with_decay_card(decay_card_modes[:SigmaPiPiPi0])

alg_sigpipipi0.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_sigpipipi0.signal_side do |s|         # recoiling Lambda_c+ -> Sigma- pi+ pi+ pi0
  s.charged(pip: 2, pim: 1)              # two pi+ and the pi- from Sigma- -> n pi-
  s.require_charge 1
  s.photons 2                            # the two photons of the pi0
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :n0                          # neutron (massive form), inferred by 4-momentum conservation
end

alg_sigpipipi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_sigpipipi0.note(key, text) }
alg_sigpipipi0.apply
alg_sigpipipi0.execute_on([data_4600, incMC_4600] + exMC_modes)
