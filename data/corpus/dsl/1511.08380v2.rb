# =============================================================================
# BESIII: first measurement of the absolute hadronic branching fractions of the
# Lambda_c+ baryon at the Lambda_c+ Lambda_c- production threshold
# arXiv:1511.08380v2
#
# 567 pb^-1 of e+e- collisions at sqrt(s) = 4.599 GeV, at the Lambda_c+ Lambda_c-
# threshold where no additional hadrons accompany the pair. Twelve
# Cabibbo-favored Lambda_c+ hadronic decay modes are analyzed with the
# double-tag technique:
#   N_j^ST  = N_LLbar * B_j * eps_j
#   N_ij^DT = N_LLbar * B_i * B_j * eps_ij
#   B_i     = (N_ij^DT / N_j^ST) * (eps_j / eps_ij)
# One Lambda_c baryon is reconstructed as the tag (single tag, ST) from BOSS's
# pre-stored tag collection (DTagAlg/DTagTool), and the recoiling baryon is
# reconstructed in one of the twelve modes (double tag, DT), so the whole
# Lambda_c+ Lambda_c- event is accounted for.
#
# The tag side is declared once, with all twelve modes; DTagTool pairs the two
# tag sides but rejects Lambda_c modes in findDTag (decayMode() >= 1000), so a
# Lambda_c double tag cannot be declared as two tag_side blocks. The DT content
# is therefore expressed as: one Lambda_c tag side (all twelve modes, charm
# pinned to the Lambda_c- side) + a fully reconstructed signal side carrying
# the lambda_c+ final state of mode i. Because the twelve signal modes have
# different final states and kinematic-fit constraints (Rule T1), each mode gets
# its own TagAnalysis algorithm. Charge-conjugate modes are implicit in the
# paper; they are covered by the same selection with all charges reversed (see
# the charge_conjugate_modes note).
# =============================================================================

### Dataset description ###
# 4.599 GeV threshold point (nominal 4600 MeV sample of the 703 round-07 scan).
data_4599  = DatasetManager.real_data.find("703_4600")
incMC_4599 = DatasetManager.inclusive_mc.find("703_4600")

### Decay cards for the exclusive signal MC ###
# One card per signal mode: Lambda_c+ decays to mode i, the opposite
# Lambda_c- decays generically ("anything"), as in the paper's signal samples.

signal_cards = {}

signal_cards["LcToPKsDT"] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-            PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  p+  K_S0                             PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                             PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                             PHSP;
    Enddecay

    End
DECAYCARD

signal_cards["LcToPKPiDT"] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-            PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  p+  K-  pi+                          PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                             PHSP;
    Enddecay

    End
DECAYCARD

signal_cards["LcToPKsPi0DT"] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-            PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  p+  K_S0  pi0                        PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                             PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                             PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                         PHSP;
    Enddecay

    End
DECAYCARD

signal_cards["LcToPKsPiPiDT"] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-            PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  p+  K_S0  pi+  pi-                   PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                             PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-                             PHSP;
    Enddecay

    End
DECAYCARD

signal_cards["LcToPKPiPi0DT"] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-            PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  p+  K-  pi+  pi0                     PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                             PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                         PHSP;
    Enddecay

    End
DECAYCARD

signal_cards["LcToLambdaPiDT"] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-            PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  Lambda0  pi+                         PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                             PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                              PHSP;
    Enddecay

    End
DECAYCARD

signal_cards["LcToLambdaPiPi0DT"] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-            PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  Lambda0  pi+  pi0                    PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                             PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                              PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                         PHSP;
    Enddecay

    End
DECAYCARD

signal_cards["LcToLambdaPiPiPiDT"] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-            PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  Lambda0  pi+  pi-  pi+               PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                             PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                              PHSP;
    Enddecay

    End
DECAYCARD

signal_cards["LcToSigma0PiDT"] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-            PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  Sigma0  pi+                          PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                             PHSP;
    Enddecay

    Decay Sigma0
    1.0000  gamma  Lambda0                       PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+  pi-                              PHSP;
    Enddecay

    End
DECAYCARD

signal_cards["LcToSigmaPPi0DT"] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-            PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  Sigma+  pi0                          PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                             PHSP;
    Enddecay

    Decay Sigma+
    1.0000  p+  pi0                              PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                         PHSP;
    Enddecay

    End
DECAYCARD

signal_cards["LcToSigmaPPiPiDT"] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-            PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  Sigma+  pi+  pi-                     PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                             PHSP;
    Enddecay

    Decay Sigma+
    1.0000  p+  pi0                              PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                         PHSP;
    Enddecay

    End
DECAYCARD

signal_cards["LcToSigmaPOmegaDT"] = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+  anti-Lambda_c-            PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000  Sigma+  omega                        PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anything                             PHSP;
    Enddecay

    Decay Sigma+
    1.0000  p+  pi0                              PHSP;
    Enddecay

    Decay omega
    1.0000  pi+  pi-  pi0                        OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma                         PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC, one sample per Lambda_c+ decay mode (phase-space Lambda_c+
# decays, matching the paper's phase-space-generated signal samples).
exMC_signal = signal_cards.map do |sample_name, card|
  DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = sample_name
    config.related_dataset = data_4599
    config.events          = 200_000
    config.decay_card      = card
    config.cross_section   = :default
  end
end
exMC_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Common tag-side mode list (all twelve Cabibbo-favored modes) ###
tag_modes = [:LambdacPtoKsP,          # Lambda_c -> p K_S0
             :LambdacPtoKPiP,         # Lambda_c -> p K- pi+
             :LambdacPtoKsPi0P,       # Lambda_c -> p K_S0 pi0
             :LambdacPtoKsPiPiP,      # Lambda_c -> p K_S0 pi+ pi-
             :LambdacPtoKPiPi0P,      # Lambda_c -> p K- pi+ pi0
             :LambdacPtoLambdaPi,     # Lambda_c -> Lambda pi+
             :LambdacPtoLambdaPiPi0,  # Lambda_c -> Lambda pi+ pi0
             :LambdacPtoLambdaPiPiPi, # Lambda_c -> Lambda pi+ pi- pi+
             :LambdacPtoSigma0Pi,     # Lambda_c -> Sigma0 pi+
             :LambdacPtoSigmaPPi0,    # Lambda_c -> Sigma+ pi0
             :LambdacPtoSigmaPPiPi,   # Lambda_c -> Sigma+ pi+ pi-
             :LambdacPtoSigmaPOmega]  # Lambda_c -> Sigma+ omega

### Notes shared by the twelve algorithms ###
common_notes = [
  [:charge_conjugate_modes,
   "Throughout the analysis charge-conjugate modes are implicitly assumed: the " \
   "tagged baryon may be Lambda_c- (charm pinned) or Lambda_c+, and the recoil " \
   "baryon correspondingly Lambda_c+ or Lambda_c-. Only the Lambda_c- tag side " \
   "with the Lambda_c+ recoil final state is declared here; the charge-conjugate " \
   "selection differs only by the sign of every charge."],
  [:tag_side_selection,
   "The tagged Lambda_c candidates are taken from the pre-stored tag collection. " \
   "The underlying charged-track quality requirements are |cos(theta)| < 0.93 and a " \
   "distance of closest approach to the interaction point below 10 cm along the " \
   "beam axis and below 1 cm in the perpendicular plane (tracks used for K_S0 and " \
   "Lambda reconstruction are exempt). Protons are identified when the PID " \
   "likelihoods satisfy L(p) > L(K) and L(p) > L(pi); charged kaons and pions are " \
   "separated by L(K) > L(pi) or L(pi) > L(K)."],
  [:photon_and_pi0_selection,
   "EMC showers not associated with any charged track are photon candidates with " \
   "E > 25 MeV in the barrel (|cos(theta)| < 0.8) and E > 50 MeV in the end-cap " \
   "(0.84 < |cos(theta)| < 0.92), and with an EMC time deviation from the event " \
   "start time within (0, 700) ns. pi0 candidates are built from photon pairs with " \
   "0.115 < M(gamma gamma) < 0.150 GeV/c^2 and are then mass-constrained to the " \
   "nominal pi0 mass, the fitted four-momentum being used in the further analysis."],
  [:ks_and_lambda_selection,
   "K_S0 and Lambda candidates are formed from oppositely charged track pairs " \
   "(pi+ pi- and p pi-). The two daughter tracks must have a distance of closest " \
   "approach to the interaction point within +/- 20 cm along the beam direction (no " \
   "transverse-plane requirement); the charged pion is not subject to PID while " \
   "proton PID is applied. The daughters are constrained to a common decay vertex " \
   "requiring the vertex-fit chi^2 < 100, and the decay vertex must be separated " \
   "from the interaction point by at least twice the fitted vertex resolution. The " \
   "fitted momenta are used. Mass windows: 0.487 < M(pi+ pi-) < 0.511 GeV/c^2 " \
   "(K_S0) and 1.111 < M(p pi-) < 1.121 GeV/c^2 (Lambda), about 3 sigma wide."],
  [:sigma_and_omega_selection,
   "Sigma0, Sigma+ and omega candidates are selected through invariant-mass " \
   "windows: 1.179 < M(Lambda gamma) < 1.203 GeV/c^2 (Sigma0), " \
   "1.176 < M(p pi0) < 1.200 GeV/c^2 (Sigma+) and " \
   "0.760 < M(pi+ pi- pi0) < 0.800 GeV/c^2 (omega)."],
  [:background_veto,
   "For the modes p K_S0 pi0, p K_S0 pi+ pi- and Sigma+ pi+ pi-, background from " \
   "Lambda -> p pi- in the final state is rejected by requiring M(p pi-) outside " \
   "(1.110, 1.120) GeV/c^2. For the mode p K_S0 pi0, candidate events with " \
   "1.170 < M(p pi0) < 1.200 GeV/c^2 are excluded to suppress Sigma+ background. " \
   "For the modes Lambda pi+ pi- pi+, Sigma+ pi0 and Sigma+ pi+ pi-, K_S0 " \
   "candidates are removed by rejecting any pi+ pi- or pi0 pi0 pair whose mass " \
   "falls in (0.480, 0.520) GeV/c^2."],
  [:deltaE_requirements,
   "Mode-dependent Delta E = E - E_beam windows, corresponding to about 3 times " \
   "the resolutions, are applied: p K_S0, p K- pi+, p K_S0 pi+ pi-, Lambda pi+, " \
   "Lambda pi+ pi- pi+ and Sigma0 pi+ use (-0.020, 0.020) GeV; p K_S0 pi0, " \
   "p K- pi+ pi0, Lambda pi+ pi0, Sigma+ pi+ pi- and Sigma+ omega use " \
   "(-0.030, 0.020) GeV; Sigma+ pi0 uses (-0.050, 0.030) GeV. The windows are " \
   "mode-dependent, so they are applied on the stored Delta E in the ROOT analysis " \
   "rather than as a single tag-side window in BOSS."],
  [:mbc_signal_region,
   "The beam-constrained mass is defined as M_BC c^2 = sqrt(E_beam^2 - p^2 c^2) " \
   "with E_beam the beam energy and p the momentum of the Lambda_c candidate " \
   "(nominal Lambda_c mass 2286.46 MeV/c^2). The signal region " \
   "2.276 < M_BC < 2.300 GeV/c^2 is applied on the stored M_BC in the ROOT analysis."],
  [:yield_extraction_and_global_fit,
   "Single-tag and double-tag yields are extracted from unbinned extended maximum " \
   "likelihood fits to the M_BC distributions: the signal shape is taken from " \
   "signal MC convolved with a Gaussian whose parameters are floated, and the " \
   "background is described by an ARGUS function. The 12 x 12 double-tag " \
   "efficiencies eps_ij are evaluated from double-tag signal MC. A global " \
   "least-squares fitter with thirteen free parameters (the twelve branching " \
   "fractions and N(Lambda_c+ Lambda_c-)) accounts for the statistical and " \
   "systematic correlations among the modes; the fit gives " \
   "N(Lambda_c+ Lambda_c-) = (105.9 +/- 4.8 +/- 0.5) x 10^3 and " \
   "chi^2/ndf = 9.9/(24-13) = 0.9. Peaking backgrounds and cross feeds among the " \
   "twelve tag modes are negligible and are not considered in the fit."],
  [:signal_model_reweighting,
   "Phase-space-generated Lambda_c+ decays in the MC are reweighted according to " \
   "the observed behaviour in data; the reweighting factors for the twelve signal " \
   "models are varied within their statistical uncertainties when evaluating the " \
   "systematic uncertainty on the efficiencies."],
  [:mc_statistics_and_intermediate_branching_fractions,
   "Systematic uncertainties include the limited statistics of the MC samples and " \
   "the uncertainties of the intermediate-state branching fractions taken from the " \
   "PDG. The reconstruction uncertainties of the intermediate states are 1.0% per " \
   "pi0, 1.2% for K_S0 and 2.5% for Lambda; tracking and PID uncertainties are " \
   "obtained from control samples of e+e- -> pi+ pi+ pi- pi-, K+ K- pi+ pi- and " \
   "p pbar pi+ pi- above sqrt(s) = 4.0 GeV. Delta E and M_BC requirement " \
   "uncertainties are negligible because the resolutions in MC are corrected to " \
   "accord with those in data."],
  [:photon_arity_limitation,
   "Modes with two pi0 in the final state (Sigma+ pi0 and Sigma+ omega) have four " \
   "photons; the tag-fit signal side consumes at most two photons in this version, " \
   "so for these two modes only the four-momentum constraint and the tag-side " \
   "Lambda_c mass constraint are applied and the photons are declared with an " \
   "event-level shower cap. The remaining photon pairs enter as mass-constrained " \
   "pi0 candidates in the ROOT analysis."]
]

### Event selection (BOSS) — TagAnalysis, one algorithm per signal mode ###

# ---------------------------------------------------------------------------
# Mode 1: Lambda_c+ -> p K_S0
# ---------------------------------------------------------------------------
alg_pks = TagAnalysis.new("LcToPKsDT")
alg_pks.set_header(["LcToPKsDTAlg/LcToPKsDT.h"])
       .set_constant({ "ECMS" => [:double, 4.599] })
       .set_alias({ "std::vector<double>" => "Vdouble" })
       .with_decay_card(signal_cards["LcToPKsDT"])

alg_pks.tag_side(:Lambdac) do |t|   # tagged Lambda_c- (charm pinned)
  t.modes(*tag_modes)
  t.charm -1
end

alg_pks.signal_side do |s|          # recoiling Lambda_c+ -> p K_S0
  s.charged(prp: 1, pip: 1, pim: 1) # p, pi+ and pi- from K_S0
  s.require_charge 1                # +1 (p) + 1 (pi+) - 1 (pi-) = +1
  s.photons 0
end

alg_pks.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.invariant_mass_of(:prp, :pip, :pim).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_pks.note(key, text) }
alg_pks.apply
alg_pks.execute_on([data_4599, incMC_4599] + exMC_signal)

# ---------------------------------------------------------------------------
# Mode 2: Lambda_c+ -> p K- pi+
# ---------------------------------------------------------------------------
alg_pkpi = TagAnalysis.new("LcToPKPiDT")
alg_pkpi.set_header(["LcToPKPiDTAlg/LcToPKPiDT.h"])
        .set_constant({ "ECMS" => [:double, 4.599] })
        .set_alias({ "std::vector<double>" => "Vdouble" })
        .with_decay_card(signal_cards["LcToPKPiDT"])

alg_pkpi.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_pkpi.signal_side do |s|          # recoiling Lambda_c+ -> p K- pi+
  s.charged(prp: 1, km: 1, pip: 1)
  s.require_charge 1                 # +1 (p) - 1 (K-) + 1 (pi+) = +1
  s.photons 0
end

alg_pkpi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:prp, :km, :pip).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_pkpi.note(key, text) }
alg_pkpi.apply
alg_pkpi.execute_on([data_4599, incMC_4599] + exMC_signal)

# ---------------------------------------------------------------------------
# Mode 3: Lambda_c+ -> p K_S0 pi0
# ---------------------------------------------------------------------------
alg_pkspi0 = TagAnalysis.new("LcToPKsPi0DT")
alg_pkspi0.set_header(["LcToPKsPi0DTAlg/LcToPKsPi0DT.h"])
          .set_constant({ "ECMS" => [:double, 4.599] })
          .set_alias({ "std::vector<double>" => "Vdouble" })
          .with_decay_card(signal_cards["LcToPKsPi0DT"])

alg_pkspi0.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_pkspi0.signal_side do |s|        # recoiling Lambda_c+ -> p K_S0 pi0
  s.charged(prp: 1, pip: 1, pim: 1)  # p, pi+ and pi- from K_S0
  s.require_charge 1
  s.photons 2                        # the two photons of the pi0
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_pkspi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.invariant_mass_of(:prp, :pip, :pim, :gamma, :gamma)
   .constrain_to_nominal_mass_of(:"Lambda_c+")
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_pkspi0.note(key, text) }
alg_pkspi0.apply
alg_pkspi0.execute_on([data_4599, incMC_4599] + exMC_signal)

# ---------------------------------------------------------------------------
# Mode 4: Lambda_c+ -> p K_S0 pi+ pi-
# ---------------------------------------------------------------------------
alg_pkspipi = TagAnalysis.new("LcToPKsPiPiDT")
alg_pkspipi.set_header(["LcToPKsPiPiDTAlg/LcToPKsPiPiDT.h"])
           .set_constant({ "ECMS" => [:double, 4.599] })
           .set_alias({ "std::vector<double>" => "Vdouble" })
           .with_decay_card(signal_cards["LcToPKsPiPiDT"])

alg_pkspipi.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_pkspipi.signal_side do |s|       # recoiling Lambda_c+ -> p K_S0 pi+ pi-
  s.charged(prp: 1, pip: 2, pim: 2)  # p, the K_S0 daughters and pi+ pi-
  s.require_charge 1                 # +1 + 2 - 2 = +1
  s.photons 0
end

alg_pkspipi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_pkspipi.note(key, text) }
alg_pkspipi.apply
alg_pkspipi.execute_on([data_4599, incMC_4599] + exMC_signal)

# ---------------------------------------------------------------------------
# Mode 5: Lambda_c+ -> p K- pi+ pi0
# ---------------------------------------------------------------------------
alg_pkpipi0 = TagAnalysis.new("LcToPKPiPi0DT")
alg_pkpipi0.set_header(["LcToPKPiPi0DTAlg/LcToPKPiPi0DT.h"])
           .set_constant({ "ECMS" => [:double, 4.599] })
           .set_alias({ "std::vector<double>" => "Vdouble" })
           .with_decay_card(signal_cards["LcToPKPiPi0DT"])

alg_pkpipi0.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_pkpipi0.signal_side do |s|       # recoiling Lambda_c+ -> p K- pi+ pi0
  s.charged(prp: 1, km: 1, pip: 1)
  s.require_charge 1
  s.photons 2                        # the two photons of the pi0
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_pkpipi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:prp, :km, :pip, :gamma, :gamma)
   .constrain_to_nominal_mass_of(:"Lambda_c+")
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_pkpipi0.note(key, text) }
alg_pkpipi0.apply
alg_pkpipi0.execute_on([data_4599, incMC_4599] + exMC_signal)

# ---------------------------------------------------------------------------
# Mode 6: Lambda_c+ -> Lambda pi+
# ---------------------------------------------------------------------------
alg_lpi = TagAnalysis.new("LcToLambdaPiDT")
alg_lpi.set_header(["LcToLambdaPiDTAlg/LcToLambdaPiDT.h"])
       .set_constant({ "ECMS" => [:double, 4.599] })
       .set_alias({ "std::vector<double>" => "Vdouble" })
       .with_decay_card(signal_cards["LcToLambdaPiDT"])

alg_lpi.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_lpi.signal_side do |s|           # recoiling Lambda_c+ -> Lambda pi+
  s.charged(prp: 1, pim: 1, pip: 1)  # p and pi- from Lambda, plus pi+
  s.require_charge 1
  s.photons 0
end

alg_lpi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
  f.invariant_mass_of(:prp, :pim, :pip).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_lpi.note(key, text) }
alg_lpi.apply
alg_lpi.execute_on([data_4599, incMC_4599] + exMC_signal)

# ---------------------------------------------------------------------------
# Mode 7: Lambda_c+ -> Lambda pi+ pi0
# ---------------------------------------------------------------------------
alg_lpipi0 = TagAnalysis.new("LcToLambdaPiPi0DT")
alg_lpipi0.set_header(["LcToLambdaPiPi0DTAlg/LcToLambdaPiPi0DT.h"])
          .set_constant({ "ECMS" => [:double, 4.599] })
          .set_alias({ "std::vector<double>" => "Vdouble" })
          .with_decay_card(signal_cards["LcToLambdaPiPi0DT"])

alg_lpipi0.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_lpipi0.signal_side do |s|        # recoiling Lambda_c+ -> Lambda pi+ pi0
  s.charged(prp: 1, pim: 1, pip: 1)
  s.require_charge 1
  s.photons 2                        # the two photons of the pi0
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_lpipi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
  f.invariant_mass_of(:prp, :pim, :pip, :gamma, :gamma)
   .constrain_to_nominal_mass_of(:"Lambda_c+")
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_lpipi0.note(key, text) }
alg_lpipi0.apply
alg_lpipi0.execute_on([data_4599, incMC_4599] + exMC_signal)

# ---------------------------------------------------------------------------
# Mode 8: Lambda_c+ -> Lambda pi+ pi- pi+
# ---------------------------------------------------------------------------
alg_lpipipi = TagAnalysis.new("LcToLambdaPiPiPiDT")
alg_lpipipi.set_header(["LcToLambdaPiPiPiDTAlg/LcToLambdaPiPiPiDT.h"])
           .set_constant({ "ECMS" => [:double, 4.599] })
           .set_alias({ "std::vector<double>" => "Vdouble" })
           .with_decay_card(signal_cards["LcToLambdaPiPiPiDT"])

alg_lpipipi.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_lpipipi.signal_side do |s|       # recoiling Lambda_c+ -> Lambda pi+ pi- pi+
  s.charged(prp: 1, pim: 2, pip: 2)  # p and pi- from Lambda, plus pi+ pi-
  s.require_charge 1                 # +1 + 2 - 2 = +1
  s.photons 0
end

alg_lpipipi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_lpipipi.note(key, text) }
alg_lpipipi.apply
alg_lpipipi.execute_on([data_4599, incMC_4599] + exMC_signal)

# ---------------------------------------------------------------------------
# Mode 9: Lambda_c+ -> Sigma0 pi+   (Sigma0 -> Lambda gamma, Lambda -> p pi-)
# ---------------------------------------------------------------------------
alg_s0pi = TagAnalysis.new("LcToSigma0PiDT")
alg_s0pi.set_header(["LcToSigma0PiDTAlg/LcToSigma0PiDT.h"])
        .set_constant({ "ECMS" => [:double, 4.599] })
        .set_alias({ "std::vector<double>" => "Vdouble" })
        .with_decay_card(signal_cards["LcToSigma0PiDT"])

alg_s0pi.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_s0pi.signal_side do |s|          # recoiling Lambda_c+ -> Sigma0 pi+
  s.charged(prp: 1, pim: 1, pip: 1)  # p, pi- (from Lambda) and pi+
  s.require_charge 1
  s.photons 1                        # the photon from Sigma0 -> Lambda gamma
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_s0pi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
  f.invariant_mass_of(:prp, :pim, :gamma).constrain_to_nominal_mass_of(:"Sigma0")
  f.invariant_mass_of(:prp, :pim, :gamma, :pip)
   .constrain_to_nominal_mass_of(:"Lambda_c+")
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_s0pi.note(key, text) }
alg_s0pi.apply
alg_s0pi.execute_on([data_4599, incMC_4599] + exMC_signal)

# ---------------------------------------------------------------------------
# Mode 10: Lambda_c+ -> Sigma+ pi0   (Sigma+ -> p pi0; two pi0 -> four photons)
# ---------------------------------------------------------------------------
alg_sppi0 = TagAnalysis.new("LcToSigmaPPi0DT")
alg_sppi0.set_header(["LcToSigmaPPi0DTAlg/LcToSigmaPPi0DT.h"])
         .set_constant({ "ECMS" => [:double, 4.599] })
         .set_alias({ "std::vector<double>" => "Vdouble" })
         .with_decay_card(signal_cards["LcToSigmaPPi0DT"])

alg_sppi0.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_sppi0.signal_side do |s|         # recoiling Lambda_c+ -> Sigma+ [-> p pi0] pi0
  s.charged(prp: 1)
  s.require_charge 1                 # proton from Sigma+ -> p pi0
  s.photons 2..4                     # two pi0 -> up to four showers
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_sppi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_sppi0.note(key, text) }
alg_sppi0.apply
alg_sppi0.execute_on([data_4599, incMC_4599] + exMC_signal)

# ---------------------------------------------------------------------------
# Mode 11: Lambda_c+ -> Sigma+ pi+ pi-   (Sigma+ -> p pi0)
# ---------------------------------------------------------------------------
alg_sppipi = TagAnalysis.new("LcToSigmaPPiPiDT")
alg_sppipi.set_header(["LcToSigmaPPiPiDTAlg/LcToSigmaPPiPiDT.h"])
          .set_constant({ "ECMS" => [:double, 4.599] })
          .set_alias({ "std::vector<double>" => "Vdouble" })
          .with_decay_card(signal_cards["LcToSigmaPPiPiDT"])

alg_sppipi.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_sppipi.signal_side do |s|        # recoiling Lambda_c+ -> Sigma+ pi+ pi-
  s.charged(prp: 1, pip: 1, pim: 1)  # the proton of Sigma+ plus pi+ pi-
  s.require_charge 1
  s.photons 2                        # the two photons of the pi0 from Sigma+
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_sppipi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:prp, :gamma, :gamma).constrain_to_nominal_mass_of(:"Sigma+")
  f.invariant_mass_of(:prp, :gamma, :gamma, :pip, :pim)
   .constrain_to_nominal_mass_of(:"Lambda_c+")
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_sppipi.note(key, text) }
alg_sppipi.apply
alg_sppipi.execute_on([data_4599, incMC_4599] + exMC_signal)

# ---------------------------------------------------------------------------
# Mode 12: Lambda_c+ -> Sigma+ omega   (Sigma+ -> p pi0, omega -> pi+ pi- pi0;
# two pi0 -> four photons)
# ---------------------------------------------------------------------------
alg_spomega = TagAnalysis.new("LcToSigmaPOmegaDT")
alg_spomega.set_header(["LcToSigmaPOmegaDTAlg/LcToSigmaPOmegaDT.h"])
           .set_constant({ "ECMS" => [:double, 4.599] })
           .set_alias({ "std::vector<double>" => "Vdouble" })
           .with_decay_card(signal_cards["LcToSigmaPOmegaDT"])

alg_spomega.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_spomega.signal_side do |s|       # recoiling Lambda_c+ -> Sigma+ omega
  s.charged(prp: 1, pip: 1, pim: 1)  # p from Sigma+ and pi+ pi- from omega
  s.require_charge 1
  s.photons 2..4                     # two pi0 -> up to four showers
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_spomega.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

common_notes.each { |key, text| alg_spomega.note(key, text) }
alg_spomega.apply
alg_spomega.execute_on([data_4599, incMC_4599] + exMC_signal)
