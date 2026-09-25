# =============================================================================
# e+e- -> Lambda_c+ anti-Lambda_c-   (psi(4260) energy region)
# Single-tag (ST) analysis of the anti-Lambda_c- reconstructed in twelve
# hadronic decay modes; the recoiling Lambda_c+ signal side is left to the
# downstream ROOT analysis.  Tag-based surface -> TagAnalysis (no Selection).
# =============================================================================

### Dataset preparation ###

# Real data: 4.5 fb^-1 collected at seven c.m. energy points (4599.53-4698.82 MeV)
data_points = [
  DatasetManager.real_data.find("703_4600"),   # 4599.53 MeV
  DatasetManager.real_data.find("706_4610"),   # 4611.86 MeV
  DatasetManager.real_data.find("706_4620"),   # 4628.00 MeV
  DatasetManager.real_data.find("706_4640"),   # 4640.91 MeV
  DatasetManager.real_data.find("706_4660"),   # 4661.24 MeV
  DatasetManager.real_data.find("706_4680"),   # 4681.92 MeV
  DatasetManager.real_data.find("706_4700"),   # 4698.82 MeV
]

# Matching inclusive MC, one sample per energy point
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
]

# Decay card for the signal process (EvtGen format, KKMC top mother psi(4260))
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000  Lambda_c+  anti-Lambda_c-           PHSP;
  Enddecay

  Decay Lambda_c+
  0.08333 p+  K_S0                     PHSP;
  0.08333 p+  K-  pi+                  PHSP;
  0.08333 p+  K_S0  pi0                PHSP;
  0.08333 p+  K_S0  pi+  pi-           PHSP;
  0.08333 p+  K-  pi+  pi0             PHSP;
  0.08333 Lambda0  pi+                 PHSP;
  0.08333 Lambda0  pi+  pi0            PHSP;
  0.08333 Lambda0  pi+  pi-  pi+       PHSP;
  0.08333 Sigma0  pi+                  PHSP;
  0.08333 Sigma+  pi0                  PHSP;
  0.08333 Sigma+  pi+  pi-             PHSP;
  0.08333 p+  pi+  pi-                 PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.000  anti-p-  K+  pi-               PHSP;
  Enddecay

  Decay Sigma0
  1.000  Lambda0  gamma                 PHSP;
  Enddecay

  Decay Sigma+
  0.516  p+  pi0                        PHSP;
  0.484  n0  pi+                        PHSP;
  Enddecay

  Decay Lambda0
  1.000  p+  pi-                        HypWK;
  Enddecay

  Decay K_S0
  1.000  pi+  pi-                       VSS;
  Enddecay

  Decay pi0
  1.000  gamma  gamma                   PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive signal MC: 500k events per energy point, shared decay card, one MC per point
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lambdac_st"   # auto-suffixed per energy point
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Tag-based event selection (TagAnalysis replaces Algorithm + Selection) ###
alg = TagAnalysis.new("LambdacST")
alg.set_header(["LambdacSTAlg/LambdacST.h"])
   .set_constant({"ECMS" => [:double, 4.640]})  # nominal energy (per-run beam E from DB)
   .with_decay_card(decay_card_signal)

# Single tag: the anti-Lambda_c- reconstructed in the twelve hadronic channels
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP, :LambdacPtoKPiP, :LambdacPtoKsPi0P, :LambdacPtoKsPiPiP,
          :LambdacPtoKPiPi0P, :LambdacPtoLambdaPi, :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi, :LambdacPtoSigma0Pi, :LambdacPtoSigmaPi0,
          :LambdacPtoSigmaPiPi, :LambdacPtoPPiPi
  t.charm -1   # pin the tagged side to the anti-Lambda_c-
end

# Signal side: recoiling Lambda_c+ (not reconstructed -> ROOT) plus 0-48 photons
alg.signal_side do |s|
  s.photons 0..48
  s.missing :"Lambda_c+"    # untagged recoil, stored for the downstream ROOT analysis
end

# 4C kinematic fit: four-momentum conservation + tag mass constrained to m(Lambda_c)
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

# Per-mode tag-side intermediate-resonance vetoes (not expressible in the tag DSL)
alg.note(:background_veto,
         "Tag-side intermediate-resonance vetoes applied per channel: " \
         "reject M(anti-p pi+) in [1.17, 1.20] GeV/c^2 for the anti-p K_S0 pi0, " \
         "anti-p K_S0 pi- pi+, anti-Sigma- pi- pi+ and anti-p pi- pi+ modes; " \
         "reject M(pi+pi-) or M(pi0pi0) in [0.48, 0.52] GeV/c^2 for the " \
         "anti-Lambda0 pi- pi+ pi-, anti-Sigma- pi0, anti-Sigma- pi- pi+ and " \
         "anti-p pi- pi+ modes; require M(anti-p pi0) outside [1.17, 1.20] GeV/c^2 " \
         "for the anti-p K_S0 pi0 mode.")
    .note(:tag_candidate_selection,
          "For events with more than one tag candidate the one with the smallest " \
          "|DeltaE| is kept; single-tag candidate ranking is not expressible in the " \
          "tag DSL. mBC, DeltaE and the missing mass are stored and windowed in ROOT.")

alg.apply                                    # no Selection argument for TagAnalysis
root_files = alg.execute_on(data_points + incMC_points + exMCs)