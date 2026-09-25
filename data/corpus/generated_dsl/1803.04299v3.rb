# =====================================================================
# BOSS DSL spec: absolute branching fractions of
#   Lambda_c+ -> Xi0 K+        and
#   Lambda_c+ -> Xi(1530)0 K+
# via the double-tag technique at sqrt(s) = 4.5995 GeV (567 pb^-1)
#
# Tag side  : ST anti-Lambda_c- in 10 hadronic modes
# Signal side: exactly one K+ + one missing Xi0 / Xi(1530)0
# =====================================================================

### ------------------------------------------------------------------
### Datasets  (BOSS 7.0.3, CMS energy 4599.53 MeV -> sample 703_4600)
### ------------------------------------------------------------------
data_4600  = DatasetManager.real_data.find("703_4600")        # 567 pb^-1 of real data
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")     # inclusive MC at 4.6 GeV

### ------------------------------------------------------------------
### Decay cards (EvtGen syntax)
### ------------------------------------------------------------------
# Signal mode A: Lambda_c+ -> Xi0 K+
decay_card_xi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 Xi0 K+ PHSP;
    Enddecay

    Decay Xi0
    1.0000 Lambda0 pi0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Signal mode B: Lambda_c+ -> Xi(1530)0 K+
decay_card_xi1530 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 Xi(1530)0 K+ PHSP;
    Enddecay

    Decay Xi(1530)0
    1.0000 Xi0 pi0 PHSP;
    Enddecay

    Decay Xi0
    1.0000 Lambda0 pi0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

### ------------------------------------------------------------------
### Exclusive MC (200k events per signal mode)
### ------------------------------------------------------------------
exMC_xi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_4600_LcToXi0K"
  config.related_dataset = data_4600
  config.events         = 200_000
  config.decay_card     = decay_card_xi0
  config.cross_section  = :default
end

exMC_xi1530 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_4600_LcToXi1530K"
  config.related_dataset = data_4600
  config.events         = 200_000
  config.decay_card     = decay_card_xi1530
  config.cross_section  = :default
end

### ------------------------------------------------------------------
### Shared ST tag-side mode list (anti-Lambda_c- tag)
### (Lambda_c+ channel names; charm -1 selects the anti-Lambda_c- tag)
### The Sigma_bar- pi0 and Sigma_bar- pi- pi+ modes have no DTagAlg channel.
### ------------------------------------------------------------------
tag_modes = [
  :LambdacPtoPKs,          # anti-Lc- -> anti-p K_S0
  :LambdacPtoPKPi,         # anti-Lc- -> anti-p K+ pi-
  :LambdacPtoPKsPi0,       # anti-Lc- -> anti-p K_S0 pi0
  :LambdacPtoPKsPiPi,      # anti-Lc- -> anti-p K_S0 pi- pi+
  :LambdacPtoPKPiPi0,      # anti-Lc- -> anti-p K+ pi- pi0
  :LambdacPtoPPiPi,        # anti-Lc- -> anti-p pi- pi+
  :LambdacPtoLambdaPi,     # anti-Lc- -> anti-Lambda0 pi-
  :LambdacPtoLambdaPiPi0,  # anti-Lc- -> anti-Lambda0 pi- pi0
  :LambdacPtoLambdaPiPiPi, # anti-Lc- -> anti-Lambda0 pi- pi+ pi-
  :LambdacPtoSigma0Pi      # anti-Lc- -> anti-Sigma0 pi-
]

### ==================================================================
### Signal mode A: Lambda_c+ -> Xi0 K+
### ==================================================================
alg_xi0 = TagAnalysis.new("LcToXi0K")
alg_xi0.set_header(["LcToXi0KAlg/LcToXi0K.h"])
       .set_constant({"ECMS" => [:double, 4.5995]})   # E_CMS = 4.5995 GeV
       .with_decay_card(decay_card_xi0)

# ---- ST tag side: anti-Lambda_c- in the 10 hadronic modes ----
alg_xi0.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm(-1)                                     # tagged side = anti-Lambda_c-
  t.window :deltaE, min: -0.050, max: 0.030       # ~3 sigma DeltaE window
  t.window :mBC,   min:  2.282, max: 2.291        # ST yield region in mBC
end

# ---- Signal side: exactly one K+ opposite the tag + missing Xi0 ----
alg_xi0.signal_side do |s|
  s.charged(kp: 1)              # exactly one K+ on the signal side
  s.require_charge(1)           # total signal-side charge = +1
  s.missing :Xi0, mass: 1.31486 # undetected Xi0 inferred from the missing mass
end

# ---- Kinematic fit: 4-momentum conservation ----
alg_xi0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# ---- BOSS-side procedures that cannot be expressed in the tag DSL ----
alg_xi0
  .note(:tag_reconstruction_cuts,
        "Tag-side selections handled inside DTagAlg and not expressible through the tag DSL: " \
        "charged tracks |cos(theta)|<0.93, |Vz|<10 cm, Vr<1 cm; likelihood-ratio PID " \
        "(protons L(p)>L(K) and L(p)>L(pi); kaons L(K)>L(pi); pions L(pi)>L(K)); " \
        "photons E>25 MeV (barrel) / 50 MeV (endcap) with EMC time 0-700 ns; " \
        "pi0 -> gamma gamma within 115-150 MeV/c^2 using a 1C mass-constraint fit; " \
        "K_S0 / Lambda secondary vertices with fit chi^2<100 and decay vertex >2 sigma from the IP, " \
        "K_S0 mass 487-511 MeV/c^2, anti-Lambda0 mass 1111-1121 MeV/c^2.")
  .note(:mode_dependent_deltaE,
        "The DeltaE requirement is ~3 sigma and mode dependent (-50 to +30 MeV); a single " \
        "tag-side window :deltaE (-0.050, 0.030) GeV is declared here, the per-mode windows are " \
        "applied in the ROOT analysis.")
  .note(:multiple_dt_candidates,
        "Events containing more than one double-tag candidate are rejected; candidate ranking and " \
        "selection are delegated to DTagTool.")

alg_xi0.apply                       # no Selection argument for a tag analysis
alg_xi0.execute_on([data_4600, incMC_4600, exMC_xi0])

### ==================================================================
### Signal mode B: Lambda_c+ -> Xi(1530)0 K+   (identical selection chain)
### ==================================================================
alg_xi1530 = TagAnalysis.new("LcToXi1530K")
alg_xi1530.set_header(["LcToXi1530KAlg/LcToXi1530K.h"])
          .set_constant({"ECMS" => [:double, 4.5995]})
          .with_decay_card(decay_card_xi1530)

# ---- ST tag side: same 10 hadronic anti-Lambda_c- modes ----
alg_xi1530.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm(-1)
  t.window :deltaE, min: -0.050, max: 0.030
  t.window :mBC,   min:  2.282, max: 2.291
end

# ---- Signal side: same chain, Xi(1530)0 missing-mass hypothesis ----
alg_xi1530.signal_side do |s|
  s.charged(kp: 1)
  s.require_charge(1)
  s.missing :Xi1530, mass: 1.5318   # undetected Xi(1530)0 inferred from the missing mass
end

# ---- Kinematic fit: identical 4-momentum constraint ----
alg_xi1530.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_xi1530
  .note(:tag_reconstruction_cuts,
        "Tag-side selections handled inside DTagAlg and not expressible through the tag DSL: " \
        "charged tracks |cos(theta)|<0.93, |Vz|<10 cm, Vr<1 cm; likelihood-ratio PID " \
        "(protons L(p)>L(K) and L(p)>L(pi); kaons L(K)>L(pi); pions L(pi)>L(K)); " \
        "photons E>25 MeV (barrel) / 50 MeV (endcap) with EMC time 0-700 ns; " \
        "pi0 -> gamma gamma within 115-150 MeV/c^2 using a 1C mass-constraint fit; " \
        "K_S0 / Lambda secondary vertices with fit chi^2<100 and decay vertex >2 sigma from the IP, " \
        "K_S0 mass 487-511 MeV/c^2, anti-Lambda0 mass 1111-1121 MeV/c^2.")
  .note(:mode_dependent_deltaE,
        "The DeltaE requirement is ~3 sigma and mode dependent (-50 to +30 MeV); a single " \
        "tag-side window :deltaE (-0.050, 0.030) GeV is declared here, the per-mode windows are " \
        "applied in the ROOT analysis.")
  .note(:multiple_dt_candidates,
        "Events containing more than one double-tag candidate are rejected; candidate ranking and " \
        "selection are delegated to DTagTool.")

alg_xi1530.apply
alg_xi1530.execute_on([data_4600, incMC_4600, exMC_xi1530])