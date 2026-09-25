# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# √s = 4.600 GeV : real data and the matching inclusive MC (BOSS 7.0.3, 4599.53 MeV)
data_4600  = DatasetManager.real_data.find("703_4600")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

# Decay card for the signal: e+e- -> Lambda_c+ Lambda_c-,
# with the single-tag side Lambda_c- -> pbar K+ pi- and the signal side Lambda_c+ -> K_S0 X
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c-          PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K+ pi-                    PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 K_S0 X0                           PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                           PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive signal MC
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "lambdac_ks_inclusive_4600"
  config.related_dataset = data_4600
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) — tag-based (single-tag Lambda_c- + signal side) ###
alg = TagAnalysis.new("LambdacIncKs")
alg.set_header(["LambdacIncKsAlg/LambdacIncKs.h"])
   .set_constant({ "ECMS" => [:double, 4.600] })
   .with_decay_card(decay_card_signal)
   # Signal-side K_S0 is built from the two remaining pions with a secondary vertex fit;
   # this reconstruction step has no dedicated tag-DSL block.
   .note(:ks0_reconstruction,
         "signal-side K_S0 formed from the remaining pi+ pi- pair through a secondary
          vertex fit requiring chi2 < 100, decay vertex displaced from the interaction
          point by more than 2 sigma, |cos(theta)| < 0.93, |Vz| < 20 cm and no Vr
          constraint; when several combinations pass, the smallest-chi2 candidate is kept")
   # Mode-dependent tag-side background vetoes.
   .note(:background_veto,
         "reject Lambda via M(pbar pi+) in [1.110, 1.120] GeV for the pbar K_S0 pi0 and
          pbar K_S0 pi+ pi- modes; reject K_S0 via M(pi+ pi-) or M(pi0 pi0) in
          [0.480, 0.520] GeV for the Lambda pi- pi+ pi- mode; reject Sigma- via
          M(pbar pi0) in [1.170, 1.200] GeV for the pbar K_S0 pi0 mode")

# Single-tag side: reconstruct the Lambda_c- (charm -1) through the eight available
# hadronic tag modes (the pbar Sigma0 pi-, Sigma- pi0 and Sigma- pi+ pi- modes are
# unavailable and therefore omitted).
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP,        # pbar K_S0
          :LambdacPtoKPiP,       # pbar K+ pi-
          :LambdacPtoKsPPi0,     # pbar K_S0 pi0
          :LambdacPtoKsPPiPi,    # pbar K_S0 pi+ pi-
          :LambdacPtoKPiPPi0,    # pbar K+ pi- pi0
          :LambdacPtoLambdaPi,   # Lambda pi-
          :LambdacPtoLambdaPiPi0,# Lambda pi- pi0
          :LambdacPtoLambdaPiPiPi # Lambda pi- pi+ pi-
  t.charm -1
end

# Signal side: exactly one pi+ and one pi- with net charge 0, plus a massless
# missing X0 accounting for the undetected remainder of the Lambda_c+ decay.
alg.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.require_charge(0)
  s.missing :X0, mass: nil   # massless missing particle
end

# 4C kinematic fit. The K_S0 mass is not constrained here — it is selected later by a
# two-dimensional unbinned maximum-likelihood fit to beam-constrained mass vs M(pi+ pi-),
# so the chi2 cut is kept effectively unbounded.
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 999999
end

alg.apply

# Execute on real data, inclusive MC and the exclusive signal MC
root_files = alg.execute_on([data_4600, incMC_4600, exMC_signal])