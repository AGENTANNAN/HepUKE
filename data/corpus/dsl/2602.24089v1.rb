# Analysis: Precision studies and searches for CP asymmetries in the
# inclusive decay Lambda_c+ -> Lambda X (BESIII).
# Method: double-tag technique with a hadronic tag on the Lambda_c-bar
# ("tag side") reconstructed in one of eleven single-tag modes.  The signal
# side is inclusive Lambda_c+ -> Lambda X, with Lambda -> p pi- selected from
# the remaining tracks.  Data: 4.5 fb-1 of e+e- annihilation at
# sqrt(s) = 4.600 - 4.699 GeV.

### Dataset description ###
# Seven c.m. energy points where Lambda_c+ Lambda_c-bar pairs are produced
xyz_datasets = [
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("706_4610"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
]

xyz_incMCs = [
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
]

# ---------------------------------------------------------------------------
# Signal decay card: Lambda_c+ -> Lambda X (inclusive), Lambda_c-bar tagged
# via one of eleven hadronic modes.  Below only one representative tag mode
# is written into the card; the exclusive-MC decay card is provided to
# generate the signal-side inclusive Lambda_c+ -> Lambda X and any tag mode.
# ---------------------------------------------------------------------------
decay_card_signal = <<~DECAYCARD
    Decay Lambda_c+
    1.0000 Lambda0 anti-Lambda_c-              PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K+ pi-                      PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                              PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: signal MC with tag-side Lambda_c-bar decaying into any of the
# eleven ST modes and signal-side Lambda_c+ decaying inclusively into modes
# containing a Lambda in the final state.
exMC_signal = DatasetManager.create_exclusive_mc_for(xyz_datasets) do |config|
  config.sample_name   = "LcpToLambdaX_DT_signalMC"
  config.events        = 500000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMC_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) — Tag-based analysis ###
alg_name = "LcpLambdaInclusiveTag"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.682] })
   .set_alias({ "std::vector<double>" => "Vdouble" })

# ---------------------------------------------------------------------------
# Tag side: Lambda_c-bar reconstructed in one of 11 hadronic modes.
#   pbar K_S0                pbar K+ pi-           pbar K_S0 pi0
#   pbar K_S0 pi- pi+        pbar K+ pi- pi0       Lambdabar pi-
#   Lambdabar pi- pi0        Lambdabar pi- pi+ pi- Sigmabar0 pi-
#   Sigmabar- pi0            Sigmabar- pi- pi+
# ---------------------------------------------------------------------------
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP,
          :LambdacPtoKPiP,
          :LambdacPtoKsPi0P,
          :LambdacPtoKsPiPiP,
          :LambdacPtoKPiPi0P,
          :LambdacPtoLambdaPi,
          :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi,
          :LambdacPtoSigma0Pi,
          :LambdacPtoSigmaMPi0,
          :LambdacPtoSigmaMPiPi
  t.charm(-1)     # tag the Lambda_c-bar
end

# Signal side: inclusive Lambda_c+ -> Lambda X, with Lambda -> p pi-.  Only the
# two charged tracks from the Lambda decay (proton + pi-) are declared here;
# the remaining hadronic content of X is captured as an inclusive massless
# "X0" missing placeholder so that a 4C constraint remains meaningful in the
# fit derivation.
alg.signal_side do |s|
  s.charged(prp: 1, pim: 1, at_least: true)
  s.missing :X0, mass: nil        # inclusive "X" (any allowed final state)
end

# Tag-side kinematic fit: constrain the total 4-momentum to ECMS.
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg
  .note(:tag_modes_full_list,
        "Lambda_c-bar single-tag modes used (charge-conjugation implicit): " \
        "pbar K_S0, pbar K+ pi-, pbar K_S0 pi0, pbar K_S0 pi- pi+, " \
        "pbar K+ pi- pi0, Lambdabar pi-, Lambdabar pi- pi0, " \
        "Lambdabar pi- pi+ pi-, Sigmabar0 pi-, Sigmabar- pi0, " \
        "Sigmabar- pi- pi+.  ST selection follows the procedure of BESIII " \
        "Ref. [37] (Phys. Rev. D 106, 072002).")
  .note(:signal_side_lambda_reconstruction,
        "Signal-side Lambda candidate reconstructed from the remaining " \
        "tracks via Lambda -> p pi-: proton and pi- from the tag's " \
        "otherTracks() combined by a secondary-vertex fit; same selection " \
        "criteria as for the tag-side Lambda candidates.")
  .note(:truth_matching,
        "MC truth-matching angles: theta_match(p) < 15 deg and " \
        "theta_match(pi-) < 25 deg on both ST and DT sides; " \
        "theta_match < 10 deg for other particles on the ST side.  Failing " \
        "the requirement classifies the event as 'unmatched'.")
  .note(:bdt_reweighting,
        "Signal MC samples reweighted with a BDT to correct the recoil mass " \
        "against Lambda and the Lambda momentum distributions.")
  .note(:lambda_reconstruction_correction,
        "Lambda reconstruction efficiency corrected using J/psi -> pbar K+ " \
        "Lambda and J/psi -> Lambda Lambdabar control samples, with per-bin " \
        "corrections in Lambda momentum and cos(theta).")

alg.with_decay_card(decay_card_signal).apply
alg.execute_on(xyz_datasets + xyz_incMCs + exMC_signal)
