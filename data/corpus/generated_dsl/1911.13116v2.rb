# ============================================================
# BOSS part — psi(3770) data
# Single-tag analysis: D- -> K+ pi- pi-  (tag side)
# BNV signal on the other D+: D+ -> Lambda e+, Sigma0 e+,
#                             Lambdabar e+, Sigmabar0 e+
# with Lambda -> p pi-, Lambdabar -> pbar pi+
# The tag is taken from the pre-stored EvtRecDTag collection (DTagTool);
# signal side = tracks/showers not used by the tag.
# ============================================================

### Dataset preparation ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # 2.93 fb^-1 e+e- data at sqrt(s) = 3.773 GeV
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC sample

# ------------------- decay cards (EvtGen) -------------------
# The tagged D- always decays to K+ pi- pi-

decay_card_Lambda_e = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay D+
    1.0000 Lambda0 e+ PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    End
DECAYCARD

decay_card_Sigma0_e = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay D+
    1.0000 Sigma0 e+ PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    End
DECAYCARD

decay_card_Lambdabar_e = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay D+
    1.0000 anti-Lambda0 e+ PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

decay_card_Sigmabar0_e = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay D+
    1.0000 anti-Sigma0 e+ PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 gamma anti-Lambda0 PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# ---------- exclusive MC: 100k events for each of the four modes ----------
exMC_Lambda_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_DpToLambdaE"
  config.related_dataset = psi3770_data
  config.events          = 100_000
  config.decay_card      = decay_card_Lambda_e
  config.cross_section   = :default
end

exMC_Sigma0_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_DpToSigma0E"
  config.related_dataset = psi3770_data
  config.events          = 100_000
  config.decay_card      = decay_card_Sigma0_e
  config.cross_section   = :default
end

exMC_Lambdabar_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_DpToLambdabarE"
  config.related_dataset = psi3770_data
  config.events          = 100_000
  config.decay_card      = decay_card_Lambdabar_e
  config.cross_section   = :default
end

exMC_Sigmabar0_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_DpToSigmabar0E"
  config.related_dataset = psi3770_data
  config.events          = 100_000
  config.decay_card      = decay_card_Sigmabar0_e
  config.cross_section   = :default
end

### Event selection (BOSS) — tag-based ###
# Common tag side for all four modes: D- -> K+ pi- pi-
# (charm = -1 pins the tagged side to the D-).

# ---- Mode 1: D+ -> Lambda e+, Lambda -> p pi- ----
alg_Lambda_e = TagAnalysis.new("DpToLambdaE")
alg_Lambda_e.set_header(["DpToLambdaEAlg/DpToLambdaE.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            .note(:track_quality_cuts, "signal-side charged tracks required |cos(theta)|<0.93, "
                                       "|Vz|<20 cm (loose, for the Lambda daughters), Vr<1 cm")
            .note(:signal_pid, "probability-method PID with CL>0.001: proton identified against "
                               "K and pi (exactly one p); e+ identified against pi, K, p "
                               "(at least one e+)")
            .note(:secondary_vertex_fit, "Lambda built from p pi- by a secondary-vertex fit with "
                                         "mass-difference minimisation; the tag layer consumes the "
                                         "remaining tracks directly, so this step is not emitted")
alg_Lambda_e.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi      # DTagAlg channel name for D+ -> K- pi+ pi+ (charge conjugate taken via charm)
  t.charm -1              # tag the D-
end
alg_Lambda_e.signal_side do |s|
  s.charged(ep: 1, prp: 1, pim: 1, at_least: true)   # e+, p, pi- (>= 2 positive and >= 1 negative track)
  s.require_charge 1
end
alg_Lambda_e.fit do |f|
  f.constrain_four_momentum     # 4C: tag + signal = measured CMS four-momentum
  f.chi2_cut 200                # loose BOSS-level cut; tight cut applied in ROOT
end
alg_Lambda_e.apply
alg_Lambda_e.execute_on([psi3770_data, psi3770_incMC, exMC_Lambda_e])

# ---- Mode 2: D+ -> Sigma0 e+, Sigma0 -> gamma Lambda, Lambda -> p pi- ----
alg_Sigma0_e = TagAnalysis.new("DpToSigma0E")
alg_Sigma0_e.set_header(["DpToSigma0EAlg/DpToSigma0E.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            .note(:track_quality_cuts, "signal-side charged tracks required |cos(theta)|<0.93, "
                                       "|Vz|<20 cm (loose, for the Lambda daughters), Vr<1 cm")
            .note(:signal_pid, "probability-method PID with CL>0.001: proton identified against "
                               "K and pi (exactly one p); e+ identified against pi, K, p "
                               "(at least one e+)")
            .note(:secondary_vertex_fit, "Lambda built from p pi- by a secondary-vertex fit with "
                                         "mass-difference minimisation; the tag layer consumes the "
                                         "remaining tracks directly, so this step is not emitted")
alg_Sigma0_e.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi
  t.charm -1
end
alg_Sigma0_e.signal_side do |s|
  s.photons 1                                        # at least one photon (Sigma0 -> gamma Lambda)
  s.min_photon_angle 10.0                            # photon at least 10 deg from any track
  s.charged(ep: 1, prp: 1, pim: 1, at_least: true)
  s.require_charge 1
end
alg_Sigma0_e.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :prp, :pim).constrain_to_nominal_mass_of(:Sigma0)  # 1C: M(Lambda gamma) = m_Sigma0
  f.chi2_cut 200
end
alg_Sigma0_e.apply
alg_Sigma0_e.execute_on([psi3770_data, psi3770_incMC, exMC_Sigma0_e])

# ---- Mode 3: D+ -> Lambdabar e+, Lambdabar -> pbar pi+ ----
alg_Lambdabar_e = TagAnalysis.new("DpToLambdabarE")
alg_Lambdabar_e.set_header(["DpToLambdabarEAlg/DpToLambdabarE.h"])
               .set_constant({"ECMS" => [:double, 3.773]})
               .note(:track_quality_cuts, "signal-side charged tracks required |cos(theta)|<0.93, "
                                          "|Vz|<20 cm (loose, for the Lambdabar daughters), Vr<1 cm")
               .note(:signal_pid, "probability-method PID with CL>0.001: antiproton identified "
                                  "against K and pi (exactly one pbar); e+ identified against pi, K, p "
                                  "(at least one e+)")
               .note(:secondary_vertex_fit, "Lambdabar built from pbar pi+ by a secondary-vertex fit "
                                            "with mass-difference minimisation; the tag layer consumes "
                                            "the remaining tracks directly, so this step is not emitted")
alg_Lambdabar_e.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi
  t.charm -1
end
alg_Lambdabar_e.signal_side do |s|
  s.charged(ep: 1, prm: 1, pip: 1, at_least: true)   # e+, pbar, pi+
  s.require_charge 1
end
alg_Lambdabar_e.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end
alg_Lambdabar_e.apply
alg_Lambdabar_e.execute_on([psi3770_data, psi3770_incMC, exMC_Lambdabar_e])

# ---- Mode 4: D+ -> Sigmabar0 e+, Sigmabar0 -> gamma Lambdabar, Lambdabar -> pbar pi+ ----
alg_Sigmabar0_e = TagAnalysis.new("DpToSigmabar0E")
alg_Sigmabar0_e.set_header(["DpToSigmabar0EAlg/DpToSigmabar0E.h"])
               .set_constant({"ECMS" => [:double, 3.773]})
               .note(:track_quality_cuts, "signal-side charged tracks required |cos(theta)|<0.93, "
                                          "|Vz|<20 cm (loose, for the Lambdabar daughters), Vr<1 cm")
               .note(:signal_pid, "probability-method PID with CL>0.001: antiproton identified "
                                  "against K and pi (exactly one pbar); e+ identified against pi, K, p "
                                  "(at least one e+)")
               .note(:secondary_vertex_fit, "Lambdabar built from pbar pi+ by a secondary-vertex fit "
                                            "with mass-difference minimisation; the tag layer consumes "
                                            "the remaining tracks directly, so this step is not emitted")
alg_Sigmabar0_e.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi
  t.charm -1
end
alg_Sigmabar0_e.signal_side do |s|
  s.photons 1                                         # at least one photon (Sigmabar0 -> gamma Lambdabar)
  s.min_photon_angle 10.0
  s.charged(ep: 1, prm: 1, pip: 1, at_least: true)
  s.require_charge 1
end
alg_Sigmabar0_e.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :prm, :pip).constrain_to_nominal_mass_of(:Sigmabar0)  # 1C: M(Lambdabar gamma) = m_Sigmabar0
  f.chi2_cut 200
end
alg_Sigmabar0_e.apply
alg_Sigmabar0_e.execute_on([psi3770_data, psi3770_incMC, exMC_Sigmabar0_e])