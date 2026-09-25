# =====================================================================
# Λc+ single-tag analysis: B(Λc+ → Λ K+) relative to the reference mode
# Λc+ → Λ π+ (Λ → p π- in both), produced via e+e- → Λc+ anti-Λc-
# with the anti-Λc- side treated inclusively.
# Tag-based surface (TagAnalysis): the tag Λc+ comes from the pre-stored
# DTag candidates; the tag carries its own selected, PID'd tracks and
# showers, so no select_track / select_photon / pid / Selection is used.
# =====================================================================

### Dataset description ###
# 13 c.m. energy points from 4.599 to 4.950 GeV (6.44 fb^-1)
data_points = [
  DatasetManager.real_data.find("703_4600"),   # 4599.53 MeV
  DatasetManager.real_data.find("706_4610"),   # 4611.86 MeV
  DatasetManager.real_data.find("706_4620"),   # 4628.00 MeV
  DatasetManager.real_data.find("706_4640"),   # 4640.91 MeV
  DatasetManager.real_data.find("706_4660"),   # 4661.24 MeV
  DatasetManager.real_data.find("706_4680"),   # 4681.92 MeV
  DatasetManager.real_data.find("706_4700"),   # 4698.82 MeV
  DatasetManager.real_data.find("707_4740"),   # 4739.70 MeV
  DatasetManager.real_data.find("707_4750"),   # 4750.05 MeV
  DatasetManager.real_data.find("707_4780"),   # 4780.54 MeV
  DatasetManager.real_data.find("707_4840"),   # 4843.07 MeV
  DatasetManager.real_data.find("707_4914"),   # 4918.02 MeV
  DatasetManager.real_data.find("707_4946")    # 4950.93 MeV
]

# Corresponding inclusive MC samples (one per energy point)
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
  DatasetManager.inclusive_mc.find("707_4740"),
  DatasetManager.inclusive_mc.find("707_4750"),
  DatasetManager.inclusive_mc.find("707_4780"),
  DatasetManager.inclusive_mc.find("707_4840"),
  DatasetManager.inclusive_mc.find("707_4914"),
  DatasetManager.inclusive_mc.find("707_4946")
]

# Decay card for the signal mode Λc+ → Λ K+ (Λ → p π-); the anti-Λc- side is inclusive
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay
  Decay Lambda_c+
  1.0000 Lambda0 K+ PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay
  End
DECAYCARD

# Decay card for the reference mode Λc+ → Λ π+ (Λ → p π-)
decay_card_reference = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay
  Decay Lambda_c+
  1.0000 Lambda0 pi+ PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay
  End
DECAYCARD

# 500k-event exclusive MC for each of the two modes, one sample per energy point
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "LambdacToLambdaK"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

exMC_reference = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "LambdacToLambdaPi"
  config.events        = 500_000
  config.decay_card    = decay_card_reference
  config.cross_section = :default
end

### Tag-based event selection (BOSS) ###
alg = TagAnalysis.new("LambdacSTTag")
alg.set_header(["LambdacSTTagAlg/LambdacSTTag.h"])
   .set_constant({"ECMS" => [:double, 4.600]})   # nominal; per-run beam energy handled in the fit
   .with_decay_card(decay_card_signal)

# Tag side: single tag of Λc+ (one tag_side call -> ST). Both the signal mode
# (Λc+ → Λ K+) and the reference mode (Λc+ → Λ π+) are declared as tag modes;
# the Λ → p π- sub-decay is reconstructed internally by DTagAlg.
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoLambdaK, :LambdacPtoLambdaPi   # Λc+ → Λ K+ / Λc+ → Λ π+
  t.charm 1                                          # tag the Λc+ (other side anti-Λc- inclusive)
  t.window :deltaE, min: -0.009, max: 0.012          # -0.009 < ΔE < 0.012 GeV (explicitly requested)
end

# Signal side: anti-Λc- is reconstructed inclusively -> nothing declared.
# M_BC / ΔE are stored unconditionally and windowed/fitted in ROOT.

# Kinematic fit over the derived tag participants
alg.fit do |f|
  f.constrain_four_momentum   # 4-momentum conservation against the measured CMS 4-vector
  f.chi2_cut 200              # loose BOSS-level χ²; optimal cut applied in ROOT
end

# BOSS-side procedures that are not expressible in the tag DSL
alg.note(:background_veto, "E/p < 0.9 applied to the bachelor K+ (Λc+ → Λ K+) and to the " \
                           "bachelor π+ (Λc+ → Λ π+ reference mode) to suppress the " \
                           "Λc+ → Λ e+ νe background")
alg.note(:efficiency_curve, "among multiple Λc+ tag candidates per event the candidate with " \
                            "minimum |ΔE| is retained; the M_BC distributions of both modes " \
                            "across all 13 energy points are fitted simultaneously (MC signal " \
                            "shape convolved with a Gaussian + ARGUS background, ROOT-level)")

alg.apply
root_files = alg.execute_on(data_points + incMC_points + exMC_signal + exMC_reference)