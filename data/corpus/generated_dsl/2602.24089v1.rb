# =============================================================================
# BOSS DSL — Double-tag study of inclusive Lambda_c+ -> Lambda X at BESIII
#   tag side   : anti-Lambda_c- in 11 hadronic single-tag (ST) modes
#   signal side: inclusive Lambda_c+ -> Lambda X,  Lambda -> p pi-  (X massless)
#   data       : 4.5 fb^-1 at 4600, 4610, 4620, 4640, 4660, 4680, 4700 MeV
# =============================================================================

### ---------------------------- Datasets ----------------------------------- ###
data_4600 = DatasetManager.real_data.find("703_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")

data_points = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]
incMC_points = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

### --------------------------- Decay card --------------------------------- ###
# Representative tag/signal card for the exclusive MC generation:
#   tag side    : anti-Lambda_c- -> anti-p- K+ pi-   (a representative ST mode)
#   signal side : Lambda_c+ -> Lambda X, with X modelled by a massless
#                 placeholder (a photon) and Lambda -> p+ pi-
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.0 anti-Lambda_c- Lambda_c+ PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0 anti-p- K+ pi- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0 Lambda gamma PHSP;
    Enddecay

    Decay Lambda
    1.0 p+ pi- HypWK;
    Enddecay

    End
DECAYCARD

### ------------------- Exclusive MC (500k events / point) ------------------ ###
# Same signal card run over each of the seven c.m. energies -> one ExclusiveMC
# per energy point, sharing cross section and event count.
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "incl_Lambdac2LambdaX_ST"   # auto-suffixed per dataset
  config.events        = 500_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

### ------------------------ Tag analysis (BOSS) --------------------------- ###
alg_name = "Lambdac2LambdaXDT"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.682] })   # sqrt(s) = 4.682 GeV
   .with_decay_card(decay_card)
   # ST tag-side selection follows the BESIII procedure of PRD 106, 072002;
   # no additional explicit PID / photon / generic track-quality cuts are imposed.
   .note(:tag_selection,
         "Single-tag side is reconstructed with the BESIII DTagAlg procedure " \
         "of Phys. Rev. D 106, 072002; no further explicit PID, photon, or " \
         "generic track-quality cuts are encoded at BOSS level. mBC/deltaE " \
         "windows are stored unconditionally and applied in the ROOT analysis.")
   # MC truth matching on both sides (unmatched events used as MC background).
   .note(:truth_matching,
         "In MC the event is truth-matched when theta(p) < 15 deg and " \
         "theta(pi-) < 25 deg for the signal-side Lambda daughters and " \
         "theta < 10 deg for the other tag-side particles; events failing " \
         "these criteria are classified as unmatched.")

# One tag side = single tag of the anti-Lambda_c- (charm = -1).
# The eleven hadronic ST modes (charge-conjugated Lambda_c+ channel names).
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP,      # p K_S0
          :LambdacPtoKPiP,     # p K- pi+
          :LambdacPtoKsPi0P,   # p K_S0 pi0
          :LambdacPtoKsPiPiP,  # p K_S0 pi+ pi-
          :LambdacPtoKPiPi0P,  # p K- pi+ pi0
          :LambdacPtoLPi,      # Lambda pi+
          :LambdacPtoLPiPi0,   # Lambda pi+ pi0
          :LambdacPtoLPiPiPi,  # Lambda pi+ pi- pi+
          :LambdacPtoS0Pi,     # Sigma0 pi+
          :LambdacPtoSPi0,     # Sigma+ pi0
          :LambdacPtoSPiPi     # Sigma+ pi+ pi-
  t.charm(-1)   # the anti-Lambda_c- carries the tag
end

# Signal side: the inclusive Lambda_c+ -> Lambda X is built from the tracks the
# tag did not use; Lambda -> p pi-, and X is a massless missing placeholder.
alg.signal_side do |s|
  s.charged(prp: 1, pim: 1, at_least: true)   # at least one proton and one pi-
  s.missing :X0, mass: nil                     # X treated as a massless placeholder
end

# 4C kinematic fit: constrain the total four-momentum to sqrt(s), chi2 < 200.
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.apply   # tag specification: takes NO Selection argument

### --------------------------- Execution ---------------------------------- ###
root_files = alg.execute_on(data_points + incMC_points + exMCs)