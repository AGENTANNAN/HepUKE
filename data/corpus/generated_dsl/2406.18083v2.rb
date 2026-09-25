### Dataset preparation — BOSS 706, threshold e+e- -> Lambda_c+ anti-Lambda_c- ###

# Six real-data points across 4.600-4.699 GeV and the corresponding inclusive MC
data_points = [
  DatasetManager.real_data.find("706_4610"),   # 4.61186 GeV
  DatasetManager.real_data.find("706_4620"),   # 4.62800 GeV
  DatasetManager.real_data.find("706_4640"),   # 4.64091 GeV
  DatasetManager.real_data.find("706_4660"),   # 4.66124 GeV
  DatasetManager.real_data.find("706_4680"),   # 4.68192 GeV
  DatasetManager.real_data.find("706_4700"),   # 4.69882 GeV
]

inclusive_mcs = [
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
]

# Decay card for the exclusive signal MC (EvtGen syntax)
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 p+ K_L0 PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K_S0 PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# 200k-event exclusive MC, one per energy point, shared by the three signal modes
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lamc_lamcbar"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

datasets = data_points + inclusive_mcs + exMCs_signal

# Twelve Cabibbo-favored tag modes used to tag the anti-Lambda_c-
tag_modes = [
  :LambdacPtoKsP,         # anti-p K_S0
  :LambdacPtoKPiP,        # anti-p K- pi+
  :LambdacPtoKsPi0P,      # anti-p K_S0 pi0
  :LambdacPtoKsPiPiP,     # anti-p K_S0 pi+ pi-
  :LambdacPtoKPiPi0P,     # anti-p K- pi+ pi0
  :LambdacPtoPiPiP,       # anti-p pi+ pi-
  :LambdacPtoLambdaPi,    # anti-Lambda pi+
  :LambdacPtoLambdaPiPi0, # anti-Lambda pi+ pi0
  :LambdacPtoLambdaPiPiPi,# anti-Lambda pi+ pi- pi+
  :LambdacPtoSigma0Pi,    # anti-Sigma0 pi+
  :LambdacPtoSigmaPi0,    # anti-Sigma- pi0
  :LambdacPtoSigmaPiPi,   # anti-Sigma- pi+ pi-
]

tag_veto_note = "tag-side sub-decay vetoes applied in ROOT: " \
  "M(p pi-) outside [1.11, 1.12] GeV/c2 (Lambda -> p pi-), " \
  "M(pi+ pi-) outside [0.48, 0.52] GeV/c2 (K_S0 -> pi+ pi-), " \
  "M(p pi0) outside [1.17, 1.20] GeV/c2 (Sigma+ -> p pi0)"

root_note_prefix = "applied in the later ROOT analysis: " \
  "signal Lambda_c+ mass constraint, pi0 1C fit and the optimized chi2 cut "

### Signal mode I: Lambda_c+ -> p K_L0, K_L0 missing ###
alg_I = TagAnalysis.new("LcDTagPKL")
alg_I.set_header(["LcDTagPKLAlg/LcDTagPKL.h"])
     .set_constant({"ECMS" => [:double, 4.68192]})   # nominal scan energy; ecms_lab read per run
     .note(:background_veto, tag_veto_note)
     .note(:root_stage_constraints, root_note_prefix + "(< 60) for p K_L0")
     .note(:track_quality_and_pid, "charged-track quality cuts and PID are delegated to the analysis headers")
     .with_decay_card(decay_card_signal)

alg_I.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm(-1)             # tag the anti-Lambda_c-
end

alg_I.signal_side do |s|
  s.charged(prp: 1)       # single proton
  s.missing :K_L0         # K_L0 left unreconstructed
  s.require_charge(1)     # signal-side net charge +1
end

alg_I.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg_I.apply

### Signal mode II: Lambda_c+ -> p K_L0 pi+ pi-, K_L0 missing ###
alg_II = TagAnalysis.new("LcDTagPKLPiPi")
alg_II.set_header(["LcDTagPKLPiPiAlg/LcDTagPKLPiPi.h"])
      .set_constant({"ECMS" => [:double, 4.68192]})
      .note(:background_veto, tag_veto_note)
      .note(:recoil_mass_cut, "M_recoil(p) > 1.0 GeV/c2 applied in ROOT for p K_L0 pi+ pi-")
      .note(:root_stage_constraints, root_note_prefix + "(< 25) for p K_L0 pi+ pi-")
      .note(:track_quality_and_pid, "charged-track quality cuts and PID are delegated to the analysis headers")
      .with_decay_card(decay_card_signal)

alg_II.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm(-1)
end

alg_II.signal_side do |s|
  s.charged(prp: 1, pip: 1, pim: 1)   # proton + pi+ pi-
  s.missing :K_L0
  s.require_charge(1)
end

alg_II.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg_II.apply

### Signal mode III: Lambda_c+ -> p K_L0 pi0 (pi0 -> gamma gamma), K_L0 missing ###
alg_III = TagAnalysis.new("LcDTagPKLPi0")
alg_III.set_header(["LcDTagPKLPi0Alg/LcDTagPKLPi0.h"])
       .set_constant({"ECMS" => [:double, 4.68192]})
       .note(:background_veto, tag_veto_note)
       .note(:recoil_mass_cut, "M_recoil(p) > 0.65 GeV/c2 applied in ROOT for p K_L0 pi0")
       .note(:pi0_window, "pi0 candidate window 0.115-0.150 GeV/c2 applied in ROOT")
       .note(:photon_isolation, "signal photon isolation > 20 deg from the proton track; " \
              "the two pi0 photons are also required to have opening angle > 10 deg")
       .note(:root_stage_constraints, root_note_prefix + "(< 20) for p K_L0 pi0")
       .note(:track_quality_and_pid, "charged-track quality cuts and PID are delegated to the analysis headers")
       .with_decay_card(decay_card_signal)

alg_III.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm(-1)
end

alg_III.signal_side do |s|
  s.charged(prp: 1)          # proton
  s.photons 2                # two photons from pi0
  s.min_photon_angle 20.0    # photon isolation from the proton track
  s.missing :K_L0
  s.require_charge(1)
end

alg_III.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # pi0 mode only
  f.chi2_cut 200
end

alg_III.apply

# Run the three signal-mode analyses over the full scan (data + inclusive MC + signal MC)
root_files_I   = alg_I.execute_on(datasets)
root_files_II  = alg_II.execute_on(datasets)
root_files_III = alg_III.execute_on(datasets)