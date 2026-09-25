# frozen_string_literal: true

### Dataset preparation ###
# Real data and inclusive MC at the seven c.m. energies:
# 4599.53, 4611.86, 4628.00, 4640.91, 4661.24, 4681.92, 4698.82 MeV
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

data_points  = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]
incMC_points = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

# Decay card: e+e- -> Lambda_c+ anti-Lambda_c- (PHSP);
# Lambda_c+ -> K_S0 X (inclusive signal); anti-Lambda_c- -> eight hadronic tag modes;
# K_S0 -> pi+ pi-; anti-Lambda -> anti-p pi+.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 K_S0 X PHSP;
    Enddecay

    Decay anti-Lambda_c-
    0.1250 anti-p- K_S0 PHSP;
    0.1250 anti-p- K+ pi- PHSP;
    0.1250 anti-p- K_S0 pi0 PHSP;
    0.1250 anti-p- K_S0 pi- pi+ PHSP;
    0.1250 anti-p- K+ pi- pi0 PHSP;
    0.1250 anti-Lambda0 pi- PHSP;
    0.1250 anti-Lambda0 pi- pi0 PHSP;
    0.1250 anti-Lambda0 pi- pi+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# 1M-event exclusive signal MC generated per energy point (one MC per dataset)
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lambdac_KSX"
  config.events        = 1_000_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMC_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) — tag-based double-tag analysis ###
alg_name = "LambdacKSXInclusive"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.650] })

# Tag side: anti-Lambda_c- reconstructed from the pre-stored tag candidates via eight
# hadronic modes; charm -1 pins the tagged (anti-Lambda_c-) side.
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP, :LambdacPtoKPiP, :LambdacPtoKsPi0P, :LambdacPtoKsPiPiP,
          :LambdacPtoKPiPi0P, :LambdacPtoPiL, :LambdacPtoPiPi0L, :LambdacPtoPiPiPiL
  t.charm -1
end

# Signal side: K_S0 -> pi+ pi- built from the tracks not used by the tag
# (the inclusive X may contribute further unused tracks).
alg.signal_side do |s|
  s.charged(pip: 1, pim: 1, at_least: true)
end

# Kinematic fit: four-momentum constraint, K_S0 mass window M(pi+pi-) in [0.487, 0.511] GeV,
# chi2 < 200. Tag mBC / deltaE and the signal mass are stored and windowed later in ROOT
# (store-not-cut: no explicit tag-side window is declared).
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).between(0.487, 0.511)
  f.chi2_cut 200
end

# Inexpressible BOSS-side procedures preserved for the systematic/ROOT stage
alg.note(:background_veto, "mode-dependent intermediate-resonance vetos applied to the anti-Lambda_c- tag: anti-Lambda via M(anti-p pi+) outside (1.110, 1.120) GeV/c^2 for the anti-p K_S0 pi0, anti-p K_S0 pi- pi+ and anti-Sigma- pi- pi+ modes; anti-Sigma- via M(anti-p pi0) outside (1.170, 1.200) GeV/c^2 for anti-p K_S0 pi0; K_S0 via M(pi+ pi-) and M(pi0 pi0) outside (0.480, 0.520) GeV/c^2 for anti-Lambda pi- pi+ pi-")
   .note(:efficiency_curve, "signal-side K_S0 -> pi+ pi- requires a common vertex fit with chi2 < 100 and decay length L > 2 sigma_L, keeping the candidate with the largest L/sigma_L; only the M(pi+ pi-) window is expressible in the tag DSL")
   .note(:truth_matching, "signal-side K_S0 candidates truth-matched with R_dp < 0.5; inclusive K_S0 reconstruction efficiency corrected by 1.8 +/- 0.3%")
   .note(:tag_mode_coverage, "the paper's three Sigma0/Sigma- tag modes are not reconstructed here because they have no matching DTagAlg tag symbol; the mode-dependent veto that referred to the anti-Sigma- pi- pi+ mode is therefore inert")

# Validate/render the tag spec (apply takes no Selection) and run over all datasets
root_files = alg.apply.execute_on(data_points + incMC_points + exMC_signal)