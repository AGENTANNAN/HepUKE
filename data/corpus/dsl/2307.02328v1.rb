# Paper: 2307.02328v1
# Analysis: e+e- -> p K- Lambda_bar + c.c. cross sections
# Energy: 37 points from 4.009 GeV to 4.951 GeV (R-scan)
# Method: Ordinary analysis with secondary vertex fit for Lambda

# --- Datasets ---
# 37 energy points across BOSS versions 703, 705, 706, 707
data_4009 = DatasetManager.real_data.find("703_4009")
data_4130 = DatasetManager.real_data.find("705_4130")
data_4160 = DatasetManager.real_data.find("705_4160")
data_4180 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")
data_4237 = DatasetManager.real_data.find("703_4237")
data_4246 = DatasetManager.real_data.find("703_4246")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4270 = DatasetManager.real_data.find("703_4270")
data_4280 = DatasetManager.real_data.find("703_4280")
data_4290 = DatasetManager.real_data.find("705_4290")
data_4315 = DatasetManager.real_data.find("705_4315")
data_4340 = DatasetManager.real_data.find("705_4340")
data_4360 = DatasetManager.real_data.find("703_4360")
data_4380 = DatasetManager.real_data.find("705_4380")
data_4400 = DatasetManager.real_data.find("705_4400")
data_4420 = DatasetManager.real_data.find("703_4420")
data_4440 = DatasetManager.real_data.find("705_4440")
data_4470 = DatasetManager.real_data.find("703_4470")
data_4530 = DatasetManager.real_data.find("703_4530")
data_4600 = DatasetManager.real_data.find("703_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")
data_4740 = DatasetManager.real_data.find("707_4740")
data_4750 = DatasetManager.real_data.find("707_4750")
data_4780 = DatasetManager.real_data.find("707_4780")
data_4840 = DatasetManager.real_data.find("707_4840")
data_4914 = DatasetManager.real_data.find("707_4914")
data_4946 = DatasetManager.real_data.find("707_4946")

all_data = [data_4009, data_4130, data_4160, data_4180, data_4190, data_4200,
            data_4210, data_4220, data_4230, data_4237, data_4246, data_4260,
            data_4270, data_4280, data_4290, data_4315, data_4340, data_4360,
            data_4380, data_4400, data_4420, data_4440, data_4470, data_4530,
            data_4600, data_4610, data_4620, data_4640, data_4660, data_4680,
            data_4700, data_4740, data_4750, data_4780, data_4840, data_4914,
            data_4946]

# Inclusive MC at 4.178 GeV (40x lumi of data)
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")

# --- Decay card for signal MC ---
# Note: The actual analysis uses a custom amplitude model with intermediate
# resonances (X(2085), K*2(1980), K*4(2045), K2(2250), Lambda(1520), etc.)
# for physics-accurate generation. The PHSP card below is a placeholder.
decay_card = <<~DECAYCARD
  Decay vpho
  1.000 p+ K- anti-Lambda- PHSP;
  Enddecay
  Decay anti-Lambda-
  1.000 anti-p- pi+ PHSP;
  Enddecay
  End
DECAYCARD

# Signal MC over all energy points
sig_mc_samples = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name = "sig_pK_Lambdabar"
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

# --- Algorithm ---
algorithm = Algorithm.new("PKLambdaBar")
algorithm.set_header(["PKLambdaBarAlg/PKLambdaBar.h"])
  .set_constant({ "ECMS" => [:double, 4.178] })
  .note(:amplitude_model, "Signal MC in the paper uses a custom amplitude model with intermediate resonances (X(2085), K* states, Lambda* states, N* states); PHSP card is a placeholder. Detection efficiency is model-dependent.")
  .note(:no_pid, "No particle identification applied for charged tracks due to high momenta; track assignment resolved via Lambda vertex fit and 4C kinematic fit")
  .note(:isr_correction, "ISR correction factor (1+delta) obtained iteratively using signal MC with a 4th-order polynomial line shape times (1-exp(-DeltaM/p0)); vacuum polarization factor 1/|1-Pi|^2 from ConExc")
  .note(:conexc, "Inclusive MC open charm production generated with ConExc; vacuum polarization factor obtained from ConExc. ConExc mode 2 (Lambda Lambda_bar) is closest predefined mode but does not match the p K Lambda_bar final state exactly.")
  .note(:lambda_mass_window, "Lambda mass window [1.10, 1.13] GeV/c^2 applied in for_each block")
  .note(:decay_length, "Lambda decay length > 2x vertex resolution required; not expressible in DSL — applied via for_each note")
  .note(:single_lambda, "Exactly one Lambda candidate required per event; best candidate selected")
  .note(:cos_theta_k_veto, "|cos_theta_K| < 0.83 veto applied in for_each block to suppress beam-induced kaon backgrounds")
  .note(:track_dz_dr, "Additional track quality cuts on non-Lambda tracks: |dz| < 10 cm and |dr| < 1 cm; tighter than the initial 20 cm requirement")
  .note(:signal_extraction, "Signal yield extracted from unbinned ML fit to M(p pi-) distribution; signal shape = MC shape convolved with Gaussian; background = linear; performed at ROOT level")
  .note(:cross_section, "Born cross section sigma_B = N_sig / (L_int * B * epsilon * (1+delta) * 1/|1-Pi|^2); dressed cross section sigma_D = sigma_B * 1/|1-Pi|^2")

selection = Selection.new
selection.select_track do
  cos_theta 0.93
  Vz 20.0
  nChrp "==2"
  nChrn "==2"
  nNet "==0"
end
.pid(method: :probability) do
  prob_cut 0.0
  identify :proton, against: [:pion, :kaon]
  identify :pion, against: [:proton, :kaon]
  identify :kaon, against: [:proton, :pion]
  nprp ">=1"
  nprm ">=1"
  npip ">=1"
  npim ">=1"
  nkp ">=0"
  nkm ">=0"
end
.secondary_vertex_fit([:prp, :pim]) do
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
.for_each(:Lambda) do
  where { mass >= 1.10 && mass <= 1.13 }
  best
end
.for_each(:km) do
  where { cos_theta.abs <= 0.83 }
end
.kinematic_fit([:Lambda, :prp, :km]) do
  nominal
  constrain_four_momentum
  chi2_cut 100
end

algorithm.with_decay_card(decay_card).apply(selection)
algorithm.execute_on(all_data + [incMC_4180] + sig_mc_samples)