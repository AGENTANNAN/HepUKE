# Paper: 2307.07316v1
# Analysis: e+e- -> Lambda_c+ Lambda_c- cross sections and electromagnetic form factors
# Energy: 12 points from 4.6119 to 4.9509 GeV (plus 4.5995 GeV for DT)
# Method: Tag-based analysis (TagAnalysis) — ST (single tag) primary;
#          DT (double tag) used for BF cancellation at ROOT level

# --- Datasets (12 energy points + 4.5995 for DT) ---
data_4612 = DatasetManager.real_data.find("706_4610")
data_4628 = DatasetManager.real_data.find("706_4620")
data_4641 = DatasetManager.real_data.find("706_4640")
data_4661 = DatasetManager.real_data.find("706_4660")
data_4682 = DatasetManager.real_data.find("706_4680")
data_4699 = DatasetManager.real_data.find("706_4700")
data_4740 = DatasetManager.real_data.find("707_4740")
data_4750 = DatasetManager.real_data.find("707_4750")
data_4781 = DatasetManager.real_data.find("707_4780")
data_4843 = DatasetManager.real_data.find("707_4840")
data_4918 = DatasetManager.real_data.find("707_4914")
data_4951 = DatasetManager.real_data.find("707_4946")
data_4600 = DatasetManager.real_data.find("703_4600")  # 4.5995 GeV for DT analysis

all_data = [data_4612, data_4628, data_4641, data_4661, data_4682, data_4699,
            data_4740, data_4750, data_4781, data_4843, data_4918, data_4951, data_4600]

# Inclusive MC at 4.682 GeV
incMC_4682 = DatasetManager.inclusive_mc.find("706_4680")

# --- Decay card for signal MC ---
decay_card = <<~DECAYCARD
  Decay anti-Lambda_c-
  1.000 p+ K- pi+ PHSP;
  Enddecay
  End
DECAYCARD

# Signal MC per energy point
sig_mc_samples = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name = "sig_Lc_Lcbar"
  config.events = 200_000
  config.decay_card = decay_card
  config.cross_section = :default
end

# --- TagAnalysis: ST (single tag) ---
# Primary measurement uses ST: tag one Lambda_c+ -> p K- pi+
# Signal side is empty (other Lambda_c- decays inclusively, unreconstructed)
algorithm = TagAnalysis.new("LcLcBar")
algorithm.set_header(["LcLcBarAlg/LcLcBar.h"])
  .set_constant({ "ECMS" => [:double, 4.682] })
  .with_decay_card(decay_card)
  .note(:st_dt_hybrid, "ST method used for cross-section and form-factor measurements; DT method (both Lambda_c+ and Lambda_c- tagged) used to cancel BF systematic uncertainty via 2D simultaneous fit to (M_BC+, M_BC-)")
  .note(:dt_fit, "DT analysis uses 2D simultaneous unbinned likelihood fit on (M_BC+, M_BC-) distributions across 9 energy points; BF shared across energies; signal = MC shape convolved with 2D Gaussian; background = ARGUS x Gaussian")
  .note(:st_fit, "ST yield extracted from M_BC distribution; signal shape from MC convolved with Gaussian; background = ARGUS function")
  .note(:cross_section_formula, "sigma = N_ST / (epsilon_ST * f_ISR * f_VP * L_int * B); f_ISR from KKMC iterative procedure; f_VP from ConExc ~1.055")
  .note(:form_factors, "|G_eff| = sqrt(sigma / (sigma0/3 * (1 + kappa/2))); polar-angle distributions fit with f(cos_theta) = N0*(1+alpha*cos^2_theta); |GE/GM| from alpha; |GM| extracted separately")
  .note(:ge_gm_oscillation, "First observation of energy-dependent |GE/GM| ratio for Lambda_c+; oscillation frequency ~32 GeV^-1, ~3.5x greater than proton")
  .note(:no_y4630, "No Y(4630) resonant structure observed around 4.63 GeV, in contrast to Belle ISR measurement")
  .note(:deltaE_window, "Asymmetric DeltaE window (-34, 20) MeV applied at ROOT level (not expressible in tag_side window which uses symmetric abs)")
  .note(:cos_theta_binning, "ST sample divided into 20 cos_theta bins for polar-angle distribution; ISR correction applied per bin")

# Tag side: Lambda_c+ -> p K- pi+
algorithm.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP
end

# Signal side: empty — other Lambda_c- decays inclusively (ST method)
# No charged tracks or photons required on signal side
algorithm.signal_side do |s|
  s.photons 0
end

# Fit: 4C kinematic fit
algorithm.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

algorithm.apply
algorithm.execute_on(all_data + [incMC_4682] + sig_mc_samples)