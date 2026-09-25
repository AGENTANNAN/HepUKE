# 2508.18594v1: Search for Lambda_c Sigma_c bound state H_c
# e+e- -> pi+ H_c-, H_c- -> pi- Lambda_c+ anti-Lambda_c-
# Lambda_c+ -> p K- pi+, partial recoil vs Lambda_c+ + missing anti-Lambda_c-
# Data at sqrt(s)=4918.02 MeV (BOSS 707_4914) and 4950.93 MeV (BOSS 707_4946)

### Dataset preparation ###
data_4914 = DatasetManager.real_data.find("707_4914")
data_4946 = DatasetManager.real_data.find("707_4946")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000 pi+ pi- Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.000 p+ K- pi+ PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.000 anti-p- K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_4914 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Hc_4914_exclusive_mc"
  config.related_dataset = data_4914
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

exMC_4946 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Hc_4946_exclusive_mc"
  config.related_dataset = data_4946
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection: reconstruct Lambda_c+ -> p K- pi+ with vertex fit ###
alg = Algorithm.new("HcSearch_LambdacTag")
alg.set_header(["HcSearch_LambdacTagAlg/HcSearch_LambdacTag.h"])
   .set_constant({"ECMS" => [:double, 4.91802]})

event_selection = Selection.new
event_selection.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
end
.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :kaon, against: [:pion, :proton]
  identify :pion, against: [:kaon, :proton]
  nprp ">=1"
  nkm ">=1"
  npip ">=1"
end
# Lambda_c+ reconstruction: vertex fit p K- pi+, chi2 < 200, mass window
.kinematic_fit([:prp, :km, :pip]) do
  vertex_fit([0, 1, 2])
  invariant_mass_of(:prp, :km, :pip).between(2.270, 2.300)
  chi2_cut 200
  nominal
end
.note(:partial_reconstruction,
  "Analysis uses partial reconstruction: RM(Lambda_c+) = sqrt((P_cms - P_Lambda_c+)^2). " \
  "Two tag types: 'Lambda_c Tag' (Lambda_c+ only, use RM(Lambda_c+) + M(Lambda_c+) - m(Lambda_c+)) " \
  "and 'Lambda_c pi Tag' (Lambda_c+ pi+, use 2D fit on RM(Lambda_c+pi) + M(Lambda_c+) - m(Lambda_c+) vs RM(pi)). " \
  "Remaining pi+ and pi- also vertex-fit with Lambda_c+ (chi2 < 200, best by min chi2).")
.note(:signal_extraction,
  "Unbinned maximum-likelihood fit on RM distributions (ROOT-level). " \
  "Signal MC for 15 H_c mass-width combinations (4715-4735 MeV/c^2, 5/10/20 MeV width). " \
  "Background: qqbar from Lambda_c sidebands [2.190,2.250]&[2.320,2.380]; " \
  "e+e- -> Sigma_c Sigma_c, Lambda_c Sigma_c pi, Lambda_c Lambda_c(2595), Lambda_c Lambda_c(2625).")
.note(:upper_limit,
  "90% C.L. Bayesian upper limits on sigma(e+e- -> pi+ H_c- + c.c.) x B(H_c- -> pi- Lambda_c+ anti-Lambda_c-). " \
  "Multiplicative systematics incorporated via likelihood convolution. " \
  "Additive systematics: most conservative limit across fit variations.")

alg.with_decay_card(decay_card).apply(event_selection)
alg.execute_on([data_4914, data_4946, incMC_4914, incMC_4946, exMC_4914, exMC_4946])