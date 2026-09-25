# Determination of Spin-parity Quantum Numbers for X(2085) in e+e- -> p K- Lambda + c.c.
# 6 energy points: 4.008, 4.178, 4.226, 4.258, 4.416, 4.682 GeV

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# Data at 6 energy points
data_4009 = DatasetManager.real_data.find("703_4009")
data_4180 = DatasetManager.real_data.find("703_4180")
data_4230 = DatasetManager.real_data.find("703_4230")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4420 = DatasetManager.real_data.find("703_4420")
data_4680 = DatasetManager.real_data.find("706_4680")
data_points = [data_4009, data_4180, data_4230, data_4260, data_4420, data_4680]

# Inclusive MC
incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4420 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_points = [incMC_4009, incMC_4180, incMC_4230, incMC_4260, incMC_4420, incMC_4680]

# Decay card: e+e- -> p K- Lambda (c.c. implied); KKMC generator via psi(4260)
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0000 p+ K- Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- HypWK;
  Enddecay

  End
DECAYCARD

# Exclusive MC for energy scan
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_pKLambda"
  config.events        = 100_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

# Algorithm: e+e- -> p K- Lambda (charge conjugated mode implied)
alg = Algorithm.new("pKLambdaAnalysis")
alg.set_header(["pKLambdaAnalysisAlg/pKLambdaAnalysis.h"])

# Event selection
event_selection = Selection.new
event_selection.select_track do
  cos_theta 0.93    # |cos(theta)| < 0.93
  Vz 20.0           # |dz| < 20 cm for all charged tracks
  nChrp ">=2"       # at least 2 positive tracks (p from Lambda + prompt p)
  nChrn ">=2"       # at least 2 negative tracks (pi- from Lambda + K-)
  nNet  "==0"       # net charge zero
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :kaon, against: [:pion, :proton]
  nprp ">=2"        # two protons: one from Lambda, one prompt
  nkm  ">=1"        # one K-
end
.remove([:prp <= :chrgp, :prm <= :chrgn, :kp <= :chrgp, :km <= :chrgn])
.assign({chrgp: :pip, chrgn: :pim})
# Lambda -> p pi- reconstruction via secondary vertex fit
.secondary_vertex_fit([:prp, :pim]) do
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
# 4C kinematic fit: e+e- -> p K- Lambda
.kinematic_fit([:prp, :km, :Lambda]) do
  nominal
  constrain_four_momentum
  chi2_cut 100       # chi2 < 100
end

alg.note(:lambda_mass_window, "|M(p pi-) - M(Lambda)| < 6 MeV/c2 applied to Lambda candidates before 4C fit")
   .note(:lambda_decay_length, "Lambda decay length > 2 sigma of vertex resolution required")
   .note(:lambda_vertex_chi2, "Secondary vertex fit chi2 < 100 for Lambda candidates (p pi- vertex)")
   .note(:prompt_track_cuts, "|dz| < 10 cm and |dr| < 1 cm for prompt proton and kaon not from Lambda decay (tighter than Lambda daughter cuts)")
   .note(:cos_theta_k_veto, "Events with |cos(theta_K)| > 0.83 rejected to suppress e+e- -> (gamma) e+e- Bhabha background")
   .note(:charge_conjugate, "Charge conjugated mode e+e- -> pbar K+ Lambda_bar implied throughout")
   .note(:amplitude_analysis, "Amplitude analysis with X(2085), K*(1980), K*(2045), K2(2250), Lambda(1520), Lambda(1890), Lambda(2350), N(1720), N(2570) performed in ROOT; spin-parity determined from unbinned maximum likelihood fit")
   .with_decay_card(decay_card)
   .apply(event_selection)

datasets = data_points + incMC_points + exMCs
alg.execute_on(datasets)