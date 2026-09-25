# Paper 2108.02410v2: e+e- -> Lambda anti-Lambda cross section at sqrt(s)=3.51-4.60 GeV
# Ordinary analysis with secondary vertex fit for Lambda / anti-Lambda

decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda0 anti-Lambda0 PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay
  End
DECAYCARD

# Multi-energy scan from 3.51 to 4.60 GeV
data_scan = DatasetManager.real_data.where(cms_energy: {value: 3500..4610})
incMC_scan = DatasetManager.inclusive_mc.where(cms_energy: {value: 3500..4610})

sigMC = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "sig_lambda_lambdabar"
  config.events        = 100_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

algorithm = Algorithm.new("LambdaLambdabar")
algorithm.set_header(["LambdaLambdabarAlg/LambdaLambdabar.h"])

event_selection = Selection.new
event_selection.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==2"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprp ">=1"
  nprm ">=1"
end
.remove([:prp <= :chrgp, :prm <= :chrgn])
.assign({chrgp: :pip, chrgn: :pim})
.secondary_vertex_fit([:prp, :pim]) do
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
.secondary_vertex_fit([:prm, :pip]) do
  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end
.kinematic_fit([:Lambda, :Lambda_bar]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

algorithm
  .note(:secondary_vertex_chi2_cut, "secondary vertex fit chi^2 < 500; Lambda decay length > 0 required")
  .note(:lambda_mass_window, "Lambda mass window: |M(ppi) - m_Lambda| < 5 MeV/c^2 in signal region; sideband subtraction for background estimation")
  .note(:signal_extraction, "signal yield from sideband subtraction: N_obs = N_S - (1/4)*sum(N_Bi); cross section extracted from dressed cross section fit")
  .with_decay_card(decay_card)
  .apply(event_selection)

algorithm.execute_on(data_scan + incMC_scan + sigMC)