# Paper 2107.06704v1: Lambda -> p mu- nu_mu_bar branching fraction at J/psi
# Double-tag technique: ST Lambda(-> p-bar pi+) + signal Lambda(-> p mu- nu_mu_bar)
# Single energy at sqrt(s)=3.097 GeV (J/psi, 10 billion events)

decay_card = <<~DECAYCARD
  Decay J/psi
  1.0000 Lambda0 anti-Lambda0 PHSP;
  Enddecay
  Decay Lambda0
  0.50 p+ pi- HypWK;
  0.50 p+ mu- anti-nu_mu PHSP;
  Enddecay
  Decay anti-Lambda0
  0.50 anti-p- pi+ HypWK;
  0.50 anti-p- mu+ nu_mu PHSP;
  Enddecay
  End
DECAYCARD

data_jpsi = DatasetManager.real_data.find("708_3097")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")

sigMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name   = "sig_lambda_pmunu"
  config.related_dataset = data_jpsi
  config.events        = 100_000
  config.decay_card    = decay_card
  config.cross_section = :default
end

algorithm = Algorithm.new("LambdaPMuNu")
algorithm.set_header(["LambdaPMuNuAlg/LambdaPMuNu.h"])
algorithm.set_constant({"ECMS" => [:double, 3.097]})

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
  identify :muon, against: [:electron, :pion, :kaon, :proton]
  nprp ">=1"
  nprm ">=1"
  nmup ">=1"
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
  constrain_four_momentum
end

algorithm
  .note(:double_tag_technique, "DT method: ST Lambda->p-bar pi+ selected by minimum |DeltaE_tag|; DeltaE_tag in [-17,13] MeV; ST yield from M_BC fit in [1.089,1.143] GeV/c^2")
  .note(:st_vertex_fit, "ST Lambda vertex fit chi^2 < 100; decay length > 2 sigma from IP; at least one Lambda required; no PID for ST side")
  .note(:dt_selection, "DT from remaining tracks: total N_track=4; muon PID L_mu>0.001 and L_mu>L_e; proton assumed from other track; missing neutrino carries U_miss")
  .note(:signal_variable, "U_miss = E_miss - c|p_miss| peaks at 0 for signal; E_miss = E_beam - E_p - E_mu; p_miss from constrained Lambda momentum")
  .note(:background_veto, "4C kinematic fit chi^2 > 20 to suppress Lambda->p pi- (this 4C fit is a veto: signal events FAIL it due to missing neutrino); M_recoil(Lambda_bar p) > 0.170 GeV/c^2; M_pmu(4C) in [1.075,1.100] GeV/c^2")
  .note(:signal_extraction, "DT yield from unbinned maximum likelihood fit to U_miss distribution; BF = (N_DT/eps_DT) / (N_ST/eps_ST)")
  .note(:cp_asymmetry, "Separate BF measurements for Lambda->p mu- nu_mu_bar and Lambda_bar->p-bar mu+ nu_mu; CP asymmetry A_CP from difference of charge-conjugated BFs")
  .with_decay_card(decay_card)
  .apply(event_selection)

algorithm.execute_on([data_jpsi, incMC_jpsi, sigMC])