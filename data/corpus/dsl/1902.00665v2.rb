# ISR measurement of e+e- -> p pbar cross section using undetected photon
# Paper: 1902.00665v2, multi-energy: 3.773, 4.008, 4.226, 4.258, 4.358, 4.416, 4.600 GeV

# Multi-energy dataset (7 energy points)
data_3773 = DatasetManager.real_data.find("712_3773")
data_4009 = DatasetManager.real_data.find("703_4009")
data_4230 = DatasetManager.real_data.find("703_4230")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4360 = DatasetManager.real_data.find("703_4360")
data_4420 = DatasetManager.real_data.find("703_4420")
data_4600 = DatasetManager.real_data.find("703_4600")

all_datasets = [data_3773, data_4009, data_4230, data_4260, data_4360, data_4420, data_4600]

# ISR decay card: e+e- -> p pbar gamma_ISR (undetected photon)
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.000 p+ anti-p- gamma PHSP;
  Enddecay
  End
DECAYCARD

exMCs = DatasetManager.create_exclusive_mc_for(all_datasets) do |config|
  config.sample_name = "ee_ppbar_gamma_ISR_exMC"
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

algorithm = Algorithm.new("EEtoPPbarISR")
algorithm
  .set_header(["EEtoPPbarISRAlg/EEtoPPbarISR.h"])

selection = Selection.new
selection.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==1"
  nChrn "==1"
  nNet "==0"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprp "==1"
  nprm "==1"
end

algorithm
  .note(:isr_technique, "ISR analysis with undetected photon: signal identified via missing momentum kinematics rather than photon detection; no explicit photon selection is performed")
  .note(:pid_emc_cut, "E_EMC/p_rec < 0.5 for positively charged track to suppress e+ background")
  .note(:theta_miss_cut, "polar angle of missing momentum required to be < 0.125 or > (pi - 0.125) rad to select small-angle ISR photon events")
  .note(:mmiss2_cut, "missing mass squared cut: [-0.1, 0.2] GeV^2/c^4 for sqrt(s)>4 GeV, [-0.02, 0.10] GeV^2/c^4 for sqrt(s)=3.773 GeV; suppresses e+e-gamma, ppbar pi0 gamma, and two-photon backgrounds")
  .note(:cos_theta_pp_cut, "proton/antiproton polar angle in pp CM system |cos(theta_p,pp^CM)| < 0.75")
  .note(:background_subtraction, "background from ppbar pi0 and two-photon channels estimated by sideband method and subtracted before cross section calculation")
  .note(:jpsi_psip_subtraction, "J/psi->ppbar and psi(3686)->ppbar resonances fitted and subtracted from M_ppbar distribution")
  .note(:radiative_correction, "signal MC generated with phokhara NLO; detection efficiency and radiative correction factor (1+delta) determined in each M_ppbar interval")
  .note(:no_kinematic_fit, "no BOSS-level kinematic fit is used in this analysis; event selection is based solely on track selection, PID, E/p cut, and missing kinematics cuts (theta_miss and M_miss^2)")
  .with_decay_card(decay_card)
  .apply(selection)

algorithm.execute_on(all_datasets + exMCs)