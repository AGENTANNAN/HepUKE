# Dataset description — energy scan e+e- -> eta hc
ds_data = DatasetManager.real_data.find("703_4180")
ds_incMC = DatasetManager.inclusive_mc.find("703_4180")

# Decay card: e+e- -> eta hc, eta -> gamma gamma, hc -> gamma eta_c, eta_c -> p anti-p-
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta hc PHSP;
  Enddecay
  Decay hc
  1.0000 gamma eta_c VSP_PWAVE;
  Enddecay
  Decay eta_c
  1.0000 p+ anti-p- HELAMP 0 0 1 0 1 0;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "ee_to_eta_hc_etac_ppbar"
  config.related_dataset = ds_data
  config.events = 500000
  config.decay_card = decay_card
  config.cross_section = :default
end

# Algorithm — eta_c → p pbar (one of sixteen hadronic modes)
alg_name = "EEToEtaHcEtacToPPbar"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
my_Algorithm.set_constant({"ECMS" => [:double, 4.178]})

# Event selection: eta -> gamma gamma, hc -> gamma eta_c, eta_c -> p pbar
event_selection = Selection.new
event_selection.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==1"
  nChrn "==1"
  nNet "==0"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=3"  # 2 from eta -> gamma gamma + 1 transition photon from hc -> gamma eta_c
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprp "==1"
  nprm "==1"
}
# Kalman fit: reconstruct eta -> gamma gamma (mass constraint)
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 25
  neta "==1"
}
# Nominal 4C kinematic fit: eta + gamma + p + anti-p constrained to CMS energy
.kinematic_fit([:eta, :gamma, :prp, :prm]) {
  nominal
  constrain_four_momentum
  chi2_cut 25
}

my_Algorithm
  .note(:multi_mode_simultaneous_fit, "16 eta_c hadronic decay modes analyzed simultaneously with weighted simultaneous fit to eta recoil mass; weights f_i = epsilon_i * B_i / sum(epsilon_i * B_i)")
  .note(:chi2_minimization, "Best PID assignment, pi0/eta/K_S0 candidates chosen by minimizing chi2 = chi2_4C + sum(chi2_1C) + sum(chi2_PID) + sum(chi2_vertex) over all combinations")
  .note(:eta_c_tag, "eta_c candidate mass closest to nominal eta_c mass selected when multiple eta candidates have recoil mass in hc signal region [3.480, 3.600] GeV/c^2")
  .note(:hc_mass_window, "hc signal window: eta recoil mass in [3.510, 3.540] GeV/c^2")
  .note(:etac_mass_window, "eta_c signal window: hadronic invariant mass in [2.944, 3.024] GeV/c^2")
  .note(:isr_correction, "ISR correction factor (1+delta) calculated iteratively using cross section line shape from threshold; VP correction factor |1+Pi|^2 applied for Born cross section extraction")
  .note(:energy_scan, "Analysis repeated at 22 energy points from sqrt(s)=4.129 to 4.600 GeV; Born cross sections extracted per energy point")
  .with_decay_card(decay_card)
  .apply(event_selection)

root_files = my_Algorithm.execute_on([ds_data, ds_incMC, exMC_signal])