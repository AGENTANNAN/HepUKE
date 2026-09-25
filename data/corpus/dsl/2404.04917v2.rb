# Dataset description
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for signal: psi(3686) -> gamma chi_c0 -> gamma 2(pi+pi-)
# Same final state shared by eta_c(2S) and chi_cJ (J=0,1,2) — single Algorithm per Rule T1
decay_card = <<~DECAYCARD
  Decay psi(3686)
  1.0000 gamma chi_c0 PHSP;
  Enddecay
  Decay chi_c0
  1.0000 pi+ pi- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive MC for signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_gamma_chic_4pi"
  config.related_dataset = psip_data
  config.events = 1000000
  config.decay_card = decay_card
  config.cross_section = :default
end

# Algorithm
alg_name = "PsiPToGammaXTo4Pi"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
my_Algorithm.set_constant({"ECMS" => [:double, 3.686]})

# Event selection: gamma + pi+ pi- pi+ pi- (exactly 4 charged tracks, >=1 photon)
event_selection = Selection.new
event_selection.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==2"
  nNet "==0"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=1"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :pion, against: [:kaon, :proton]
  npip "==2"
  npim "==2"
}
# Nominal 4C kinematic fit: photon + 4 pions constrained to CMS energy
# Vertex fit on the 4 charged pions (common vertex constraint)
.kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) {
  nominal
  vertex_fit([1, 2, 3, 4])
  constrain_four_momentum
  chi2_cut 40
}

my_Algorithm
  .note(:pid_chi2_sum, "Events retained if sum_i chi2_PID(pi,i) < any other hadron hypothesis assignment over 4 tracks; equivalent to all 4 tracks best-identified as pions")
  .note(:jpsi_veto, "Recoil mass of all pi+pi- pairs required outside J/psi mass region: M_rec(pi+pi-) < 3.0 GeV/c^2 or > 3.2 GeV/c^2, suppressing psi(3686) -> pi+pi- J/psi background")
  .note(:eta_jpsi_veto, "Recoil mass of all gamma pi+pi- combinations required outside J/psi mass region: M_rec(gamma pi+pi-) < 3.0 GeV/c^2 or > 3.2 GeV/c^2, suppressing psi(3686) -> eta J/psi with eta -> gamma pi+pi-")
  .note(:gamma_conversion_veto, "cos(theta) of all pi+pi- pairs required in [-0.999, 0.988] to veto gamma -> e+e- conversions where both e+/e- misidentified as pions")
  .note(:three_c_fit, "3C kinematic fit (omitting photon energy measurement) used for signal extraction; M_{2(pi+pi-)}^{3C} used as fit variable to separate psi(3686)->2(pi+pi-) background from eta_c(2S) signal")
  .note(:fsr_correction, "FSR correction factor f_FSR = 2.00 +/- 0.02 applied to MC; psi(3686)->gamma_FSR 2(pi+pi-) backgrounds scaled accordingly")
  .note(:helix_correction, "Helix parameter correction applied to charged tracks in MC to improve data-MC consistency in kinematic fit")
  .with_decay_card(decay_card)
  .apply(event_selection)

root_files = my_Algorithm.execute_on([psip_data, psip_incMC, exMC_signal])