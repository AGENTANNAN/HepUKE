# Dataset description — energy scan e+e- -> omega eta' at 2.125 GeV (max significance)
scan_data = DatasetManager.real_data.find("713_Rscan_2125")
scan_incMC = DatasetManager.inclusive_mc.find("713_Rscan_2125")

# Decay card: e+e- -> omega eta', omega -> pi+ pi- pi0, eta' -> gamma pi+ pi-
# Uses KKMC + psi(4260) top mother (omega eta' not a predefined ConExc mode)
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0000 omega etap PHSP;
  Enddecay
  Decay omega
  1.0000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay
  Decay etap
  1.0000 gamma pi+ pi- PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "ee_to_omega_etap_4pi3gam"
  config.related_dataset = scan_data
  config.events = 500000
  config.decay_card = decay_card
  config.cross_section = :default
end

# Algorithm
alg_name = "EEToOmegaEtap"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
my_Algorithm.set_constant({"ECMS" => [:double, 2.125]})

# Event selection: omega -> pi+ pi- pi0, eta' -> gamma pi+ pi-
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
  nGam ">=3"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :pion, against: [:kaon, :proton]
  npip "==2"
  npim "==2"
}
# Kalman fit: reconstruct pi0 from photon pairs
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0 ">=1"
}
# Nominal 4C kinematic fit: 4 pions + pi0 + remaining photon constrained to CMS energy
.kinematic_fit([:pip, :pim, :pip, :pim, :pi0, :gamma]) {
  nominal
  constrain_four_momentum
  chi2_cut 60
}
# Competing hypothesis fit: 4 pions + 1 gamma (2 photons missing, pi0 not reconstructed)
.kinematic_fit([:pip, :pim, :pip, :pim, :gamma]) {
  constrain_four_momentum
}
# Competing hypothesis fit: 4 pions + 2 gammas (1 photon missing, pi0 not reconstructed)
.kinematic_fit([:pip, :pim, :pip, :pim, :gamma, :gamma]) {
  constrain_four_momentum
}

my_Algorithm
  .note(:conexc_mc_generation, "Signal MC generated with ConExc for ISR modelling and P-wave amplitude; KKMC card used here for reconstruction-only final state since omega eta' is not a predefined ConExc mode")
  .note(:chi2_comparison, "Require chi2_4C_nominal < chi2_4C_5particle AND chi2_4C_nominal < chi2_4C_6particle to veto e+e- -> 2(pi+pi-)pi0 and e+e- -> 2(pi+pi-pi0) contamination")
  .note(:pi0_selection, "Photon pair with smallest |M(gamma gamma) - M_pi0| chosen as pi0 candidate (post-fit)")
  .note(:omega_selection, "pi+pi-pi0 combination with smallest |M(pi+pi-pi0) - M_omega| chosen; events retained if |M(pi+pi-pi0) - M_omega| < 0.029 GeV/c^2 (about 3 sigma)")
  .note(:etap_selection, "Remaining gamma + remaining pi+pi-; signal region |M(gamma pi+pi-) - M_etap| < 0.025 GeV/c^2; sideband 0.05 < |M-M_etap| < 0.10 GeV/c^2")
  .note(:energy_scan, "Analysis performed at 22 energy points from sqrt(s)=2.000 to 3.080 GeV; Born cross sections extracted with ISR and VP correction factors")
  .note(:isr_correction, "ISR correction factor (1+delta) and efficiency epsilon determined iteratively from cross section line shape; VP factor |1-Pi|^2 from Ref.[25]")
  .with_decay_card(decay_card)
  .apply(event_selection)

root_files = my_Algorithm.execute_on([scan_data, scan_incMC, exMC_signal])