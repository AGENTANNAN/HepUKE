# Paper 2004.13788v2: e+e- -> pi0 pi0 J/psi cross sections at 27 energy points (3.808-4.600 GeV)
# J/psi -> l+ l- (e/mu), pi0 -> gamma gamma, 6C kinematic fit
# Lepton ID via E/p; Born cross section measurement with ISR correction

decay_card = <<~DECAYCARD
Decay vpho
1.0 pi0 pi0 J/psi PHSP;
Enddecay
Decay J/psi
1.0 e+ e- PHSP;
Enddecay
Decay J/psi
1.0 mu+ mu- PHSP;
Enddecay
Decay pi0
1.0 gamma gamma PHSP;
Enddecay
End
DECAYCARD

# 27 energy points from 3.808 to 4.600 GeV (BOSS 703 scans)
scan_datasets = [
  DatasetManager.real_data.find("703_3810"),
  DatasetManager.real_data.find("703_3900"),
  DatasetManager.real_data.find("703_4009"),
  DatasetManager.real_data.find("703_4090"),
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4245"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("703_4310"),
  DatasetManager.real_data.find("703_4360"),
  DatasetManager.real_data.find("703_4390"),
  DatasetManager.real_data.find("703_4420"),
  DatasetManager.real_data.find("703_4470"),
  DatasetManager.real_data.find("703_4530"),
  DatasetManager.real_data.find("703_4575"),
  DatasetManager.real_data.find("703_4600"),
]

scan_incMC = scan_datasets.map { |ds| DatasetManager.inclusive_mc.find(ds) }

exMC = DatasetManager.create_exclusive_mc_for(scan_datasets) do |config|
  config.sample_name = "e_e_to_pi0_pi0_Jpsi"
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

algorithm = Algorithm.new("e_e_to_pi0_pi0_Jpsi")
algorithm.set_header(["e_e_to_pi0_pi0_JpsiAlg/e_e_to_pi0_pi0_Jpsi.h"])
# Multi-energy scan: no single ECMS constant

event_selection = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==1"
    nChrn "==1"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=4"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  }
  # Lepton PID is done via E/p cut (not probability-based):
  # electron: E/p > 0.7, muon: E/p < 0.3
  # Cannot be expressed via DSL pid block; implemented in C++ code
  .assign({ chrgp: :lp, chrgn: :lm })
  # Reconstruct two pi0 candidates via Kalman fits
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # 4C kinematic fit: 4-momentum conservation
  # (pi0 mass constraints were already applied in Kalman fits above)
  # Paper uses 6C fit — the extra 2 pi0 mass re-constraints in the main
  # kinematic fit are inexpressible in DSL (would require listing all 4
  # individual photons rather than the pi0 composites)
  .kinematic_fit([:lp, :lm, :pi0, :pi0]) {
    constrain_four_momentum
    chi2_cut 75
    nominal
  }

algorithm.note(:six_c_fit,
  "Paper uses 6C kinematic fit: 4-momentum conservation (4C) + two pi0 mass constraints (2C). " \
  "In DSL, pi0 composites from Kalman fits are already mass-constrained. " \
  "The extra pi0 mass re-constraints in the main kinematic fit are inexpressible."
)

algorithm.note(:lepton_ep_cut,
  "Lepton ID by E/p cut (not probability-based): electron E/p > 0.7, muon E/p < 0.3. " \
  "For J/psi -> mu+ mu- channel in PWA stage: at least one muon with >5 MUC layers."
)

algorithm.note(:selection_details,
  "pi0 mass window in Kalman fit: 0.11 < M(gamma gamma) < 0.15 GeV. " \
  "< 3 pi0 pi0 combinations per event. " \
  "Best combination chosen by minimum chi2_6C (< 75). " \
  "Uses ConExc generator for signal MC (Born cross section measurement with ISR correction). " \
  "27 energy points from 3.808 to 4.600 GeV. " \
  "J/psi region: 3.06 < M(l+l-) < 3.13 GeV/c^2. " \
  "Peaking backgrounds from eta J/psi and gamma psi(2S) subtracted."
)

algorithm.note(:born_cross_section,
  "Born cross section = N_obs / (L_int * (1+delta_r) * (1+delta_v) * epsilon * B_inter). " \
  "Radiative correction (1+delta_r) and vacuum polarization (1+delta_v) factors from QED. " \
  "Cross-section line shape fit for Y(4220)/Y(4320) resonant parameters. " \
  "ConExc generator used for signal MC simulation to handle ISR properly."
)

algorithm.with_decay_card(decay_card).apply(event_selection)
algorithm.execute_on(scan_datasets + scan_incMC + exMC)