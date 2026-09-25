# e+e- -> eta_c pi+ pi- pi0 cross section at 6 CMS energies (4.178-4.600 GeV)
# eta_c reconstructed via 16 hadronic decay channels simultaneously
# Also: e+e- -> eta_c pi+ pi- and e+e- -> eta_c pi0 gamma (upper limits)
# Z_c -> eta_c pi search at 4.23 GeV

# Energy scan data points for cross-section measurement
scan_data_etac = [
  DatasetManager.real_data.find("703_4180"),   # 4.178 GeV, 3194.5 pb-1
  DatasetManager.real_data.find("703_4230"),   # 4.226 GeV (4230), 1091.7 pb-1 (actually has two entries)
  DatasetManager.real_data.find("703_4260"),   # 4.258 GeV (4258-4260), 825.7 pb-1
  DatasetManager.real_data.find("703_4360"),   # 4.358 GeV, 539.8 pb-1
  DatasetManager.real_data.find("703_4420"),   # 4.416 GeV, 1073.6 pb-1
  DatasetManager.real_data.find("703_4600"),   # 4.600 GeV, 566.9 pb-1
]

# Decay card: e+e- -> eta_c pi+ pi- pi0, with eta_c -> K_S0 K+ pi- (mode 05, 2.43%)
decay_card_etac_pipipi0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta_c pi+ pi- pi0  PHSP;
  Enddecay

  Decay eta_c
  1.0000 K_S0 K+ pi-  PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-  PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay

  End
DECAYCARD

exMC_etac_pipipi0 = DatasetManager.create_exclusive_mc_for(scan_data_etac) do |config|
  config.sample_name = "sig_etac_pipipi0"
  config.events = 100_000
  config.decay_card = decay_card_etac_pipipi0
  config.cross_section = :default
end

algorithm_etac_pipipi0 = Algorithm.new("EtacPipipi0")
algorithm_etac_pipipi0.set_header(["EtacPipipi0Alg/EtacPipipi0.h"])
                       .note(:simultaneous_fit_16_channels,
                         "eta_c reconstructed via 16 hadronic decay channels simultaneously: " \
                         "3(pi+pi-), 2(pi+pi-pi0), pi+pi-pi0pi0, 2(pi+pi-), " \
                         "K_S0 K+ pi-, K+K-pi+pi-, K+K-pi0, K_S0 K+ pi- pi+ pi-, " \
                         "2(pi+pi-)eta, pi+pi-eta, K+K-eta, K+K-K+K-, " \
                         "K+K-2(pi+pi-), p pbar, p pbar pi+ pi-, p pbar pi0. " \
                         "Common cross-section fit performed in ROOT; " \
                         "only representative decay channel shown here.")
                       .note(:isr_correction,
                         "ISR radiative correction via iterative procedure using " \
                         "Y(4260) lineshape as input. Radiative correction factors " \
                         "kappa_i depend on eta_c decay channel and sqrt(s).")
                       .note(:Dmeson_veto,
                         "Events rejected if D-meson candidate reconstructed in: " \
                         "D0->K-pi+, D0->K-pi+pi0, D+->K-pi+pi+, " \
                         "D+->K_S0 pi+, D+->K_S0 pi+ pi0. " \
                         "Also veto K*(892)->Kpi, omega->pi+pi-pi0, eta->pi+pi-pi0. " \
                         "Vetoes optimized per channel via FOM = S/sqrt(B).")

event_selection_etac_pipipi0 = Selection.new
event_selection_etac_pipipi0.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=1"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=2"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  identify :pion, against: [:kaon, :proton]
  nkp "==1"
  npim ">=1"
}
.remove([:kp <= :chrgp])
# Reconstruct K_S0 -> pi+ pi-
.secondary_vertex_fit([:pip, :pim]) {
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
}
# Reconstruct pi0 from photon pairs
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0 ">=1"
}
# Main kinematic fit: eta_c candidates from K_S0 K+ pi-
# with recoil pi+ pi- pi0
.kinematic_fit([:K_S0, :kp, :pim, :pip, :pim, :pi0]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}
.note(:eta_c_mass_window,
  "Z_c search applies additional cut: 2.880 < M_eta_c < 3.080 GeV/c^2.")
.note(:Zc_search,
  "Z_c search at sqrt(s)=4.23 GeV only: " \
  "scan over 10 masses (3625-3805 MeV/c^2) x 4 widths (8-38 MeV). " \
  "Charged Z_c^+- -> eta_c pi^+- and neutral Z_c^0 -> eta_c pi0. " \
  "No significant signal found; max significance 3.2 sigma " \
  "(neutral, m=3685 MeV/c^2, Gamma=28 MeV), reduced to ~2 sigma " \
  "after Look-elsewhere effect correction.")

algorithm_etac_pipipi0.with_decay_card(decay_card_etac_pipipi0)
                       .apply(event_selection_etac_pipipi0)
algorithm_etac_pipipi0.execute_on(scan_data_etac + exMC_etac_pipipi0)


# e+e- -> eta_c pi+ pi- (upper limits only)
decay_card_etac_pipi = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta_c pi+ pi-  PHSP;
  Enddecay

  Decay eta_c
  1.0000 K_S0 K+ pi-  PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-  PHSP;
  Enddecay

  End
DECAYCARD

exMC_etac_pipi = DatasetManager.create_exclusive_mc_for(scan_data_etac) do |config|
  config.sample_name = "sig_etac_pipi"
  config.events = 100_000
  config.decay_card = decay_card_etac_pipi
  config.cross_section = :default
end

algorithm_etac_pipi = Algorithm.new("EtacPipi")
algorithm_etac_pipi.set_header(["EtacPipiAlg/EtacPipi.h"])
                    .set_constant({"ECMS" => [:double, 4.260]})
                    .note(:upper_limits_only,
                      "No significant eta_c production observed in e+e- -> eta_c pi+ pi-. " \
                      "90% CL upper limits: 3-19 pb across sqrt(s). " \
                      "Conservative ISR correction using narrow resonance assumption " \
                      "(Gamma=10 MeV) for kappa_min values.")

event_selection_etac_pipi = Selection.new
event_selection_etac_pipi.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=1"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  identify :pion, against: [:kaon, :proton]
  nkp "==1"
  npim ">=1"
}
.remove([:kp <= :chrgp])
.secondary_vertex_fit([:pip, :pim]) {
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
}
.kinematic_fit([:K_S0, :kp, :pim, :pip, :pim]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

algorithm_etac_pipi.with_decay_card(decay_card_etac_pipi)
                    .apply(event_selection_etac_pipi)
algorithm_etac_pipi.execute_on(scan_data_etac + exMC_etac_pipi)