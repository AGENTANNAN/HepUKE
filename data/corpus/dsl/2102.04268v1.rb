# Paper: 2102.04268v1
# Cross section measurement of e+e- -> p pbar eta and e+e- -> p pbar omega
# at 17 center-of-mass energies between 3.773 GeV and 4.5995 GeV
# BESIII multi-energy scan with two independent decay channels

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# 17 scan energy points
data_points = [
  DatasetManager.real_data.find("712_3773"),   # 3.7730 GeV
  DatasetManager.real_data.find("703_3810"),   # 3.8077 GeV
  DatasetManager.real_data.find("703_3872"),   # 3.8670 GeV
  DatasetManager.real_data.find("703_3900"),   # 3.8962 GeV
  DatasetManager.real_data.find("703_4009"),   # 4.0076 GeV
  DatasetManager.real_data.find("703_4090"),   # 4.0855 GeV
  DatasetManager.real_data.find("703_4180"),   # 4.1784 GeV
  DatasetManager.real_data.find("703_4190"),   # 4.1886 GeV
  DatasetManager.real_data.find("703_4200"),   # 4.1971 GeV
  DatasetManager.real_data.find("703_4210"),   # 4.2077 GeV
  DatasetManager.real_data.find("703_4220"),   # 4.2171 GeV
  DatasetManager.real_data.find("703_4230"),   # 4.2263 GeV
  DatasetManager.real_data.find("703_4237"),   # 4.2357 GeV
  DatasetManager.real_data.find("703_4260"),   # 4.2580 GeV
  DatasetManager.real_data.find("703_4360"),   # 4.3583 GeV
  DatasetManager.real_data.find("703_4420"),   # 4.4156 GeV
  DatasetManager.real_data.find("703_4600")    # 4.5995 GeV
]

incMC_points = data_points.map { |dp| DatasetManager.inclusive_mc.find(dp.sample_name) }

# ============================================================
# Channel A: e+e- -> p pbar eta, eta -> gamma gamma
# Signal MC: ConExc mode 51 (p pbar eta continuum production)
# Omit Particle vpho line for multi-energy auto-injection
# ============================================================
decay_card_chA = <<~DECAYCARD
  Decay vpho
  1 ConExc 51;
  Enddecay
  Decay vhdr
  1 p+ anti-p- eta PHSP;
  Enddecay
  Decay eta
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_chA = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_ppbar_eta"
  config.events        = 100_000
  config.decay_card    = decay_card_chA
  config.cross_section = :default
end

alg_chA = Algorithm.new("PPbarEta")
alg_chA.set_header(["PPbarEtaAlg/PPbarEta.h"])
        .note(:channel, "Channel A: e+e- -> p pbar eta, eta -> gamma gamma")
        .note(:multi_energy_scan, "Cross section measured at 17 c.m. energies from 3.773 to 4.5995 GeV. Per-point sqrt(s) from MeasuredEcmsSvc.")
        .note(:isr_correction, "ISR correction factor (1+delta) computed iteratively. Detection efficiency corrected iteratively.")

sel_chA = Selection.new
sel_chA.select_track {
          nChrp "==1"
          nChrn "==1"
          nNet  "==0"
          cos_theta 0.93
          Vz 10.0
          Vr 1.0
        }
        .select_photon {
          tdc_emc_start     0
          tdc_emc_end       14
          energyThreshold_b 0.025
          energyThreshold_e 0.050
          angle_to_track    10.0
          nGam              ">=2"
        }
        .pid(method: :probability) {
          prob_cut 0.001
          identify :proton, against: [:kaon, :pion]
          nprp "==1"
          nprm "==1"
        }
        .kalman_kinematic_fit([:gamma, :gamma]) {
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
          chi2_cut 25
          neta ">=1"
        }
        .kinematic_fit([:prp, :prm, :eta]) {
          nominal
          constrain_four_momentum
          chi2_cut 200
        }
        .note(:photon_barrel_endcap_cos_theta, "Barrel: |cos(theta)| < 0.80; Endcap: 0.86 < |cos(theta)| < 0.92")
        .note(:eta_mass_window, "eta mass window applied in ROOT: |M(gammagamma) - m_eta| < 3 sigma.")

alg_chA.with_decay_card(decay_card_chA).apply(sel_chA)
alg_chA.execute_on(data_points + incMC_points + exMC_chA)

# ============================================================
# Channel B: e+e- -> p pbar omega, omega -> pi+ pi- pi0, pi0 -> gamma gamma
# Signal MC: KKMC + psi(4260) (no predefined ConExc mode for p pbar omega)
# ============================================================
decay_card_chB = <<~DECAYCARD
  Decay psi(4260)
  1 p+ anti-p- omega PHSP;
  Enddecay
  Decay omega
  1 pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

exMC_chB = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_ppbar_omega"
  config.events        = 100_000
  config.decay_card    = decay_card_chB
  config.cross_section = :default
end

alg_chB = Algorithm.new("PPbarOmega")
alg_chB.set_header(["PPbarOmegaAlg/PPbarOmega.h"])
        .note(:channel, "Channel B: e+e- -> p pbar omega, omega -> pi+ pi- pi0, pi0 -> gamma gamma")
        .note(:multi_energy_scan, "Cross section measured at 17 c.m. energies from 3.773 to 4.5995 GeV. Per-point sqrt(s) from MeasuredEcmsSvc.")
        .note(:isr_correction, "ISR correction factor (1+delta) computed iteratively. Detection efficiency corrected iteratively.")
        .note(:background_veto, "Vetoes applied in ROOT: K_S0 -> pi+pi-, Lambda -> p pi- peaking backgrounds. Mass window cuts on intermediate resonances.")

sel_chB = Selection.new
sel_chB.select_track {
          nChrp "==2"
          nChrn "==2"
          nNet  "==0"
          cos_theta 0.93
          Vz 10.0
          Vr 1.0
        }
        .select_photon {
          tdc_emc_start     0
          tdc_emc_end       14
          energyThreshold_b 0.025
          energyThreshold_e 0.050
          angle_to_track    10.0
          nGam              ">=2"
        }
        .pid(method: :probability) {
          prob_cut 0.001
          identify :proton, against: [:kaon, :pion]
          identify :pion, against: [:kaon, :proton]
          nprp "==1"
          nprm "==1"
          npip "==1"
          npim "==1"
        }
        .kalman_kinematic_fit([:gamma, :gamma]) {
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          chi2_cut 25
          npi0 ">=1"
        }
        .kinematic_fit([:prp, :prm, :pip, :pim, :pi0]) {
          nominal
          constrain_four_momentum
          chi2_cut 200
        }
        .note(:photon_barrel_endcap_cos_theta, "Barrel: |cos(theta)| < 0.80; Endcap: 0.86 < |cos(theta)| < 0.92")
        .note(:omega_mass_window, "omega mass window applied in ROOT: |M(pi+pi-pi0) - m_omega| < 30 MeV/c^2.")
        .note(:signal_extraction, "Signal yield extracted from simultaneous fit to M(ppbar) recoil mass spectra. Born cross section from N_sig / (L_int * epsilon * (1+delta) * B).")

alg_chB.with_decay_card(decay_card_chB).apply(sel_chB)
alg_chB.execute_on(data_points + incMC_points + exMC_chB)