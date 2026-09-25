DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# 18 energy points from 4.189 to 4.437 GeV (Table 1 of the paper)
scan_points = [
  DatasetManager.real_data.find("703_4190"),  # 4.189 GeV
  DatasetManager.real_data.find("703_4200"),  # 4.199 GeV
  DatasetManager.real_data.find("703_4210"),  # 4.209 GeV
  DatasetManager.real_data.find("703_4220"),  # 4.219 GeV
  DatasetManager.real_data.find("703_4230"),  # 4.226 GeV
  DatasetManager.real_data.find("703_4237"),  # 4.236 GeV
  DatasetManager.real_data.find("703_4246"),  # 4.244 GeV
  DatasetManager.real_data.find("703_4260"),  # 4.258 GeV
  DatasetManager.real_data.find("703_4270"),  # 4.267 GeV
  DatasetManager.real_data.find("703_4280"),  # 4.278 GeV
  DatasetManager.real_data.find("705_4290"),  # 4.288 GeV
  DatasetManager.real_data.find("705_4315"),  # 4.312 GeV
  DatasetManager.real_data.find("705_4340"),  # 4.338 GeV
  DatasetManager.real_data.find("703_4360"),  # 4.358 GeV
  DatasetManager.real_data.find("705_4380"),  # 4.378 GeV
  DatasetManager.real_data.find("705_4400"),  # 4.397 GeV
  DatasetManager.real_data.find("703_4420"),  # 4.416 GeV
  DatasetManager.real_data.find("705_4440"),  # 4.437 GeV
]

scan_incMCs = scan_points.map { |d| DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}") }

# Decay card: e+ e- -> pi+ pi- hc, hc -> pi0 J/psi (SIGNAL channel)
# Use KKMC + psi(4260) top mother; hc as intermediate
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 pi+ pi- hc PHSP;
    Enddecay
    Decay hc
    1.000 pi0 J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Decay card: e+ e- -> pi+ pi- hc, hc -> gamma eta_c, eta_c -> K+ K- pi0 (NORMALIZATION channel)
decay_card_norm = <<~DECAYCARD
    Decay psi(4260)
    1.000 pi+ pi- hc PHSP;
    Enddecay
    Decay hc
    1.000 gamma eta_c PHSP;
    Enddecay
    Decay eta_c
    1.000 K+ K- pi0 PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "hc_pi0_jpsi_signal"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

exMC_norm = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "hc_gamma_etac_norm"
  config.events        = 100_000
  config.decay_card    = decay_card_norm
  config.cross_section = :default
end

# Signal channel: hc -> pi0 J/psi (J/psi -> l+ l-, l = e or mu)
alg_signal = Algorithm.new("HcToPi0Jpsi")
alg_signal.set_header(["HcToPi0JpsiAlg/HcToPi0Jpsi.h"])
           .note(:lepton_id,
             "Tracks with p > 1.0 GeV/c assigned as leptons from J/psi decay. " \
             "Leptons with EMC energy > 1.0 GeV identified as electrons, < 0.4 GeV as muons. " \
             "Other tracks assigned as pions.")
           .note(:eta_omega_veto,
             "Events with M(pi+ pi- pi0) in [0.51, 0.58] GeV/c^2 (eta region) or " \
             "[0.75, 0.81] GeV/c^2 (omega region) are excluded to suppress background.")
           .note(:jpsi_mass_window,
             "J/psi signal region: M(l+ l-) in [3.085, 3.115] GeV/c^2. " \
             "Sidebands: [3.00, 3.06] and [3.14, 3.20] GeV/c^2.")
           .note(:hc_signal_extraction,
             "Simultaneous unbinned maximum-likelihood fit to RM(pi+ pi-) in signal and normalization " \
             "channels. Signal shape from MC convolved with Gaussian. Background: 1st-order polynomial.")
           .note(:multi_energy_weighting,
             "Branching ratio computed using luminosity, cross section, radiative correction factor " \
             "and efficiency at each energy point: sum_i L_i * sigma_i * (1+delta_i) * epsilon_i.")

sel_signal = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nChrp       "==2"
    nChrn       "==2"
    nNet        "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .kinematic_fit([:pip, :pim, :ep, :em, :pi0]) do
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 30
  end

alg_signal.with_decay_card(decay_card_signal).apply(sel_signal)

# Normalization channel: hc -> gamma eta_c, eta_c -> K+ K- pi0
alg_norm = Algorithm.new("HcToGammaEtac")
alg_norm.set_header(["HcToGammaEtacAlg/HcToGammaEtac.h"])
          .note(:pid,
            "TOF and dE/dx combined to form PID likelihoods for pi, K, p hypotheses. " \
            "Each track assigned the highest-likelihood hypothesis. " \
            "Exactly two oppositely charged pions and two oppositely charged kaons required.")
          .note(:etac_mass_window,
            "eta_c signal region: M(K+ K- pi0) in [2.92, 3.04] GeV/c^2. " \
            "Sidebands: [2.59, 2.71] and [3.25, 3.37] GeV/c^2.")

sel_norm = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nChrp       "==2"
    nChrn       "==2"
    nNet        "==0"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=3"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp  "==1"
    nkm  "==1"
    npip "==1"
    npim "==1"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .kinematic_fit([:pip, :pim, :kp, :km, :gamma, :pi0]) do
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 40
  end

alg_norm.with_decay_card(decay_card_norm).apply(sel_norm)

all_datasets = scan_points + scan_incMCs + exMC_signal + exMC_norm
alg_signal.execute_on(all_datasets)
alg_norm.execute_on(all_datasets)