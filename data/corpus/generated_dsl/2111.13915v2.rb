### Dataset description ###
# 18 BESIII energy points of the e+e- -> pi+pi- h_c scan, 4.189 - 4.437 GeV
data_points = [
  DatasetManager.real_data.find("703_4190"),  # 4.1888 GeV
  DatasetManager.real_data.find("703_4200"),  # 4.1989 GeV
  DatasetManager.real_data.find("703_4210"),  # 4.2092 GeV
  DatasetManager.real_data.find("703_4220"),  # 4.2187 GeV
  DatasetManager.real_data.find("703_4230"),  # 4.2263 GeV
  DatasetManager.real_data.find("703_4237"),  # 4.2357 GeV
  DatasetManager.real_data.find("703_4245"),  # 4.2417 GeV
  DatasetManager.real_data.find("703_4246"),  # 4.2438 GeV
  DatasetManager.real_data.find("703_4260"),  # 4.2580 GeV
  DatasetManager.real_data.find("703_4270"),  # 4.2668 GeV
  DatasetManager.real_data.find("703_4280"),  # 4.2777 GeV
  DatasetManager.real_data.find("705_4290"),  # 4.2879 GeV
  DatasetManager.real_data.find("703_4310"),  # 4.3079 GeV
  DatasetManager.real_data.find("705_4315"),  # 4.3121 GeV
  DatasetManager.real_data.find("705_4340"),  # 4.3374 GeV
  DatasetManager.real_data.find("705_4380"),  # 4.3774 GeV
  DatasetManager.real_data.find("705_4400"),  # 4.3965 GeV
  DatasetManager.real_data.find("705_4440")   # 4.4362 GeV
]

# Corresponding inclusive MC samples, one per energy point
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230"),
  DatasetManager.inclusive_mc.find("703_4237"),
  DatasetManager.inclusive_mc.find("703_4245"),
  DatasetManager.inclusive_mc.find("703_4246"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4270"),
  DatasetManager.inclusive_mc.find("703_4280"),
  DatasetManager.inclusive_mc.find("705_4290"),
  DatasetManager.inclusive_mc.find("703_4310"),
  DatasetManager.inclusive_mc.find("705_4315"),
  DatasetManager.inclusive_mc.find("705_4340"),
  DatasetManager.inclusive_mc.find("705_4380"),
  DatasetManager.inclusive_mc.find("705_4400"),
  DatasetManager.inclusive_mc.find("705_4440")
]

# Decay card for the signal: e+e- -> pi+pi- h_c, h_c -> pi0 J/psi, J/psi -> e+e-
decay_card_signal_ee = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0000 pi0 J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for the signal: e+e- -> pi+pi- h_c, h_c -> pi0 J/psi, J/psi -> mu+mu-
decay_card_signal_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0000 pi0 J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 mu+ mu- PHOTOS VLL;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for the normalization channel: h_c -> gamma eta_c, eta_c -> K+ K- pi0
decay_card_normalization = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- h_c PHSP;
  Enddecay

  Decay h_c
  1.0000 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.0000 K+ K- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 100k-event exclusive MC for every decay mode at every energy point of the scan
exMCs_signal_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_pipihc_pi0jpsi_ee"
  config.events = 100_000
  config.decay_card = decay_card_signal_ee
  config.cross_section = :default
end

exMCs_signal_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_pipihc_pi0jpsi_mumu"
  config.events = 100_000
  config.decay_card = decay_card_signal_mumu
  config.cross_section = :default
end

exMCs_normalization = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "norm_pipihc_gammaetac"
  config.events = 100_000
  config.decay_card = decay_card_normalization
  config.cross_section = :default
end

### Event selection (BOSS) ###
# ------------------------------------------------------------------
# Signal channel: e+e- -> pi+pi- h_c, h_c -> pi0 J/psi, J/psi -> l+l-
# ------------------------------------------------------------------
alg_name_signal = "PiPiHcPi0Jpsi"
alg_signal = Algorithm.new(alg_name_signal)
alg_signal.set_header(["#{alg_name_signal}Alg/#{alg_name_signal}.h"])
          .set_constant({"ECMS" => [:double, 4.258]}) # nominal beam energy (each dataset runs at its own point)
          .set_alias({"std::vector<double>" => "Vdouble"})

# The extra lepton-definition detail (muon < 0.4 GeV EMC in the 0.4-1.0 GeV gap)
# cannot be encoded by identify_high_momentum_leptons' single energy threshold.
alg_signal.note(:lepton_pid_definition,
  "J/psi leptons are taken from charged tracks with p > 1.0 GeV/c; a lepton is an " \
  "electron if its EMC energy is > 1.0 GeV and a muon if < 0.4 GeV, all remaining " \
  "tracks being treated as pions. The DSL lepton helper supports only one " \
  "electron/muon energy threshold (set to 1.0 GeV), so the additional muon < 0.4 GeV " \
  "requirement inside the 0.4-1.0 GeV band is recorded here.")

sel_signal = Selection.new
  .select_track {          # charged track selection
    cos_theta 0.93         # |cos(theta)| < 0.93
    Vz 10.0                # |Vz| < 10 cm
    Vr 1.0                 # Vr < 1 cm
    nChrp  "==2"           # exactly two positive tracks
    nChrn  "==2"           # exactly two negative tracks
    nNet   "==0"           # net charge zero
  }
  .select_photon {         # photon selection
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025 # 25 MeV barrel
    energyThreshold_e 0.050 # 50 MeV endcap
    nGam ">=2"              # at least two photons (pi0 -> gamma gamma)
  }
  .pid(method: :probability) {
    # high-momentum tracks -> J/psi leptons (e or mu); rest -> pions
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 1.0
    identify :pion, against: [:kaon, :proton]
    nlp  "==1"              # one positive lepton
    nlm  "==1"              # one negative lepton
    npip "==1"              # one pi+
    npim "==1"              # one pi-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {  # reconstruct pi0 (1-C mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"              # at least one pi0 candidate
  }
  .kinematic_fit([:pip, :pim, :lp, :lm, :pi0]) { # 4C fit + pi0 mass constraint
    nominal
    constrain_four_momentum
    invariant_mass_of(:lp, :lm).within(3.085, 3.115)        # J/psi mass window
    invariant_mass_of(:pip, :pim, :pi0).out_of(0.51, 0.58)  # veto eta region
    invariant_mass_of(:pip, :pim, :pi0).out_of(0.75, 0.81)  # veto omega region
    chi2_cut 30
  }

# ------------------------------------------------------------------
# Normalization channel: h_c -> gamma eta_c, eta_c -> K+ K- pi0
# ------------------------------------------------------------------
alg_name_norm = "PiPiHcGammaEtaC"
alg_normalization = Algorithm.new(alg_name_norm)
alg_normalization.set_header(["#{alg_name_norm}Alg/#{alg_name_norm}.h"])
                 .set_constant({"ECMS" => [:double, 4.258]})
                 .set_alias({"std::vector<double>" => "Vdouble"})

sel_normalization = Selection.new
  .select_track {          # same charged track selection as the signal
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp  "==2"           # pi+ K+ (two positive)
    nChrn  "==2"           # pi- K- (two negative)
    nNet   "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=3"             # at least three photons (gamma + pi0 -> gamma gamma)
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]   # K+ and K-
    identify :pion, against: [:kaon, :proton]   # pi+ and pi-
    nkp  "==1"
    nkm  "==1"
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 reconstruction (1-C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pip, :pim, :kp, :km, :gamma, :pi0]) { # 4C fit + pi0 constraint
    nominal
    constrain_four_momentum
    invariant_mass_of(:kp, :km, :pi0).within(2.92, 3.04)   # eta_c signal region
    chi2_cut 40
  }

# Generate the BOSS algorithms (one decay card per algorithm)
alg_signal.with_decay_card(decay_card_signal_ee).apply(sel_signal)
alg_normalization.with_decay_card(decay_card_normalization).apply(sel_normalization)

# Execute on real data, inclusive MC and the whole scan's exclusive MC samples
root_files_signal = alg_signal.execute_on(
  data_points + incMC_points + exMCs_signal_ee + exMCs_signal_mumu
)
root_files_normalization = alg_normalization.execute_on(
  data_points + incMC_points + exMCs_normalization
)