# DSL for 2003.03705v1: e+e- -> eta J/psi cross section measurement
# Energies: 3.81-4.60 GeV, two signal modes
# Mode I: J/psi -> l+l-, eta -> gamma gamma
# Mode II: J/psi -> l+l-, eta -> pi+ pi- pi0

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# Energy scan points
scan_points = [
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
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
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

incMC_points = scan_points.map do |dp|
  DatasetManager.inclusive_mc.find("#{dp.boss}_#{dp.sample_name}")
end

# Decay card: ConExc mode 80 = eta J/psi
decay_card_eta_jpsi = <<~DECAYCARD
  Decay vpho
  1 ConExc 80;
  Enddecay
  Decay J/psi
  1 e+ e- PHSP;
  Enddecay
  End
DECAYCARD

exMC_scan = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_eta_jpsi"
  config.events        = 100_000
  config.decay_card    = decay_card_eta_jpsi
  config.cross_section = :default
end

# --- Mode I: J/psi -> e+e-, eta -> gamma gamma ---

alg_modeI = Algorithm.new("EtaJpsiModeI")
alg_modeI.set_header(["EtaJpsiModeIAlg/EtaJpsiModeI.h"])

decay_card_modeI = <<~DECAYCARD
  Decay vpho
  1 ConExc 80;
  Enddecay
  Decay J/psi
  1 e+ e- VLL;
  Enddecay
  Decay eta
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

sel_modeI = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==1"
    nChrn "==1"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=2"
  end
  .pid(method: :probability) do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  end
  .kinematic_fit([:ep, :em, :gamma, :gamma]) do
    nominal
    constrain_four_momentum
    invariant_mass_of(:ep, :em).constrain_to_nominal_mass_of(:jpsi)
    chi2_cut 200
  end

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on(scan_points + incMC_points + exMC_scan)

# --- Mode II: J/psi -> e+e-, eta -> pi+ pi- pi0 ---

alg_modeII = Algorithm.new("EtaJpsiModeII")
alg_modeII.set_header(["EtaJpsiModeIIAlg/EtaJpsiModeII.h"])

decay_card_modeII = <<~DECAYCARD
  Decay vpho
  1 ConExc 80;
  Enddecay
  Decay J/psi
  1 e+ e- VLL;
  Enddecay
  Decay eta
  1 pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

sel_modeII = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==2"
    nChrn "==2"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=2"
  end
  .pid(method: :probability) do
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  end
  .remove([:lp <= :chrgp, :lm <= :chrgn])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=1"
    npim ">=1"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .kinematic_fit([:ep, :em, :pip, :pim, :pi0]) do
    nominal
    constrain_four_momentum
    invariant_mass_of(:ep, :em).constrain_to_nominal_mass_of(:jpsi)
    chi2_cut 200
  end

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on(scan_points + incMC_points + exMC_scan)