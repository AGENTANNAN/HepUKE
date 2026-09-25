# Paper: 2310.03361v1 - e+e- -> eta J/psi at 44 energies from 3.808 to 4.951 GeV
# Mode I: eta -> gamma gamma, J/psi -> l+ l- (l = e/mu), 4C kinematic fit
# Mode II: eta -> pi0 pi+ pi-, J/psi -> l+ l- (l = e/mu), 5C kinematic fit

DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

# 44 energy scan points from 3.808 to 4.951 GeV
scan_energies = [
  "703_3810", "703_3872", "703_3900", "703_4009", "703_4090",
  "705_4130", "705_4160", "703_4180", "703_4190", "703_4200",
  "703_4210", "703_4220", "703_4230", "703_4237", "703_4246",
  "703_4260", "703_4270", "703_4280", "705_4290", "703_4310",
  "705_4315", "705_4340", "703_4360", "705_4380", "703_4390",
  "705_4400", "703_4420", "705_4440", "703_4470", "703_4530",
  "703_4575", "703_4600", "706_4610", "706_4620", "706_4640",
  "706_4660", "706_4680", "706_4700", "707_4740", "707_4750",
  "707_4780", "707_4840", "707_4946"
]
# Note: "707_4946" covers sqrt(s) ~ 4951 MeV

scan_data = scan_energies.map { |s| DatasetManager.real_data.find(s) }
scan_incMC = scan_energies.map { |s| DatasetManager.inclusive_mc.find(s) }

decay_card_eta_jpsi = <<~DECAYCARD
  Decay vpho
  1.0000 eta J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0000 e+ e- PHSP;
  Enddecay
  Decay J/psi
  1.0000 mu+ mu- PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  Decay eta
  1.0000 pi0 pi+ pi- PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

signal_mc_mode1 = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name     = "sig_eta_jpsi_mode1"
  config.events          = 100_000
  config.decay_card      = decay_card_eta_jpsi
  config.cross_section   = :default
end

signal_mc_mode2 = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name     = "sig_eta_jpsi_mode2"
  config.events          = 100_000
  config.decay_card      = decay_card_eta_jpsi
  config.cross_section   = :default
end

# ============================================================================
# Mode I: eta -> gamma gamma, J/psi -> l+ l-, 4C kinematic fit
# ============================================================================
algorithm_mode1 = Algorithm.new("EtaJpsi_ModeI", "00-00-01")

selection_mode1 = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nTot        "==2"    # two leptons
    nNet        "==0"    # zero net charge
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons do
      treat_as_lepton_if_momentum_above 1.0
      treat_as_electron_if_energy_above 0.8
      treat_as_muon_if_energy_below 0.4
    end
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    20.0
    nGam              ">=2"
  end
  .kinematic_fit([:gamma, :gamma, :lep_p, :lep_m]) do
    constrain_four_momentum
    chi2_cut 40
    nominal
  end
  .for_each(:gamma) do
    define(:e_gamma) { energy }
    where { e_gamma > 0.08 }
  end

algorithm_mode1
  .set_header(["EtaJpsiAlg/EtaJpsiAlg.h"])
  .note(:mode_I, "Mode I: eta -> gamma gamma, J/psi -> l+ l-. 4C kinematic fit under e+e- -> gamma gamma l+ l- hypothesis. chi2_4C < 40. Photon energy > 0.08 GeV after fit (suppresses radiative Bhabha/dimuon). Lepton identification: momentum > 1.0 GeV/c assigned as lepton; muons E <= 0.4 GeV, electrons E/pc > 0.8. J/psi mass window [3.067, 3.127] GeV/c2 applied in ROOT.")
  .note(:mode_II, "Mode II: eta -> pi0 pi+ pi-, J/psi -> l+ l-. 5C kinematic fit with pi0 mass constraint. chi2_5C < 80.")
  .note(:cross_section, "Cross section measured via simultaneous unbinned ML fit to M(gamma gamma) and M(pi0 pi+ pi-) spectra from J/psi -> e+e- and J/psi -> mu+mu- modes. ISR correction, vacuum polarization, and iterative weighting applied. Upper limits at 90% CL for low-statistics points.")
  .note(:scan_points, "44 energy points from 3.808 to 4.951 GeV. Dataset sample names: 703_3810 through 707_4946.")
  .note(:lepton_separation, "High-momentum lepton PID: tracks with p > 1.0 GeV/c assigned as leptons. Electron: E/pc > 0.8. Muon: E <= 0.4 GeV. Only same-flavor opposite-charge lepton pairs accepted.")
  .note(:Jpsi_mass_window, "J/psi mass window [3.067, 3.127] GeV/c2 applied on M(l+l-). Sideband regions [3.027, 3.057] and [3.137, 3.167] for non-J/psi background estimation (applied in ROOT).")
  .note(:photon_angle, "Photon-track opening angle > 20 degrees to suppress photons from charged track interactions.")
  .with_decay_card(decay_card_eta_jpsi)
  .apply(selection_mode1)

algorithm_mode1.execute_on(scan_data + scan_incMC + signal_mc_mode1)

# ============================================================================
# Mode II: eta -> pi0 pi+ pi-, J/psi -> l+ l-, 5C kinematic fit
# ============================================================================
algorithm_mode2 = Algorithm.new("EtaJpsi_ModeII", "00-00-01")

selection_mode2 = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nTot        "==4"     # two pions + two leptons
    nNet        "==0"     # zero net charge
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons do
      treat_as_lepton_if_momentum_above 1.0
      treat_as_electron_if_energy_above 0.8
      treat_as_muon_if_energy_below 0.4
    end
    identify :pion, against: [:kaon, :proton]
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    20.0
    nGam              ">=2"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .kinematic_fit([:pi0, :pip, :pim, :lep_p, :lep_m]) do
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 80
    nominal
  end

algorithm_mode2
  .set_header(["EtaJpsiAlg/EtaJpsiAlg.h"])
  .with_decay_card(decay_card_eta_jpsi)
  .apply(selection_mode2)

algorithm_mode2.execute_on(scan_data + scan_incMC + signal_mc_mode2)