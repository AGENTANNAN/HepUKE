# Paper 2107.09210v2: e+e- -> pi+ pi- psi(3686) cross section at sqrt(s)=4.0076-4.6984 GeV
# Mode I: psi(3686) -> pi+ pi- J/psi (charged decay)
# Mode II: psi(3686) -> neutrals + J/psi (neutral decay)

decay_card_modeI = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- psi(3686) PHSP;
  Enddecay
  Decay psi(3686)
  1.0000 pi+ pi- J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0000 e+ e- PHSP;
  Enddecay
  End
DECAYCARD

decay_card_modeII = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- psi(3686) PHSP;
  Enddecay
  Decay psi(3686)
  1.0000 pi+ pi- J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0000 e+ e- PHSP;
  Enddecay
  End
DECAYCARD

# Multi-energy scan from 4.0076 to 4.6984 GeV
data_scan = DatasetManager.real_data.where(cms_energy: {value: 4007..4700})
incMC_scan = DatasetManager.inclusive_mc.where(cms_energy: {value: 4007..4700})

sigMC_modeI = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "sig_pipi_psi3686_charged"
  config.events        = 100_000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

sigMC_modeII = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "sig_pipi_psi3686_neutral"
  config.events        = 100_000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

# === Mode I: psi(3686) -> pi+ pi- J/psi (charged decay, 5 or 6 charged tracks) ===
alg_modeI = Algorithm.new("PiPiPsi3686Charged")
alg_modeI.set_header(["PiPiPsi3686ChargedAlg/PiPiPsi3686Charged.h"])

sel_modeI = Selection.new
sel_modeI.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
end
.pid(method: :probability) do
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.7
  nep ">=1"
  nem ">=1"
end
.kinematic_fit([:pip, :pim, :pip, :pim, :ep, :em]) do
  nominal
  constrain_four_momentum
  chi2_cut 60
end

alg_modeI
  .note(:track_multiplicity, "5 charged tracks with net charge +-1 (allowing one missed pion), or 6 charged tracks with net charge 0")
  .note(:momentum_separation, "tracks with p > 1.0 GeV/c assigned as leptons; tracks with p < 0.65 (0.80) GeV/c for Ecms below (above) 4.465 GeV assigned as pions")
  .note(:lepton_separation, "electron/muon separation: E/p > 0.7 for electrons, E < 0.45 GeV for muons; lepton pair must be same flavor, opposite charge")
  .note(:kinematic_fits, "4C fit (6 tracks): e+e- -> pi+pi-pi+pi-l+l-; 1C fit (5 tracks): missing pion mass constraint, chi2_1C < 15; then J/psi mass window [3.05,3.15] GeV/c^2; 5C/2C fit after J/psi mass constraint")
  .note(:psi3686_signal, "psi(3686) signal from M(l+l-pi+pi-); 4 combinations, closest to nominal psi(3686) mass selected")
  .note(:signal_extraction, "simultaneous unbinned maximum likelihood fits to M(l+l-pi+pi-) and M_rec(pi+pi-) spectra; signal shape from MC convolved with Gaussian")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

datasets_modeI = data_scan + incMC_scan + sigMC_modeI
alg_modeI.execute_on(datasets_modeI)

# === Mode II: psi(3686) -> neutrals + J/psi (neutral decay, 4 charged tracks + photons) ===
alg_modeII = Algorithm.new("PiPiPsi3686Neutral")
alg_modeII.set_header(["PiPiPsi3686NeutralAlg/PiPiPsi3686Neutral.h"])

sel_modeII = Selection.new
sel_modeII.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
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
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.7
  nep ">=1"
  nem ">=1"
end

alg_modeII
  .note(:track_photon_multiplicity, "4 charged tracks with zero net charge; at least 2 good photon candidates; no kinematic fit performed in Mode II")
  .note(:momentum_separation, "tracks with p > 1.0 GeV/c assigned as leptons; tracks with p < 0.65 (0.80) GeV/c assigned as pions")
  .note(:lepton_separation, "electron/muon separation via E/p: E/p > 0.7 for e, E < 0.45 GeV for mu; same flavor, opposite charge lepton pair")
  .note(:background_veto, "cos(theta_pi+pi-) < 0.9 to suppress radiative Bhabha/dimuon; |M_rec(pi+pi-l+l-)| > 63 MeV/c^2 to suppress pi+pi-J/psi; |M_corr(psi(3686)) - M(psi(3686))| > 8 MeV/c^2 to suppress neutrals+psi(3686); |M(gamma gamma pi+pi-) - M(eta)| > 50 MeV/c^2 to suppress eta J/psi")
  .note(:jpsi_mass_window, "J/psi mass window [3.05, 3.15] GeV/c^2 on M(l+l-)")
  .note(:psi3686_signal, "psi(3686) signal from M_rec(pi+pi-) distribution")
  .note(:signal_extraction, "simultaneous unbinned maximum likelihood fits to M(l+l-pi+pi-) and M_rec(pi+pi-) spectra; shared cross section between two modes")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

datasets_modeII = data_scan + incMC_scan + sigMC_modeII
alg_modeII.execute_on(datasets_modeII)