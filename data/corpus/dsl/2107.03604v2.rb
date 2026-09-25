# Paper 2107.03604v2: e+e- -> gamma chi_c1,c2 cross sections at sqrt(s)=3.773-4.600 GeV
# Mode I: J/psi -> e+ e-

decay_card_modeI = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma chi_c1 PHSP;
  Enddecay
  Decay chi_c1
  1.0000 gamma J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0000 e+ e- PHSP;
  Enddecay
  End
DECAYCARD

decay_card_modeII = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma chi_c1 PHSP;
  Enddecay
  Decay chi_c1
  1.0000 gamma J/psi PHSP;
  Enddecay
  Decay J/psi
  1.0000 mu+ mu- PHSP;
  Enddecay
  End
DECAYCARD

# Data: multi-energy scan from 3.773 to 4.600 GeV
data_scan = DatasetManager.real_data.where(cms_energy: {value: 3770..4600})
incMC_scan = DatasetManager.inclusive_mc.where(cms_energy: {value: 3770..4600})

sigMC_modeI = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "sig_gam_chic1_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

sigMC_modeII = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name   = "sig_gam_chic1_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

# === Mode I: J/psi -> e+ e- ===
alg_modeI = Algorithm.new("GamChic1EE")
alg_modeI.set_header(["GamChic1EEAlg/GamChic1EE.h"])

sel_modeI = Selection.new
sel_modeI.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=1"
  nChrn ">=1"
  nNet "==0"
end
.for_each(:charged) do
  define(:momentum) { p }
  where { momentum < 1.0 }
  remove
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
                                 treat_as_electron_if_energy_above: 1.0
  nep ">=1"
  nem ">=1"
end
.kinematic_fit([:gamma, :gamma, :ep, :em]) do
  nominal
  constrain_four_momentum
  chi2_cut 40
end

alg_modeI
  .note(:background_veto, "radiative Bhabha veto: cos(theta_e_gamma) < 0.86; |cos(theta_gamma)| < 0.8 for both photons; |M(gamma_H,gamma_L) - m(eta)| > 0.03 GeV/c^2 for ee mode")
  .note(:jpsi_mass_window, "J/psi mass window [3.08, 3.12] GeV/c^2 on M(l+l-) after kinematic fit")
  .note(:signal_extraction, "chi_c1,c2 signal extracted from M(gamma J/psi) distribution; 2D fit at 4.009 GeV")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

datasets_modeI = data_scan + incMC_scan + sigMC_modeI
alg_modeI.execute_on(datasets_modeI)

# === Mode II: J/psi -> mu+ mu- ===
alg_modeII = Algorithm.new("GamChic1MuMu")
alg_modeII.set_header(["GamChic1MuMuAlg/GamChic1MuMu.h"])

sel_modeII = Selection.new
sel_modeII.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=1"
  nChrn ">=1"
  nNet "==0"
end
.for_each(:charged) do
  define(:momentum) { p }
  where { momentum < 1.0 }
  remove
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
                                 treat_as_electron_if_energy_above: 0.4
  nmup ">=1"
  nmum ">=1"
end
.kinematic_fit([:gamma, :gamma, :mup, :mum]) do
  nominal
  constrain_four_momentum
  chi2_cut 40
end

alg_modeII
  .note(:background_veto, "|M(gamma_H,gamma_L) - m(pi0)| > 0.015 GeV/c^2 for mumu mode to suppress pi+pi-pi0 background; |cos(theta_gamma)| < 0.8 for both photons")
  .note(:jpsi_mass_window, "J/psi mass window [3.08, 3.12] GeV/c^2 on M(l+l-) after kinematic fit")
  .note(:signal_extraction, "chi_c1,c2 signal extracted from M(gamma J/psi) distribution")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

datasets_modeII = data_scan + incMC_scan + sigMC_modeII
alg_modeII.execute_on(datasets_modeII)