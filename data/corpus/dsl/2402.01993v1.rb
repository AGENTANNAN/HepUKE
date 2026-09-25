# J/psi -> gamma eta', eta' -> pi+ pi- l+ l- (l=e,mu)
# BESIII: arXiv:2402.01993v1
# TFF measurement, CP-violation search, ALP search
# Two independent decay modes: e+e- mode and mu+mu- mode (Rule T1)

### Dataset ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for Mode I: eta' -> pi+ pi- e+ e-
decay_card_ee = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta' HELAMP 1 0 1 0 -1 0 -1 0 0 0;
  Enddecay
  Decay eta'
  1.000 pi+ pi- e+ e- VMD;
  Enddecay
  End
DECAYCARD

# Decay card for Mode II: eta' -> pi+ pi- mu+ mu-
decay_card_mumu = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta' HELAMP 1 0 1 0 -1 0 -1 0 0 0;
  Enddecay
  Decay eta'
  1.000 pi+ pi- mu+ mu- VMD;
  Enddecay
  End
DECAYCARD

exMC_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_etap_etap_pipi_ee"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_ee
  config.cross_section   = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_etap_etap_pipi_mumu"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_mumu
  config.cross_section   = :default
end

### Mode I Algorithm: eta' -> pi+ pi- e+ e- ###
alg_modeI = Algorithm.new("EtapPipiEE")
alg_modeI.set_header(["EtapPipiEEAlg/EtapPipiEE.h"])

sel_modeI = Selection.new
sel_modeI.select_track do
  cos_theta 0.93
  Vz        10.0
  Vr        1.0
  nChrp     "==2"
  nChrn     "==2"
  nNet      "==0"
end
.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    15.0
  nGam              ">=1"
end
.pid(method: :probability) do
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 1.0
  identify :pion, against: [:kaon, :proton]
  nlp  "==1"
  nlm  "==1"
  npip "==1"
  npim "==1"
end
.remove([:lp <= :chrgp, :lm <= :chrgn])
.remove([:pip <= :chrgp, :pim <= :chrgn])
# 4C kinematic fit: gamma pi+ pi- e+ e- (no PID in fit; PID-based best combination selection)
.kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

alg_modeI
  .note(:combined_chi2_cut, "chi2_4C+PID = chi2_4C + sum(chi2_PID_i) < 60; PID chi2 combined with 4C kinematic fit chi2 for best candidate selection; tight cut applied at ROOT level")
  .note(:photon_conversion_veto, "photon conversion veto: 2D cut on Rxy vs M_BP(e+e-) (vertex distance vs beam-pipe invariant mass) and Rxy vs Phi_ee (opening angle); rejects ~98.6% of conversion bg while keeping ~86% signal; veto curve parameters: (0.004 GeV/c2,0 cm), (0.004, 2 cm), (0.03, 3 cm), (0.07, 10 cm) in Rxy vs M_BP; reject Phi_ee < 75 deg when 2 < Rxy < 7.5 cm")
  .note(:etap_mass_window, "eta' mass window: 0.91 < M(pi+pi-e+e-) < 1.0 GeV/c2 for branching fraction; 0.945-0.97 for TFF/CPV/ALP")
  .note(:jpsi_events, "(10087+/-44) x 10^6 J/psi events used; J/psi -> gamma eta' radiative decay")
  .note(:vmd_generator, "signal MC generated with VMD model for eta'->pi+pi-l+l- decay")
  .with_decay_card(decay_card_ee)
  .apply(sel_modeI)

### Mode II Algorithm: eta' -> pi+ pi- mu+ mu- ###
alg_modeII = Algorithm.new("EtapPipiMuMu")
alg_modeII.set_header(["EtapPipiMuMuAlg/EtapPipiMuMu.h"])

sel_modeII = Selection.new
sel_modeII.select_track do
  cos_theta 0.93
  Vz        10.0
  Vr        1.0
  nChrp     "==2"
  nChrn     "==2"
  nNet      "==0"
end
.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    15.0
  nGam              ">=1"
end
.pid(method: :probability) do
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 1.0
  identify :pion, against: [:kaon, :proton]
  nlp  "==1"
  nlm  "==1"
  npip "==1"
  npim "==1"
end
.remove([:lp <= :chrgp, :lm <= :chrgn])
.remove([:pip <= :chrgp, :pim <= :chrgn])
# Nominal 4C fit: gamma pi+ pi- mu+ mu-
.kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end
# Competing-hypothesis veto: gamma pi+ pi- pi+ pi- (Rule T2)
# Assign muon slots to pion hypothesis, reuse track indices from nominal fit
.assign({lp: :pip, lm: :pim})
.kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) do
  use_track_index_from_nominal_kmfit
  constrain_four_momentum
end

alg_modeII
  .note(:chi2_4c_cut, "chi2_4C(pi+pi-mu+mu-) < 25 applied at BOSS level (tight cut from paper); chi2_4C+PID(gamma pi+pi-mu+mu-) < chi2_4C+PID(gamma 4pi) for background suppression")
  .note(:muon_mass_veto, "|M(mu+mu-) - 0.548| > 0.02 GeV/c2 to suppress eta'->pi+pi-eta, eta->mu+mu- background")
  .note(:etap_mass_window, "eta' mass window: 0.9 < M(pi+pi-mu+mu-) < 1.0 for BF; 0.945-0.975 for TFF/CPV/ALP")
  .note(:jpsi_events, "(10087+/-44) x 10^6 J/psi events used; J/psi -> gamma eta' radiative decay")
  .with_decay_card(decay_card_mumu)
  .apply(sel_modeII)

alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_ee])
alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_mumu])