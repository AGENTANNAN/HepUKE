# Search for chi_cJ -> mu+ mu- J/psi at psi(3686)
# Paper: 1901.06627v2, using 4.48e8 psi(3686) events

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Two J/psi decay modes require separate algorithms per Rule T1
# Mode A: J/psi -> e+ e-
# Mode B: J/psi -> mu+ mu-

# Mode A: J/psi -> e+ e-
decay_card_A = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_cJ PHSP;
  Enddecay

  Decay chi_cJ
  1.000 mu+ mu- J/psi PHSP;
  Enddecay

  Decay J/psi
  1.000 e+ e- PHSP;
  Enddecay

  End
DECAYCARD

exMC_A = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "chicJ_mumu_Jpsi_ee_exMC"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_A
  config.cross_section = :default
end

alg_A = Algorithm.new("ChicJMumuJpsiEE")
alg_A.set_header(["ChicJMumuJpsiEEAlg/ChicJMumuJpsiEE.h"])
     .set_constant({"ECMS" => [:double, 3.686]})

sel_A = Selection.new
sel_A.select_track do
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
  angle_to_track 20.0
  nGam ">=1"
end
.pid(method: :probability) do
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 1.0
  identify :muon, against: [:pion, :kaon, :proton]
  nlp "==1"
  nlm "==1"
  nmup ">=1"
  nmum ">=1"
end
.kinematic_fit([:gamma, :mup, :mum, :ep, :em]) do
  nominal
  constrain_four_momentum
  chi2_cut 40
end

alg_A
  .note(:momentum_sorting, "tracks with p>1 GeV/c assigned as J/psi daughters (e/mu by E/p: E>1GeV->e, E<0.3GeV->mu); low-p tracks assigned as chi_cJ muons")
  .note(:jpsi_mass_window, "J/psi lepton pair invariant mass required in [3.085, 3.110] GeV/c^2")
  .note(:vertex_fit_prefit, "vertex fit on all 4 charged tracks performed before 4C kinematic fit")
  .note(:photon_conversion_veto, "R_xy > 8.5 cm to veto photon-conversion background; R_xy is distance from reconstructed conversion vertex to z-axis")
  .note(:eta_jpsi_veto, "M(gamma mu+ mu-) < 0.535 or > 0.560 GeV/c^2 to suppress psi(3686)->eta J/psi background")
  .note(:helix_correction, "helix parameter correction applied to charged tracks before 4C kinematic fit; determined from psi(3686)->pi+pi-J/psi control sample")
  .note(:photon_best, "if multiple photon candidates, the one with smallest 4C fit chi^2 is retained")
  .with_decay_card(decay_card_A)
  .apply(sel_A)

alg_A.execute_on([psip_data, psip_incMC, exMC_A])


# Mode B: J/psi -> mu+ mu-
decay_card_B = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_cJ PHSP;
  Enddecay

  Decay chi_cJ
  1.000 mu+ mu- J/psi PHSP;
  Enddecay

  Decay J/psi
  1.000 mu+ mu- PHSP;
  Enddecay

  End
DECAYCARD

exMC_B = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "chicJ_mumu_Jpsi_mumu_exMC"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_B
  config.cross_section = :default
end

alg_B = Algorithm.new("ChicJMumuJpsiMuMu")
alg_B.set_header(["ChicJMumuJpsiMuMuAlg/ChicJMumuJpsiMuMu.h"])
     .set_constant({"ECMS" => [:double, 3.686]})

sel_B = Selection.new
sel_B.select_track do
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
  angle_to_track 20.0
  nGam ">=1"
end
.pid(method: :probability) do
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.3
  identify :muon, against: [:pion, :kaon, :proton]
  nmup ">=2"
  nmum ">=2"
end
.kinematic_fit([:gamma, :mup, :mum, :mup, :mum]) do
  nominal
  constrain_four_momentum
  chi2_cut 40
end

alg_B
  .note(:momentum_sorting, "tracks with p>1 GeV/c assigned as J/psi daughters (identified as muons by E/p<0.3 GeV); low-p tracks assigned as chi_cJ muons")
  .note(:jpsi_mass_window, "J/psi muon pair invariant mass required in [3.085, 3.110] GeV/c^2")
  .note(:vertex_fit_prefit, "vertex fit on all 4 charged tracks performed before 4C kinematic fit")
  .note(:photon_conversion_veto, "R_xy > 8.5 cm to veto photon-conversion background; R_xy is distance from reconstructed conversion vertex to z-axis")
  .note(:eta_jpsi_veto, "M(gamma mu+ mu-) < 0.535 or > 0.560 GeV/c^2 to suppress psi(3686)->eta J/psi background")
  .note(:helix_correction, "helix parameter correction applied to charged tracks before 4C kinematic fit; determined from psi(3686)->pi+pi-J/psi control sample")
  .note(:photon_best, "if multiple photon candidates, the one with smallest 4C fit chi^2 is retained")
  .with_decay_card(decay_card_B)
  .apply(sel_B)

alg_B.execute_on([psip_data, psip_incMC, exMC_B])