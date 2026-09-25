# =============================================================================
#  e+e- -> eta eta J/psi  (J/psi -> e+e- / mu+mu- , eta -> gamma gamma / pi+pi-pi0)
#  XYZ scan ~4.226-4.950 GeV : eight energy points anchored on the 4260 point
# =============================================================================

### ------------------------------- Datasets ------------------------------- ###
# Eight scan points (4.260 -> 4.950 GeV); 4260 is the anchor point.
energy_points = [
  DatasetManager.real_data.find("703_4260"),  # 4.260 GeV (anchor)
  DatasetManager.real_data.find("703_4360"),  # 4.360 GeV
  DatasetManager.real_data.find("703_4420"),  # 4.420 GeV
  DatasetManager.real_data.find("703_4600"),  # 4.600 GeV
  DatasetManager.real_data.find("707_4740"),  # 4.740 GeV
  DatasetManager.real_data.find("707_4780"),  # 4.780 GeV
  DatasetManager.real_data.find("707_4840"),  # 4.840 GeV
  DatasetManager.real_data.find("707_4946")   # 4.950 GeV
]

# Matching inclusive MC for each scan point
incmc_points = [
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4360"),
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("707_4740"),
  DatasetManager.inclusive_mc.find("707_4780"),
  DatasetManager.inclusive_mc.find("707_4840"),
  DatasetManager.inclusive_mc.find("707_4946")
]

### ------------------------------ Decay cards ----------------------------- ###
# Four-photon mode, J/psi -> e+e-   (psi(4260) is the KKMC top mother)
decay_4gamma_ee = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta eta J/psi PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Four-photon mode, J/psi -> mu+mu-
decay_4gamma_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta eta J/psi PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  1.0000 mu+ mu- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# eta -> pi+pi-pi0 mode, J/psi -> e+e-
decay_pipipi0_ee = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta eta J/psi PHSP;
  Enddecay

  Decay eta
  1.0000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# eta -> pi+pi-pi0 mode, J/psi -> mu+mu-
decay_pipipi0_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta eta J/psi PHSP;
  Enddecay

  Decay eta
  1.0000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  1.0000 mu+ mu- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

### --------------- Exclusive MC : 100k events per energy point ------------ ###
exMC_4gamma_ee = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_etaeajpsi_4gamma_ee"
  config.events        = 100_000
  config.decay_card    = decay_4gamma_ee
  config.cross_section = :default
end

exMC_4gamma_mumu = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_etaeajpsi_4gamma_mumu"
  config.events        = 100_000
  config.decay_card    = decay_4gamma_mumu
  config.cross_section = :default
end

exMC_pipipi0_ee = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_etaeajpsi_pipipi0_ee"
  config.events        = 100_000
  config.decay_card    = decay_pipipi0_ee
  config.cross_section = :default
end

exMC_pipipi0_mumu = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_etaeajpsi_pipipi0_mumu"
  config.events        = 100_000
  config.decay_card    = decay_pipipi0_mumu
  config.cross_section = :default
end

### ------------------- Common part of the event selection ----------------- ###
# Charged tracks / photons / PID are identical for all eta eta J/psi topologies.
common_selection = Selection.new
common_selection
  .select_track {                       # charged track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz 10.0                             # |Vz| < 10 cm
    Vr 1.0                              # |Vr| < 1 cm
    nChrp ">=1"                         # at least one positive track
    nChrn ">=1"                         # at least one negative track
    nNet  "==0"                         # net charge zero
  }
  .select_photon {                      # photon selection
    tdc_emc_start 0                     # TDC start 0 (700 ns)
    tdc_emc_end 14                      # TDC end 14 (700 ns)
    angle_to_track 10.0                 # angle to nearest charged track > 10 deg
    energyThreshold_b 0.025             # 25 MeV barrel
    energyThreshold_e 0.050             # 50 MeV endcap
    nGam ">=4"                          # at least four photons
  }
  .pid(method: :probability) {          # probability PID
    prob_cut 0.001                      # probability > 0.001
    # tracks with p > 1.0 GeV treated as leptons; electron if EMC E > 0.6 GeV else muon
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon, :proton]   # pions (p < 1.0 GeV) vs K/p
  }

### ======================================================================== ###
###  Mode I : eta -> gamma gamma  (four-photon mode),  J/psi -> l+l-
### ======================================================================== ###
alg_name_4g = "EtaEtaJpsi4Gamma"
alg_4g = Algorithm.new(alg_name_4g)
alg_4g.set_header(["#{alg_name_4g}Alg/#{alg_name_4g}.h"])
      .set_constant({"ECMS" => [:double, 4.260]})   # 4260 anchor (per-run energy used at execution)
      .set_alias({"std::vector<double>" => "Vdouble"})
      .note(:efficiency_curve,
            "ECMS varies across the 4.226-4.950 GeV scan; the declared constant is the 4260 anchor and the per-run measured CMS energy is used at execution")
      .note(:background_veto,
            "four-photon mode: veto events favouring two pi0 over two eta by comparing the summed chi2 of the two 1C mass-constrained hypotheses (gamma-gamma -> eta versus gamma-gamma -> pi0); applied in ROOT from the stored competing chi2")

sel_4g = common_selection.dup
sel_4g
  # pi0 -> gamma-gamma candidates via 1C Kalman fit (window 80-180 MeV)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).between(0.08, 0.18)         # pre-selection window
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=0"
  }
  # eta -> gamma-gamma candidates via 1C Kalman fit (window 400-700 MeV)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).between(0.40, 0.70)         # pre-selection window
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  # Nominal 7C fit : 4C + M(gamma gamma)=M(eta) x2 + M(l+l-)=M(J/psi)
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :lp, :lm]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    chi2_cut 200                                                  # loose cut; tight cut in ROOT
  }
  # Competing pi0 pi0 J/psi hypothesis : stores chi2 only (no cut / not nominal)
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :lp, :lm]) {
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
  }

alg_4g.with_decay_card(decay_4gamma_ee).apply(sel_4g)

### ======================================================================== ###
###  Mode II : eta -> pi+pi-pi0 ,  J/psi -> l+l-
### ======================================================================== ###
alg_name_pp0 = "EtaEtaJpsiPiPiPi0"
alg_pp0 = Algorithm.new(alg_name_pp0)
alg_pp0.set_header(["#{alg_name_pp0}Alg/#{alg_name_pp0}.h"])
       .set_constant({"ECMS" => [:double, 4.260]})
       .set_alias({"std::vector<double>" => "Vdouble"})
       .note(:efficiency_curve,
             "ECMS varies across the 4.226-4.950 GeV scan; the declared constant is the 4260 anchor and the per-run measured CMS energy is used at execution")

sel_pp0 = common_selection.dup
sel_pp0
  # pi0 -> gamma-gamma candidates via 1C Kalman fit (window 80-180 MeV); two pi0 expected
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).between(0.08, 0.18)         # pre-selection window
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  # Nominal 7C fit : 4C + M(pi+pi-pi0)=M(eta) x2 + M(l+l-)=M(J/psi)
  .kinematic_fit([:pip, :pim, :pi0, :pip, :pim, :pi0, :lp, :lm]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:pip, :pim, :pi0).between(0.40, 0.70)       # eta -> pi+pi-pi0 pre-selection window
    invariant_mass_of(:pip, :pim, :pi0).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:pip, :pim, :pi0).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    chi2_cut 200                                                  # loose cut; tight cut in ROOT
  }

alg_pp0.with_decay_card(decay_pipipi0_ee).apply(sel_pp0)

### ------------------------------ Execution ------------------------------- ###
root_files_4g  = alg_4g.execute_on(energy_points + incmc_points +
                                   exMC_4gamma_ee + exMC_4gamma_mumu)

root_files_pp0 = alg_pp0.execute_on(energy_points + incmc_points +
                                    exMC_pipipi0_ee + exMC_pipipi0_mumu)