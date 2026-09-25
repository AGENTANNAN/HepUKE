# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Real data at 4.009 GeV (BOSS 703, scan point 4009)
data_4009  = DatasetManager.real_data.find("703_4009")
# Corresponding inclusive MC sample at 4.009 GeV
incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")
# psi(3686) data control sample
psip_data  = DatasetManager.real_data.find("709_3686")

# Decay card: psi(4040) -> eta J/psi ; J/psi -> l+ l- (e / mu 50/50) ; eta -> gamma gamma
decay_card_eta_jpsi = <<~DECAYCARD
  Decay psi(4040)
  1.0000 eta J/psi PHSP;
  Enddecay

  Decay J/psi
  0.5000 mu+ mu- PHOTOS VLL;
  0.5000 e+  e-  PHOTOS VLL;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Decay card: psi(4040) -> pi0 J/psi ; J/psi -> mu+ mu- ; pi0 -> gamma gamma
decay_card_pi0_jpsi = <<~DECAYCARD
  Decay psi(4040)
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

# Exclusive MC for psi(4040) -> eta J/psi (200k events)
exMC_eta_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4009_etajpsi"
  config.related_dataset = data_4009
  config.events          = 200000
  config.decay_card      = decay_card_eta_jpsi
  config.cross_section   = :default
end

# Exclusive MC for psi(4040) -> pi0 J/psi (200k events)
exMC_pi0_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4009_pi0jpsi"
  config.related_dataset = data_4009
  config.events          = 200000
  config.decay_card      = decay_card_pi0_jpsi
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Two independent final states -> two Algorithm objects, one decay card each
alg_name_eta = "EtaJpsi"
alg_eta = Algorithm.new(alg_name_eta)
alg_eta.set_header(["#{alg_name_eta}Alg/#{alg_name_eta}.h"])
       .set_constant({"ECMS" => [:double, 4.009]})
       .set_alias({"std::vector<double>" => "Vdouble"})

alg_name_pi0 = "Pi0Jpsi"
alg_pi0 = Algorithm.new(alg_name_pi0)
alg_pi0.set_header(["#{alg_name_pi0}Alg/#{alg_name_pi0}.h"])
       .set_constant({"ECMS" => [:double, 4.009]})
       .set_alias({"std::vector<double>" => "Vdouble"})

# Selection chain shared by both final states (l+ l- gamma gamma)
event_selection_common = Selection.new
  .select_track {             # Charged track selection
    cos_theta 0.93            # |cos(theta)| < 0.93
    Vz        10.0            # |Vz| < 10 cm
    Vr        1.0             # Vr < 1 cm
    nChrp     "==1"           # exactly one positive track
    nChrn     "==1"           # exactly one negative track
    nNet      "==0"           # net charge zero
  }
  .select_photon {            # Photon selection
    tdc_emc_start     0       # TDC start
    tdc_emc_end       14      # TDC end
    angle_to_track    20.0    # at least 20 deg from any charged track
    energyThreshold_b 0.025   # E > 25 MeV in the barrel
    energyThreshold_e 0.050   # E > 50 MeV in the endcap
    nGam              "==2"   # exactly two good photons
  }
  .pid(method: :probability) {  # Lepton identification, probability method
    prob_cut 0.001                                       # PID probability > 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6  # p>1.0 -> lepton; EMC E>0.6 -> electron, else muon
    nlp "==1"                                            # one positive lepton
    nlm "==1"                                            # one negative lepton
  }

# eta J/psi: nominal 4C fit + stored 3C fit on the same l+ l- gamma gamma final state
eta_selection = event_selection_common.dup
  .kinematic_fit([:lp, :lm, :gamma, :gamma]) {
    nominal                    # nominal fit: four-momenta taken from this fit
    constrain_four_momentum    # 4C energy-momentum conservation
    chi2_cut 200               # loose chi2 < 200 (tight cut decided in ROOT)
  }
  .kinematic_fit([:lp, :lm, :gamma, :gamma]) {
    constrain_three_momentum   # 3C fit: only the chi2 is stored for the ROOT-level radiative-Bhabha/dimuon veto
  }

# pi0 J/psi: same selection chain as eta J/psi
pi0_selection = event_selection_common.dup
  .kinematic_fit([:lp, :lm, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  .kinematic_fit([:lp, :lm, :gamma, :gamma]) {
    constrain_three_momentum
  }

# eta J/psi: BOSS-side procedures that have no formal DSL construct
alg_eta
  .note(:recoil_mass_window, "the two-photon recoil mass M_recoil(gamma gamma) is required to lie in [2.9, 3.4] GeV/c^2 to select J/psi candidates")
  .note(:lepton_pid_criteria, "in addition to the momentum-based lepton assignment, a muon is defined by EMC deposit < 0.4 GeV and an electron by E/p > 0.8; the two lepton tracks are required to be of the same species")
  .note(:fsr_recovery, "photons inside a 5 deg cone around a lepton are added back to that lepton four-momentum (FSR / bremsstrahlung recovery) before the kinematic fit")
  .with_decay_card(decay_card_eta_jpsi).apply(eta_selection)

# pi0 J/psi: identical selection plus the extra muon-counter depth background veto
alg_pi0
  .note(:recoil_mass_window, "the two-photon recoil mass M_recoil(gamma gamma) is required to lie in [2.9, 3.4] GeV/c^2 to select J/psi candidates")
  .note(:lepton_pid_criteria, "in addition to the momentum-based lepton assignment, a muon is defined by EMC deposit < 0.4 GeV and an electron by E/p > 0.8; the two lepton tracks are required to be of the same species")
  .note(:fsr_recovery, "photons inside a 5 deg cone around a lepton are added back to that lepton four-momentum (FSR / bremsstrahlung recovery) before the kinematic fit")
  .note(:background_veto, "at least one charged track with muon-counter (MUC) hit depth > 30 cm is required to suppress the e+e- -> pi+pi-pi0 contamination")
  .with_decay_card(decay_card_pi0_jpsi).apply(pi0_selection)

# Execute on real data, inclusive MC, the psi(3686) control sample and the signal MC
root_files_eta = alg_eta.execute_on([data_4009, incMC_4009, psip_data, exMC_eta_jpsi])
root_files_pi0 = alg_pi0.execute_on([data_4009, incMC_4009, psip_data, exMC_pi0_jpsi])