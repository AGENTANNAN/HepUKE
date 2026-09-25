### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

# Mode 1 decay card: J/psi -> Ds- rho+, Ds- -> phi e- nu_e, phi -> K+ K-, rho+ -> pi+ pi0, pi0 -> gamma gamma
decay_card_mode1 = <<~DECAYCARD
  Decay J/psi
  1.000 D_s- rho+ PHSP;
  Enddecay

  Decay D_s-
  1.000 phi e- anti-nu_e PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  Decay rho+
  1.000 pi+ pi0 VSS;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode 2 decay card: J/psi -> anti-D0 anti-K*0, anti-D0 -> K+ e- nu_e, anti-K*0 -> K- pi+
decay_card_mode2 = <<~DECAYCARD
  Decay J/psi
  1.000 anti-D0 anti-K*0 PHSP;
  Enddecay

  Decay anti-D0
  1.000 K+ e- anti-nu_e PHSP;
  Enddecay

  Decay anti-K*0
  1.000 K- pi+ VSS;
  Enddecay

  End
DECAYCARD

# Exclusive MC samples: 600k events each
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_ds_rho"
  config.related_dataset = jpsi_data
  config.events          = 600_000
  config.decay_card      = decay_card_mode1
  config.cross_section   = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_d0bar_kstar0bar"
  config.related_dataset = jpsi_data
  config.events          = 600_000
  config.decay_card      = decay_card_mode2
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# ---------------------------------------------------------------------------
# Mode 1: J/psi -> Ds- rho+, Ds- -> phi e- nu_e, phi -> K+ K-, rho+ -> pi+ pi0, pi0 -> gamma gamma
# ---------------------------------------------------------------------------
alg_name_mode1 = "DsRhoMode"
alg_mode1 = Algorithm.new(alg_name_mode1)
alg_mode1.set_header(["#{alg_name_mode1}Alg/#{alg_name_mode1}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})

selection_mode1 = Selection.new
selection_mode1
  .select_track {                     # four charged tracks, net charge zero
    cos_theta 0.93                    # |cos(theta)| < 0.93
    Vz        20.0                    # |Vz| < 20 cm
    Vr        2.0                     # Vr < 2 cm
    nChrp     "==2"                   # K+ and pi+
    nChrn     "==2"                   # K- and e-
    nNet      "==0"                   # net charge zero
  }
  .select_photon {                    # photon selection
    tdc_emc_start     0               # EMC timing 0-700 ns
    tdc_emc_end       14
    angle_to_track    20.0            # more than 20 deg from any charged track
    energyThreshold_b 0.025           # 25 MeV barrel
    energyThreshold_e 0.050           # 50 MeV endcap
    nGam              ">=2"           # at least two photons (needed for pi0)
  }
  .pid(method: :probability) {        # PID, probability method
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6  # p>1.0 -> lepton; EMC E>0.6 -> electron
    identify :kaon, against: [:pion, :proton]   # K+ and K- separated from pi and p
    identify :pion, against: [:kaon, :proton]   # pi+ and pi- separated from K and p
    nkp  ">=1"                        # at least one K+
    nkm  ">=1"                        # at least one K-
    npip ">=1"                        # at least one pi+
    nlm  ">=1"                        # at least one e-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct pi0 from photon pairs
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  }
  .kinematic_fit([:kp, :km, :pip, :pi0, :lm]) {   # K+ K- pi+ pi0 e- with missing nu
    nominal
    miss_track_of :nu                 # undetected neutrino
    constrain_four_momentum           # 4C constraint
    invariant_mass_of(:kp, :km).between(1.01, 1.03)    # phi mass window
    invariant_mass_of(:pip, :pi0).between(0.62, 0.95)  # rho+ mass window
    chi2_cut 200
  }

alg_mode1
  .note(:electron_selection, "electron candidates additionally required to satisfy E/p > 0.8 and |cos theta| < 0.8; not expressible in the current DSL selection blocks")
  .note(:missing_momentum, "missing momentum > 0.1 GeV/c and |Umiss| < 0.05 GeV required for mode 1; computed from the undetected neutrino four-momentum outside the DSL")
  .with_decay_card(decay_card_mode1)
  .apply(selection_mode1)

# ---------------------------------------------------------------------------
# Mode 2: J/psi -> anti-D0 anti-K*0, anti-D0 -> K+ e- nu_e, anti-K*0 -> K- pi+
# ---------------------------------------------------------------------------
alg_name_mode2 = "D0barKstar0barMode"
alg_mode2 = Algorithm.new(alg_name_mode2)
alg_mode2.set_header(["#{alg_name_mode2}Alg/#{alg_name_mode2}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})

selection_mode2 = Selection.new
selection_mode2
  .select_track {                     # four charged tracks, net charge zero
    cos_theta 0.93
    Vz        20.0
    Vr        2.0
    nChrp     "==2"                   # K+ and pi+
    nChrn     "==2"                   # K- and e-
    nNet      "==0"
  }
  .select_photon {                    # photon selection (used for the pi0 veto)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {        # PID, probability method
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp  ">=1"
    nkm  ">=1"
    npip ">=1"
    nlm  ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # gamma gamma -> pi0 fit used as a pi0 veto (chi2 < 20)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 "==0"                        # veto events containing a pi0 candidate
  }
  .kinematic_fit([:kp, :km, :pip, :lm]) {     # K+ K- pi+ e- with missing nu
    nominal
    miss_track_of :nu                 # undetected neutrino
    constrain_four_momentum           # 4C constraint
    invariant_mass_of(:km, :pip).between(0.82, 0.98)   # anti-K*0 mass window
    chi2_cut 200
  }

alg_mode2
  .note(:electron_selection, "electron candidates additionally required to satisfy E/p > 0.8 and |cos theta| < 0.8; not expressible in the current DSL selection blocks")
  .note(:missing_momentum, "missing momentum > 0.1 GeV/c and |Umiss| < 0.02 GeV required for mode 2; computed from the undetected neutrino four-momentum outside the DSL")
  .with_decay_card(decay_card_mode2)
  .apply(selection_mode2)

### Execute ###
root_files_mode1 = alg_mode1.execute_on([jpsi_data, jpsi_incMC, exMC_mode1])
root_files_mode2 = alg_mode2.execute_on([jpsi_data, jpsi_incMC, exMC_mode2])