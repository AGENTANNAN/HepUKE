### Dataset selection ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi real data at sqrt(s) = 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

### Decay cards (EvtGen format) ###
# Mode I: J/psi -> gamma eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta VSP_PWAVE;
    Enddecay

    Decay eta
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: J/psi -> gamma eta, eta -> pi0 pi0 pi0, pi0 -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta VSP_PWAVE;
    Enddecay

    Decay eta
    1.0000 pi0 pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples (1M events each) ###
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_eta_3pi_modeI"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_eta_3pi0_modeII"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### ---------------------------------------------------------------- ###
### Mode I: J/psi -> gamma eta, eta -> pi+ pi- pi0 (pi0 -> gamma gamma)
### ---------------------------------------------------------------- ###
alg_name_modeI = "JpsiGammaEtaModeI"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
  .select_track {                       # exactly two charged tracks, net charge zero
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {                      # at least three photons
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"
  }
  .pid(method: :probability) {          # one pi+ and one pi-, identified against kaons
    prob_cut 0.001
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit to reconstruct pi0 from a gamma gamma pair
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pi0]) {  # 6C nominal fit on gamma pi+ pi- pi0
    nominal
    constrain_four_momentum                     # 4C energy-momentum conservation
    invariant_mass_of(:pip, :pim, :pi0).constrain_to_nominal_mass_of(:eta)  # eta mass constraint
    chi2_cut 80
  }
  .kinematic_fit([:gamma, :gamma, :gamma]) {          # competing J/psi -> 3gamma hypothesis (chi2 stored for ROOT-level veto)
    constrain_four_momentum
  }
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma]) {  # competing J/psi -> 4gamma hypothesis (chi2 stored for ROOT-level veto)
    constrain_four_momentum
  }

alg_modeI
  .note(:radiative_photon_selection,
        "after the pi0 is reconstructed from one gamma gamma pair, the highest-energy remaining photon is taken as the radiative gamma")
  .note(:background_veto,
        "ROOT-level veto against J/psi -> 3gamma and J/psi -> 4gamma; competing non-nominal fits store their chi2 values and the final veto is applied in the ROOT analysis")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

root_files_modeI = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])

### ----------------------------------------------------------------- ###
### Mode II: J/psi -> gamma eta, eta -> pi0 pi0 pi0 (pi0 -> gamma gamma)
### ----------------------------------------------------------------- ###
alg_name_modeII = "JpsiGammaEtaModeII"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
  .select_track {                       # no charged tracks
    nChrp "==0"
    nChrn "==0"
    nNet  "==0"
  }
  .select_photon {                      # at least seven photons (same energy/timing/opening-angle criteria)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=7"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fits to build three pi0 from gamma gamma pairs
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=3"
  }
  .kinematic_fit([:gamma, :pi0, :pi0, :pi0]) {  # 8C nominal fit on gamma pi0 pi0 pi0
    nominal
    constrain_four_momentum                             # 4C energy-momentum conservation
    invariant_mass_of(:pi0, :pi0, :pi0).constrain_to_nominal_mass_of(:eta)  # eta mass constraint
    chi2_cut 70
  }

alg_modeII
  .note(:radiative_photon_selection,
        "after three pi0 are reconstructed from gamma gamma pairs, the highest-energy remaining photon is taken as the radiative gamma")
  .note(:pi0_decay_angle_cut,
        "each pi0 candidate is required to satisfy |cos(theta_decay)| < 0.95; cut applied on the pi0 decay angle")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

root_files_modeII = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])