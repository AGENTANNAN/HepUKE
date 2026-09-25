### Dataset preparation ###
# J/psi (3.097 GeV) real data and its corresponding inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card — Mode I: J/psi -> gamma eta, eta -> e+ e- e+ e-
decay_card_modeI = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta PHSP;
  Enddecay

  Decay eta
  1.0000 e+ e- e+ e- PHSP;
  Enddecay

  End
DECAYCARD

# Decay card — Mode II: J/psi -> gamma eta', eta' -> pi+ pi- eta, eta -> e+ e- e+ e-
decay_card_modeII = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta' PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.0000 e+ e- e+ e- PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for each mode (2M events each)
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_eta_to_4e"
  config.related_dataset = jpsi_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etap_to_pipi_eta_to_4e"
  config.related_dataset = jpsi_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ---------- Mode I: J/psi -> gamma eta, eta -> e+ e- e+ e- ----------
alg_modeI = Algorithm.new("JpsiGammaEta4e")
alg_modeI.set_header(["JpsiGammaEta4eAlg/JpsiGammaEta4e.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
sel_modeI.select_track {           # charged track selection
    cos_theta 0.93                 # |cos(theta)| < 0.93
    Vz        10.0                 # |Vz| < 10 cm
    Vr        10.0                 # Vr < 10 cm
    nNet      "==0"                # net charge zero
    nChrp     "==2"                # 2 positive tracks (the two e+)
    nChrn     "==2"                # 2 negative tracks (the two e-)
  }
  .select_photon {                 # photon selection
    tdc_emc_start     0
    tdc_emc_end       14           # TDC 0-14
    angle_to_track    15.0         # > 15 deg from nearest charged track
    energyThreshold_b 0.025        # 25 MeV in barrel
    energyThreshold_e 0.050        # 50 MeV in endcap
    nGam              ">=1"        # at least one photon
  }
  .pid(method: :probability) {     # PID by probability method, electrons via lepton primitive
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
    nlp "==2"                      # two e+ (stored as leptons)
    nlm "==2"                      # two e- (stored as leptons)
  }
  # nominal 4C kinematic fit to gamma e+ e- e+ e-
  .kinematic_fit([:gamma, :lp, :lp, :lm, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 200                   # loose cut in BOSS; tight chi2_4C+PID < 60 applied later
  }

alg_modeI
  .note(:helix_correction,
    "track helix-parameter correction applied to all charged tracks before the 4C kinematic fit")
  .note(:pid_correction_method,
    "electron identification uses L(e) > L(pi) from the probability PID (prob cut 0.001); the DSL lepton primitive identify_high_momentum_leptons uses momentum/EMC-energy thresholds instead")
  .note(:background_veto,
    "photon-conversion veto: e+e- pairs with opening angle Phi_ee < 70 deg and conversion vertex 2 cm < R_xy < 8 cm are rejected; pi->e misidentification veto: e+e- pairs with 10 deg < theta_ee1 < 30 deg and 10 deg < theta_ee2 < 60 deg are rejected")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

# ---------- Mode II: J/psi -> gamma eta', eta' -> pi+ pi- eta, eta -> e+ e- e+ e- ----------
alg_modeII = Algorithm.new("JpsiGammaEtap4e")
alg_modeII.set_header(["JpsiGammaEtap4eAlg/JpsiGammaEtap4e.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
sel_modeII.select_track {          # charged track selection
    cos_theta 0.93                 # |cos(theta)| < 0.93
    Vz        10.0                 # |Vz| < 10 cm
    Vr        10.0                 # Vr < 10 cm
    nNet      "==0"                # net charge zero
    nChrp     "==3"                # 3 positive tracks (pi+ and two e+)
    nChrn     "==3"                # 3 negative tracks (pi- and two e-)
  }
  .select_photon {                 # photon selection (same as Mode I)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    15.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon, :proton]   # pi+ and pi- (identified against K/p)
    nlp  "==2"                     # two e+
    nlm  "==2"                     # two e-
    npip ">=1"                     # at least one pi+
    npim ">=1"                     # at least one pi-
  }
  # nominal 5C kinematic fit to gamma pi+ pi- e+ e- e+ e- with M(pi+pi-eta(4e)) = M(eta')
  .kinematic_fit([:gamma, :pip, :pim, :lp, :lp, :lm, :lm]) {
    nominal                        # best combination chosen by minimum chi2 (DSL default)
    constrain_four_momentum        # 4C energy-momentum constraint
    invariant_mass_of(:pip, :pim, :lp, :lp, :lm, :lm)
      .constrain_to_nominal_mass_of(:etap)   # +1C eta' mass -> 5C fit
    chi2_cut 200                   # loose cut in BOSS; tight chi2_5C < 60 applied later
  }

alg_modeII
  .note(:helix_correction,
    "track helix-parameter correction applied to all charged tracks before the 5C kinematic fit")
  .note(:pid_correction_method,
    "electron identification uses L(e) > L(pi) from the probability PID (prob cut 0.001); the DSL lepton primitive identify_high_momentum_leptons uses momentum/EMC-energy thresholds instead")
  .note(:background_veto,
    "photon-conversion veto: e+e- pairs with opening angle Phi_ee < 40 deg and conversion vertex 2 cm < R_xy < 8 cm are rejected; pi->e misidentification veto: e+e- pairs with 20 deg < theta_ee1 < 60 deg and 40 deg < theta_ee2 < 80 deg are rejected")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

### Execute on datasets ###
root_files_modeI  = alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])
root_files_modeII = alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])