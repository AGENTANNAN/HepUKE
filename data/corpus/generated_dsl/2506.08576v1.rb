# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/ψ (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive MC sample

# Decay card for J/psi -> gamma eta', eta' -> pi+ pi- eta, eta -> gamma mu+ mu-
decay_card_mu = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma mu+ mu- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for J/psi -> gamma eta', eta' -> pi+ pi- eta, eta -> gamma e+ e-
decay_card_e = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples (1,000,000 events for each eta decay mode)
exMC_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "jpsi_gamma_etap_pipim_eta_gmumu"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_mu
  config.cross_section   = :default
end

exMC_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "jpsi_gamma_etap_pipim_eta_gee"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_e
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Common selection chain shared by both eta decay modes
# (final state gamma pi+ pi- l+ l-, with l = mu or e)
common_selection = Selection.new
  .select_track {                       # Charged track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     "==2"                     # Exactly two positively charged tracks
    nChrn     "==2"                     # Exactly two negatively charged tracks
    nNet      "==0"                     # Net charge zero
  }
  .select_photon {                      # Photon selection (requires >=2 photons)
    tdc_emc_start     0                 # TDC start time
    tdc_emc_end       14                # TDC end time
    angle_to_track    10.0              # Min 10 deg separation from any charged track
    energyThreshold_b 0.025             # 25 MeV barrel energy threshold
    energyThreshold_e 0.050             # 50 MeV endcap energy threshold
    nGam              ">=2"             # At least two photons
  }
  .pid(method: :probability) {          # PID by per-track probability, separating pi/mu/e
    prob_cut 0.001                      # Probability > 0.001
    # electron/muon cannot be separated by identify(): use the high-momentum lepton recipe
    # (fills index_lp / index_lm); a track with p > 1.0 GeV is a lepton, and a lepton with
    # EMC eraw > 0.6 GeV is an electron, otherwise a muon.
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:muon, :electron]   # pi+ and pi- against muons and electrons
    npip "==1"                          # One pi+
    npim "==1"                          # One pi-
    nlp  "==1"                          # One l+ (mu+ or e+)
    nlm  "==1"                          # One l- (mu- or e-)
  }

# Mode 1 (mu): eta -> gamma mu+ mu-; 5C fit = 4-momentum conservation + M(gamma l+l-) = m_eta
mu_selection = common_selection.dup
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :lp, :lm]) {
    nominal                                                     # Nominal fit (corrected 4-momenta used)
    constrain_four_momentum                                     # 4C energy-momentum constraint
    invariant_mass_of(:gamma, :lp, :lm).constrain_to_nominal_mass_of(:eta)  # 1C eta mass constraint
    chi2_cut 40                                                 # chi^2 < 40 for the mu mode
  }

# Mode 2 (e): eta -> gamma e+ e-; 5C fit = 4-momentum conservation + M(gamma l+l-) = m_eta
e_selection = common_selection.dup
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :lp, :lm]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:gamma, :lp, :lm).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200                                                # chi^2 < 200 for the e mode
  }
# NOTE: the full 6C fit (adding the eta' mass constraint) and the e-mode photon-conversion
# veto (2.0 cm < Rxy < 8.0 cm and cos(theta_e gamma) > 0.5) are applied in the ROOT analysis
# stage, i.e. after the kinematic fit, and are therefore out of scope here.

### Generate the algorithms ###
alg_mu = Algorithm.new("JpsiEtaPGmumu")
alg_mu.set_header(["JpsiEtaPGmumuAlg/JpsiEtaPGmumu.h"])
      .set_constant({"ECMS" => [:double, 3.097]})               # J/psi center-of-mass energy
      .set_alias({"std::vector<double>" => "Vdouble"})
alg_mu.with_decay_card(decay_card_mu).apply(mu_selection)

alg_e = Algorithm.new("JpsiEtaPGee")
alg_e.set_header(["JpsiEtaPGeeAlg/JpsiEtaPGee.h"])
     .set_constant({"ECMS" => [:double, 3.097]})
     .set_alias({"std::vector<double>" => "Vdouble"})
alg_e.with_decay_card(decay_card_e).apply(e_selection)

### Execute on datasets ###
root_files_mu = alg_mu.execute_on([jpsi_data, jpsi_incMC, exMC_mu])
root_files_e  = alg_e.execute_on([jpsi_data, jpsi_incMC, exMC_e])