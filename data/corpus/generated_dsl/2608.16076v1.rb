# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi(3097) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC

# Decay card for the signal process:
#   J/psi -> Sigma0 anti-Sigma0
#   anti-Sigma0 -> anti-Lambda0 gamma      (tag side)
#   Sigma0      -> Lambda0 e+ e-           (Dalitz, signal side)
#   Lambda0 / anti-Lambda0 -> p pi- / anti-p pi+
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 Sigma0 anti-Sigma0 PHSP;
    Enddecay

    Decay anti-Sigma0
    1.000 anti-Lambda0 gamma PHSP;
    Enddecay

    Decay Sigma0
    1.000 Lambda0 e+ e- PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Decay card for the peaking "double-radiative" background:
#   J/psi -> Sigma0 anti-Sigma0 with BOTH sigmas decaying radiatively (Sigma -> Lambda gamma)
decay_card_bkg = <<~DECAYCARD
    Decay J/psi
    1.000 Sigma0 anti-Sigma0 PHSP;
    Enddecay

    Decay anti-Sigma0
    1.000 anti-Lambda0 gamma PHSP;
    Enddecay

    Decay Sigma0
    1.000 Lambda0 gamma PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples (500k events each for signal and peaking background)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_sigma0_sigmabar0_dalitz"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

exMC_bkg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_sigma0_sigmabar0_double_radiative"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_bkg
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiSigma0Sigmabar0"
my_alg = Algorithm.new(alg_name)
my_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy of J/psi(3097)
       .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {                # Charged-track selection
    cos_theta 0.93               # |cos(theta)| < 0.93
    Vz        30.0               # |Vz| < 30 cm  (no Vr cut: long Lambda flight)
    nChrp     ">=3"              # at least 3 positive tracks
    nChrn     ">=3"              # at least 3 negative tracks
  }
  .select_photon {               # Photon selection
    tdc_emc_start     0          # TDC start (0 -> 700 ns window)
    tdc_emc_end       14
    angle_to_track    10.0       # angle to nearest charged track > 10 deg
    energyThreshold_b 0.025      # > 25 MeV in barrel
    energyThreshold_e 0.050      # > 50 MeV in endcap
    nGam              ">=1"      # at least one good photon
  }
  .pid(method: :probability) {   # PID, probability method
    prob_cut 0.001               # probability > 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.5,
                                   treat_as_electron_if_energy_above: 0.6
    identify :proton, against: [:kaon, :pion]   # p / p-bar
    identify :pion,   against: [:kaon]          # pi / K separation
    nprp "==1"                   # exactly 1 p
    nprm "==1"                   # exactly 1 p-bar
    nlp  "==1"                   # exactly 1 l+
    nlm  "==1"                   # exactly 1 l-
  }
  # The remaining (soft) tracks are the pi+ pi- from Lambda / Lambda_bar
  .remove([:prp <= :chrgp, :lp <= :chrgp, :prm <= :chrgn, :lm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  # Reconstruct Lambda -> p pi- by secondary-vertex fit (best mass difference)
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Reconstruct anti-Lambda -> anti-p pi+
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Nominal 4C kinematic fit over the full final state gamma p p-bar pi+ pi- l+ l-
  .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :lp, :lm]) {
    nominal
    constrain_four_momentum      # 4C constraint to the CMS four-momentum
    chi2_cut 200
  }

# Attach the signal decay card and generate the algorithm
my_alg.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC, signal MC and peaking-background MC
root_files = my_alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_bkg])