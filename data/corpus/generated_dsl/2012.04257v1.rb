# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi real data at sqrt(s) = 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC sample

# Decay card for the signal mode: J/psi -> gamma eta', eta' -> pi+ pi- mu+ mu-
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta'          PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- mu+ mu-     PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the competing mode: J/psi -> gamma eta', eta' -> pi+ pi- pi+ pi-
decay_card_4pi = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta'          PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- pi+ pi-     PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for the signal mode
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gamma_etap_pipimumu"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# 500k-event exclusive MC for the competing four-pion mode (used for the 4pi veto)
exMC_4pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_gamma_etap_4pi"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_4pi
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "GammaEtaPrime"
my_algorithm = Algorithm.new(alg_name)                      # Algorithm for J/psi -> gamma eta'
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])   # Header generated from the decay card
            .set_constant({"ECMS" => [:double, 3.097]})     # sqrt(s) = 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {                              # Charged track selection
                  cos_theta 0.93      # |cos(theta)| < 0.93
                  Vz        10.0      # |Vz| < 10 cm
                  Vr        1.0       # Vr < 1 cm
                  nChrp     ">=2"     # at least two positive tracks
                  nChrn     ">=2"     # at least two negative tracks
                  nTot      "==4"     # four charged tracks in total
                  nNet      "==0"     # net charge zero
                }
               .select_photon {                             # Photon selection
                  tdc_emc_start     0      # TDC start time
                  tdc_emc_end       14     # TDC end time
                  angle_to_track    15.0   # > 15 degrees from the nearest charged track
                  energyThreshold_b 0.025  # 25 MeV in the EMC barrel
                  energyThreshold_e 0.050  # 50 MeV in the EMC endcap
                  nGam              ">=1"  # at least one photon
                }
               .pid(method: :probability) {                 # Probability-based PID
                  prob_cut 0.001                            # PID probability > 0.001
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.3,
                                                 treat_as_electron_if_energy_above: 0.35
                  identify :pion, against: [:kaon, :proton] # pi+ / pi- against K and p
                  nlp ">=1"                                 # at least one positive lepton (mu+)
                  nlm ">=1"                                 # at least one negative lepton (mu-)
                }
                # Nominal 4C kinematic fit to gamma pi+ pi- mu+ mu-
               .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {
                  nominal                                   # nominal fit: fitted four-momenta are saved
                  constrain_four_momentum                   # 4C energy-momentum constraint
                  chi2_cut 30                               # chi2 < 30
                }
               .assign({:chrgp => :pip, :chrgn => :pim})    # treat all charged tracks as pions
                # Competing 4C hypothesis gamma pi+ pi- pi+ pi- (no chi2 cut: chi2 stored for the ROOT-level 4pi veto)
               .kinematic_fit([:gamma, :pip, :pip, :pim, :pim]) {
                  constrain_four_momentum
                }

my_algorithm
  .note(:electron_veto, "the high-momentum lepton branch flags tracks with p > 0.3 GeV as leptons; electrons are rejected and muons are required to have EMC energy below 0.35 GeV (upper edge of the 0.05-0.35 GeV muon EMC window). The DSL lepton candidate lists merge e and mu, so the explicit 'no electron' requirement and the 0.05 GeV lower edge of the muon EMC-energy window are enforced at BOSS level outside the DSL expressions.")
  .with_decay_card(decay_card_signal)                        # Kinematic variables from the signal decay card
  .apply(event_selection)

# Execute the algorithm on real data, inclusive MC, and both exclusive MC samples
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_4pi])