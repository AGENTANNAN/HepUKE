# Paper: 2308.13725v1
# J/psi → gamma pi0(eta, eta') → gamma gamma gamma
# Pure photon analysis using J/psi data (708_3097), KKMC generator

### Dataset preparation ###
jpsi_data = DatasetManager.load_real_data.find("708_3097")
jpsi_incMC = DatasetManager.load_inclusive_mc.find("708_3097")

decay_card_for_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma pi0    HELAMP 1.0 0.0 1.0 0.0 0.0 0.0;
    Enddecay

    Decay pi0
    1.0000 gamma gamma  PHSP;
    Enddecay

    Decay J/psi
    1.0000 gamma eta    HELAMP 1.0 0.0 1.0 0.0 0.0 0.0;
    Enddecay

    Decay eta
    1.0000 gamma gamma  PHSP;
    Enddecay

    Decay J/psi
    1.0000 gamma eta'   HELAMP 1.0 0.0 1.0 0.0 0.0 0.0;
    Enddecay

    Decay eta'
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_gamma_P_to_3gamma"
  config.related_dataset = jpsi_data
  config.events = 2300000
  config.decay_card = decay_card_for_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
# Analysis selects J/psi → gamma pi0/eta/eta' with pi0/eta/eta' → gamma gamma
# Zero charged tracks required; all tracks vetoed by nTot cut

alg_name = "JpsiGammaPTo3Gamma"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
              cos_theta   0.93
              Vz   10.0
              Vr   1.0
              nTot   "==0"   # Zero charged tracks
            }
           .select_photon {
              tdc_emc_start   -500    # -500 ns EMC time difference
              tdc_emc_end     500    # +500 ns EMC time difference
              energyThreshold_b   0.080   # Min 80 MeV in barrel
              energyThreshold_e   0.080   # Min 80 MeV in endcap
              nGam   ">=3"      # At least 3 photons
            }
           # 4C kinematic fit: J/psi → gamma gamma gamma
           .kinematic_fit([:gamma, :gamma, :gamma]) {
              nominal
              constrain_four_momentum
              chi2_cut 50   # chi2_4C < 50
           }
           # Competing hypothesis: 2-gamma fit for e+e- → gamma gamma veto
           .kinematic_fit([:gamma, :gamma]) {
              constrain_four_momentum
              # In ROOT: require chi2_4C(3gamma) < chi2_4C(2gamma)
           }
           # Competing hypothesis: 4-gamma fit for pi0pi0 background veto
           .kinematic_fit([:gamma, :gamma, :gamma, :gamma]) {
              constrain_four_momentum
              # In ROOT: require chi2_4C(3gamma) < chi2_4C(4gamma) for pi0 case
           }

my_Algorithm.with_decay_card(decay_card_for_signal).apply(event_selection)
my_Algorithm.note(:qed_background, "QED background estimated from data at 3.080 GeV and normalized to J/psi luminosity")
my_Algorithm.note(:signal_extraction, "Signal yields extracted via unbinned ML fits to M_gammagamma distributions in pi0/eta/eta' mass regions")
my_Algorithm.note(:branching_fraction, "B(J/psi → gamma P) = N_obs / (N_Jpsi × B(P → gamma gamma) × epsilon)")

root_files = my_Algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])