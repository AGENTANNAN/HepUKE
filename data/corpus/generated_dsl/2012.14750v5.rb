# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset preparation ###
# Centre-of-mass energy scan, 4130-4440 MeV (BOSS-705 round-12 samples)
scan_names = ["705_4130", "705_4160", "705_4290", "705_4315",
              "705_4340", "705_4380", "705_4400", "705_4440"]

data_points  = scan_names.map { |name| DatasetManager.real_data.find(name) }     # real data at each point
incMC_points = scan_names.map { |name| DatasetManager.inclusive_mc.find(name) }  # inclusive MC at each point

# Decay card for the signal process: psi(4260) -> mu+ mu- (phase space).
# psi(4260) is the top mother used with the KKMC generator (BESIII convention).
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0 mu+ mu- PHSP;
    Enddecay
    End
DECAYCARD

# One 1-million-event exclusive MC sample per energy point, all sharing the same decay card,
# cross section and event count (only the related dataset differs).
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "mumu_scan_exclusive_mc"
  config.events        = 1_000_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "MuMuScan"
mu_mu_alg = Algorithm.new(alg_name)
mu_mu_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.26]})   # nominal CMS energy constant (per-point value resolved from the dataset at run time)
         .set_alias({"std::vector<double>" => "Vdouble"})

# Common event selection applied at every energy point.
event_selection = Selection.new
event_selection
    .select_track {
        cos_theta 0.8    # |cos(theta)| < 0.8
        Vz 10.0          # |Vz| < 10 cm along the beam axis
        Vr 1.0           # Vr < 1 cm in the transverse plane
        nChrp "==1"      # exactly one positively charged track
        nChrn "==1"      # exactly one negatively charged track
        nNet  "==0"      # net charge equal to zero
    }
    .pid(method: :probability) {
        prob_cut 0.001   # PID probability > 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.6   # p > 1.0 GeV/c -> lepton; EMC energy > 0.6 GeV -> electron, else muon
        identify :pion, against: [:kaon]     # pion rejection (separate pions from kaons)
        nmup "==1"       # require one mu+
        nmum "==1"       # require one mu-
    }
    .kinematic_fit([:mup, :mum]) {
        nominal                  # nominal kinematic fit (only its corrected four-momenta are kept)
        constrain_four_momentum  # 4C constraint of the mu+ mu- system to the CMS energy
        chi2_cut 200             # loose chi^2 < 200 (tight cut determined later in ROOT)
    }

mu_mu_alg.with_decay_card(decay_card_signal).apply(event_selection)
root_files = mu_mu_alg.execute_on(data_points + incMC_points + exMCs)