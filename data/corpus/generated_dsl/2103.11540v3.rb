### Dataset description ###
# J/psi (3.097 GeV) real data and corresponding inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal process J/psi -> e- tau+ ; tau+ -> pi+ pi0 anti-nu_tau ; pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 e- tau+    VLL;
    Enddecay

    Decay tau+
    1.0000 pi+ pi0 anti-nu_tau    TAUHADNU;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal process (100k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "jpsi_to_e_tau"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiToETau"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # J/psi center-of-mass energy (GeV)

event_selection = Selection.new
    .select_track {                 # Charged track selection
        cos_theta 0.8               # |cos(theta)| < 0.8
        Vz        100.0             # |Vz| < 100 cm
        Vr        10.0              # Vr < 10 cm
        nChrp     "==1"             # exactly one positive track
        nChrn     "==1"             # exactly one negative track
        nNet      "==0"             # net charge zero
    }
    .select_photon {                # Photon selection
        tdc_emc_start     0         # TDC start
        tdc_emc_end       14        # TDC end
        angle_to_track    10.0      # angle to nearest charged track > 10 degrees
        energyThreshold_b 0.025     # 25 MeV in EMC barrel
        energyThreshold_e 0.050     # 50 MeV in EMC endcap
        nGam              ">=2"     # at least two photons (for pi0 -> gamma gamma)
    }
    .pid(method: :probability) {    # Probability-method PID
        prob_cut 0.001              # PID probability > 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.6  # p>1.0 GeV -> lepton; EMC E>0.6 GeV -> electron
        identify :pion, against: [:kaon, :proton]  # pi+ and pi- vs K and p
        nlp  ">=1"                  # at least one e+
        nlm  ">=1"                  # at least one e-
        npip ">=1"                  # at least one pi+
        npim ">=1"                  # at least one pi-
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # Reconstruct pi0 from two photons (mass constrained)
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25                 # chi2 < 25
        npi0 ">=1"                  # at least one pi0 candidate
    }
    .kinematic_fit([:lm, :pip, :pi0]) {         # Kinematic fit: e-, pi+, pi0 with missing anti-nu_tau
        nominal                     # nominal fit (four-momentum output taken from here)
        miss_track_of :nu_tau_bar   # anti-nu_tau allowed to be missing
        constrain_four_momentum     # constrain total four-momentum to CMS
        chi2_cut 200                # loose chi2 cut (tight cut applied in ROOT)
    }

# Generate the algorithm for the process defined in the decay card and apply the selection
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])