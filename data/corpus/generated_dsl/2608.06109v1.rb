### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi real data at 3.097 GeV (~9e9 events)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive MC sample

# Decay card for J/psi -> gamma eta', eta' -> e mu (charge conjugate included)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'          PHSP;
    Enddecay

    Decay eta'
    0.5000 e+ mu-              PHSP;
    0.5000 e- mu+              PHSP;
    Enddecay

    End
DECAYCARD

# Generate 600k-event exclusive MC for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_gamma_etap_emu"
  config.related_dataset = jpsi_data      # Matched to the real data sample
  config.events          = 600000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "GammaEtaPrimeEMu"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})      # CMS energy = J/psi mass (GeV)
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
    .select_track {                       # Charged track selection
        cos_theta 0.93                    # |cos(theta)| < 0.93
        Vz        10.0                    # |Vz| < 10 cm
        Vr        1.0                     # Vr < 1 cm in the transverse plane
        nChrp     "==1"                   # exactly one positive track
        nChrn     "==1"                   # exactly one negative track
        nNet      "==0"                   # net charge zero
    }
    .select_photon {                      # Photon selection
        tdc_emc_start     0               # shower time window start (units of 50 ns)
        tdc_emc_end       14              # shower time window end  -> [0, 700] ns
        angle_to_track    10.0            # angle to closest charged track > 10 degrees
        energyThreshold_b 0.025           # barrel energy threshold 25 MeV (|cos(theta)| < 0.80)
        energyThreshold_e 0.050           # endcap energy threshold 50 MeV (0.86 < |cos(theta)| < 0.92)
        nGam              ">=1"           # at least one photon
    }
    .pid(method: :probability) {          # PID: probability method; e/mu cannot be split by identify()
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.5   # p > 0.5 GeV/c -> lepton
        nlp "==1"                          # exactly one positive lepton
        nlm "==1"                          # exactly one negative lepton
    }
    .kinematic_fit([:gamma, :lp, :lm]) {   # 4C kinematic fit over (gamma, e, mu)
        nominal
        constrain_four_momentum
        chi2_cut 200                       # loose chi2 < 200 in BOSS; tightened to < 100 in ROOT
    }

# Capture BOSS-side procedures that the DSL cannot express formally
my_algorithm
  .note(:pid_selection, "electron candidates require P_e(e)/(P_e(e)+P_e(pi)+P_e(K)) > 0.8, "\
        "P_e(e) > 0.001, E_EMC/p > 0.8 together with dE/dx cuts; muon candidates require dE/dx, "\
        "0.1 < E_EMC < 0.3 GeV and MUC penetration depth. The DSL's identify_high_momentum_leptons "\
        "only exposes the momentum threshold (0.5 GeV/c); the electron/muon discrimination cuts are "\
        "fixed defaults in the generated code and must be edited to match these criteria.")
  .note(:radiative_photon_selection, "the highest-energy photon (~1.4 GeV) is taken as the radiative "\
        "photon entering the 4C fit; multi-photon candidate ranking is not expressible in the current "\
        "DSL and is applied manually. The photon recoil-mass window (0.932, 0.982) GeV/c^2 (5 sigma "\
        "around the eta' mass) is applied in the ROOT analysis, not here.")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC, and signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])