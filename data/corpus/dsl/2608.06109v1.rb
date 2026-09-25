### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # (8998 +/- 40) x 10^6 J/psi events
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # ~9 billion J/psi inclusive MC

# Decay card for signal process: J/psi -> gamma eta', eta' -> e mu (charge-conjugate included)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000  gamma  eta'                          HELAMP 1.0 0.0 1.0 0.0;
    Enddecay

    Decay eta'
    0.5000  e+   mu-                             PHSP;
    0.5000  e-   mu+                             PHSP;
    Enddecay

    End
DECAYCARD

# Create exclusive MC sample for the signal
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gamma_etap_emu"
  config.related_dataset = jpsi_data
  config.events          = 600000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "EtapEmu"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93          # |cos(theta)| < 0.93
                  Vz        10.0          # |Vz| < 10 cm
                  Vr        1.0           # |Vr| < 1 cm
                  nChrp     "==1"         # exactly one positive track
                  nChrn     "==1"         # exactly one negative track
                  nNet      "==0"         # net charge zero
                }
               .select_photon {
                  energyThreshold_b 0.025 # >25 MeV in barrel (|cos theta| < 0.80)
                  energyThreshold_e 0.050 # >50 MeV in end cap (0.86 < |cos theta| < 0.92)
                  tdc_emc_start     0
                  tdc_emc_end       14    # shower time in [0, 700] ns
                  angle_to_track    10.0  # >10 deg to closest charged-track extrapolation
                  nGam              ">=1" # at least one good photon
                }
               .pid(method: :probability) {
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.5,
                                                 treat_as_electron_if_energy_above: 0.8
                  prob_cut 0.001         # P_e(e) > 0.001 for electron hypothesis
                  nlp "==1"              # exactly one positive lepton
                  nlm "==1"              # exactly one negative lepton
                }
               # 4C kinematic fit on gamma e mu system
               .kinematic_fit([:gamma, :lp, :lm]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200            # loose cut in BOSS; tight cut (chi2 < 100) applied in ROOT
                }

my_algorithm
  .note(:pid_correction_method,
        "Electron PID: P_e(e)/(P_e(e)+P_e(pi)+P_e(K)) > 0.8, P_e(e) > 0.001, " \
        "E_EMC/p > 0.8, and chi_dE/dx^e(e)+chi_dE/dx^e(pi) > 2. " \
        "Muon PID: chi_dE/dx^mu(mu)+chi_dE/dx^mu(e) < 0, " \
        "0.1 GeV < E_EMC^mu < 0.3 GeV, and a momentum-dependent MUC penetration " \
        "depth cut (only for P_mu > 0.5 GeV/c).")
  .note(:signal_photon,
        "The highest-energy photon (~1.4 GeV from J/psi -> gamma eta') is selected " \
        "as the radiative photon candidate for the 4C kinematic fit.")
  .note(:signal_region,
        "Photon recoil mass M_gamma^recoil in (0.932, 0.982) GeV/c^2 " \
        "(5-sigma window around nominal eta' mass) applied in ROOT.")

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
