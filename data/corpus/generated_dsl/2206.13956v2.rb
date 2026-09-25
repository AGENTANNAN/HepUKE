# =====================================================================
# Lepton-flavour-violating decay J/psi -> e mu
# Dataset preparation + BOSS event selection
# =====================================================================

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # corresponding inclusive J/psi MC

# Off-resonance continuum data used for background estimation
cont_3773 = DatasetManager.real_data.find("712_3773")             # continuum at 3.773 GeV
cont_3510 = DatasetManager.real_data.find("703_chi_c1_scan_4")    # continuum at ~3.510 GeV
cont_3080 = DatasetManager.real_data.find("708_3080")             # continuum at 3.080 GeV

# Decay card (EvtGen format) for the single signal mode J/psi -> e+ mu-
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 e+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Exclusive MC sample for the signal process (300k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_jpsi_emu"
  config.related_dataset = jpsi_data
  config.events          = 300000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "JpsiToEMu"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

# Counting experiment: no kinematic fit is performed on the BOSS side; the signal
# region is formed in ROOT.
my_algorithm.note(:no_kinematic_fit,
                  "counting experiment for J/psi -> e mu: no BOSS-side kinematic fit is run; " \
                  "the signal region is defined entirely in ROOT.")
# Lepton-PID probability cuts that have no dedicated DSL surface.
my_algorithm.note(:lepton_pid_selection,
                  "electron identified against mu/pi/K with PID probability > 0.8; muon identified " \
                  "against e/pi/K with PID probability > 0.001. Represented here by " \
                  "identify_high_momentum_leptons (the DSL-supported lepton-ID path).")

event_selection = Selection.new
event_selection.select_track {               # exactly two charged tracks
                  cos_theta  0.93            # |cos(theta)| < 0.93
                  Vz         10.0            # |Vz| < 10 cm
                  Vr         1.0             # Vr < 1 cm
                  nChrp      "==1"           # one positive track  (electron)
                  nChrn      "==1"           # one negative track  (muon)
                  nNet       "==0"           # net charge zero
                }
               .select_photon {              # veto any event containing a good photon
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    20.0     # opening angle > 20 deg
                  energyThreshold_b 0.025    # E > 25 MeV in the barrel (|cos(theta)| < 0.80)
                  energyThreshold_e 0.050    # endcap region (0.86 < |cos(theta)| < 0.92)
                  nGam              "==0"    # no photon allowed -> veto
                }
               .pid(method: :probability) {  # lepton identification (one e and one mu)
                  prob_cut 0.001
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  nlp "==1"                  # one positive lepton
                  nlm "==1"                  # one negative lepton
                }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC,
                                      cont_3773, cont_3510, cont_3080,
                                      exMC_signal])