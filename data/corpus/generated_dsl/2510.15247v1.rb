### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # Real J/psi data at sqrt(s) = 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC sample

# Decay card for the signal process: J/psi -> gamma eta_c, eta_c -> p pbar
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0000 p+ anti-p- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the dominant background process: J/psi -> p pbar pi0
decay_card_bkg = <<~DECAYCARD
    Decay J/psi
    1.0000 p+ anti-p- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: signal J/psi -> gamma eta_c, eta_c -> p pbar (200k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gamma_etac_ppbar"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Exclusive MC: background J/psi -> p pbar pi0 (200k events)
exMC_bkg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_ppbar_pi0"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_bkg
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiGammaEtacToPPbar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy = 3.097 GeV

# Charged tracks -> photons -> PID -> isolated photon -> 4C kinematic fit
event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93    # |cos(theta)| < 0.93
                  Vz        10.0    # |Vz| < 10 cm
                  Vr        1.0     # Vr < 1 cm in the transverse plane
                  nChrp     "==1"   # exactly one positive charged track
                  nChrn     "==1"   # exactly one negative charged track
                  nNet      "==0"   # net charge zero
                }
               .select_photon {
                  tdc_emc_start     0      # EMC timing window start (0)
                  tdc_emc_end       14     # EMC timing window end (14 x 50 ns = 700 ns)
                  angle_to_track    10.0   # min angle to nearest charged track (deg)
                  energyThreshold_b 0.025  # E > 25 MeV in the barrel
                  energyThreshold_e 0.050  # E > 50 MeV in the endcap
                }
               .pid(method: :probability) {
                  prob_cut 0.001                              # PID probability > 0.001
                  identify :proton, against: [:pion, :kaon]   # p+ and anti-p- vs pi and K
                  nprp "==1"                                  # exactly one proton
                  nprm "==1"                                  # exactly one anti-proton
                }
               .select_isolated_photon {
                  angle_to_prp_track 20.0   # reject clusters within 20 deg of the proton
                  angle_to_prm_track 30.0   # reject clusters within 30 deg of the anti-proton
                  nGam ">=1"                # at least one photon survives isolation
                }
               # 4C kinematic fit to gamma p pbar; best (smallest chi2) combination chosen automatically
               .kinematic_fit([:gamma, :prp, :prm]) {
                  nominal                  # nominal fit -> corrected four-momenta are saved
                  constrain_four_momentum  # 4C energy-momentum conservation
                  chi2_cut 200             # retain events with chi2 < 200
                }

# Attach the decay card and render the selection into the BOSS algorithm
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Apply the same selection chain to data, inclusive MC and both exclusive MC samples
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_bkg])