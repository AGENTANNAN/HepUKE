### Dataset description ###
# J/psi(3.097 GeV) real data and inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal process J/psi -> gamma eta_c, eta_c -> omega omega,
# omega -> pi+ pi- pi0, pi0 -> gamma gamma (EvtGen format)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta_c    PHSP;
    Enddecay

    Decay eta_c
    1.000 omega omega    PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0    OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal decay chain (2M events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etac_omegamega"
  config.related_dataset = jpsi_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "JpsiGammaEtacOmegaOmega"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy 3.097 GeV

event_selection = Selection.new
  .select_track {                       # Charged track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm in the transverse plane
    nChrp     "==2"                     # exactly 2 positive tracks
    nChrn     "==2"                     # exactly 2 negative tracks
    nNet      "==0"                     # net charge zero
  }
  .select_photon {                      # Photon selection
    tdc_emc_start     0                 # TDC start time
    tdc_emc_end       14                # TDC end time
    angle_to_track    10.0              # min angle to nearest charged track (deg)
    energyThreshold_b 0.025            # barrel energy threshold 25 MeV
    energyThreshold_e 0.050            # endcap energy threshold 50 MeV
    nGam              ">=5"             # at least 5 photons
  }
  .pid(method: :probability) {          # PID (probability method)
    prob_cut 0.001                      # PID probability > 0.001
    identify :pion, against: [:kaon, :proton]   # identify pi+/pi- vs K and p
    npip "==2"                          # exactly 2 pi+
    npim "==2"                          # exactly 2 pi-
  }
  .remove([:pip <= :chrgp, :pim <= :chrgn])     # remove identified pions from charged track lists
  .kalman_kinematic_fit([:gamma, :gamma]) {      # 1C Kalman fit: reconstruct pi0 from photon pairs
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25                          # chi2 < 25
    npi0     ">=2"                       # at least two pi0 candidates
  }
  .kinematic_fit([:gamma, :pip, :pip, :pim, :pim, :pi0, :pi0]) {   # 6C fit: gamma pi+pi-pi0 pi+pi-pi0
    nominal                              # nominal fit (corrected four-momenta used downstream)
    constrain_four_momentum              # four-momentum constraint
    chi2_cut 200                         # loose chi2 cut; tight chi2_6C < 60 applied in ROOT
  }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])