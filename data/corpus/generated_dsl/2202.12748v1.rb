### Dataset description ###
# Real data at the energy points spanning 3.872 - 4.700 GeV (BOSS 703 XYZ scan)
data_3872 = DatasetManager.real_data.find("703_3872")
data_4180 = DatasetManager.real_data.find("703_4180")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4360 = DatasetManager.real_data.find("703_4360")
data_4420 = DatasetManager.real_data.find("703_4420")
data_4600 = DatasetManager.real_data.find("703_4600")

# Corresponding inclusive MC samples at the same energy points
incMC_3872 = DatasetManager.inclusive_mc.find("703_3872")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360")
incMC_4420 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

data_points   = [data_3872, data_4180, data_4260, data_4360, data_4420, data_4600]
incMC_samples = [incMC_3872, incMC_4180, incMC_4260, incMC_4360, incMC_4420, incMC_4600]

# Decay card for the signal process e+e- -> psi(4260) -> eta pi+ pi-, eta -> gamma gamma (EvtGen format)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta pi+ pi-    PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal, mapped to the same energy points as the real data
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "etapipi_signal_mc"
    config.events        = 100000
    config.decay_card    = decay_card_signal
    config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "EtaPiPi"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.26]})  # representative CMS energy of the scan

event_selection = Selection.new
event_selection.select_track {          # Exactly one positive and one negative charged track
                    cos_theta 0.93      # |cos(theta)| < 0.93
                    Vz        10.0      # |Vz| < 10 cm
                    Vr        1.0       # Vr < 1 cm
                    nChrp     "==1"     # exactly one positive track
                    nChrn     "==1"     # exactly one negative track
                    nNet      "==0"     # net charge zero
                }
                .select_photon {        # Exactly two good photons
                    tdc_emc_start     0
                    tdc_emc_end       14
                    angle_to_track    10.0    # photon-track opening angle > 10 degrees
                    energyThreshold_b 0.025   # barrel energy > 25 MeV
                    energyThreshold_e 0.050   # endcap energy > 50 MeV
                    nGam              "==2"   # exactly two photons
                }
                .pid(method: :probability) {   # Identify the two tracks as pi+ and pi-
                    prob_cut 0.001                        # PID probability > 0.001
                    identify :pion, against: [:kaon, :proton]  # pi+ and pi- vs K and p
                    npip "==1"
                    npim "==1"
                }
                .kalman_kinematic_fit([:gamma, :gamma]) {   # Reconstruct eta from the two photons
                    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)  # eta mass constraint
                    chi2_cut 25          # chi^2 < 25 for the Kalman fit
                    neta ">=1"           # at least one eta candidate
                }
                .kinematic_fit([:eta, :pip, :pim]) {   # 4C kinematic fit on eta pi+ pi-
                    nominal                 # mark as the nominal fit
                    constrain_four_momentum # constrain total four-momentum to the CMS energy
                    chi2_cut 200            # loose chi^2 cut (tight cut applied later in ROOT)
                }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on(data_points + incMC_samples + exMCs_signal)