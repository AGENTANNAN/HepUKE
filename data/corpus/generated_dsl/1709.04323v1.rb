# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset preparation ###
# Real data samples at the eight center-of-mass energies (3.686 - 4.600 GeV)
data_3686 = DatasetManager.real_data.find("709_3686")   # psi(3686), 3.686 GeV
data_3773 = DatasetManager.real_data.find("712_3773")   # 3.773 GeV
data_4009 = DatasetManager.real_data.find("703_4009")   # ~4.008 GeV
data_4230 = DatasetManager.real_data.find("703_4230")   # ~4.226 GeV
data_4260 = DatasetManager.real_data.find("703_4260")   # ~4.258 GeV
data_4360 = DatasetManager.real_data.find("703_4360")   # ~4.358 GeV
data_4420 = DatasetManager.real_data.find("703_4420")   # ~4.416 GeV
data_4600 = DatasetManager.real_data.find("703_4600")   # ~4.600 GeV

# Corresponding inclusive MC samples
incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360")
incMC_4420 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

data_points  = [data_3686, data_3773, data_4009, data_4230, data_4260, data_4360, data_4420, data_4600]
incMC_points = [incMC_3686, incMC_3773, incMC_4009, incMC_4230, incMC_4260, incMC_4360, incMC_4420, incMC_4600]

# Decay card for e+e- -> eta Y(2175), Y(2175) -> phi f0(980),
# f0(980) -> pi+ pi-, phi -> K+ K-, eta -> gamma gamma
# (Direct production with no intermediate resonance -> KKMC top mother psi(4260))
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta Y(2175) PHSP;
    Enddecay

    Decay Y(2175)
    1.000 phi f0(980) PHSP;
    Enddecay

    Decay f0(980)
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive signal MC at each energy point (shared decay card)
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_eta_Y2175_phi_f0"
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "EtaY2175PhiF0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # nominal beam energy (per-run value read from DB)
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93    # |cos(theta)| < 0.93
                  Vz        10.0    # |Vz| < 10 cm
                  Vr        1.0     # Vr < 1 cm
                  nChrp     "==2"   # exactly four tracks ...
                  nChrn     "==2"   # ... (2 positive + 2 negative)
                  nNet      "==0"   # net charge zero
                }
               .select_photon {
                  tdc_emc_start     0      # EMC timing window
                  tdc_emc_end       14
                  angle_to_track    10.0   # angle to nearest charged track > 10 deg
                  energyThreshold_b 0.025  # barrel energy threshold 25 MeV
                  energyThreshold_e 0.050  # endcap energy threshold 50 MeV
                  nGam              ">=2"  # at least two photons
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :kaon, against: [:pion, :proton]   # K+ and K- (pi/K separation)
                  identify :pion, against: [:kaon, :proton]   # pi+ and pi-
                  nkp  "==1"
                  nkm  "==1"
                  npip "==1"
                  npim "==1"
                }
               # Reconstruct eta -> gamma gamma (1C Kalman fit constraining m(gamma gamma) to m(eta))
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                  chi2_cut 25
                  neta  ">=1"   # at least one eta candidate
                }
               # Nominal 4C kinematic fit to eta K+ K- pi+ pi-
               .kinematic_fit([:eta, :kp, :km, :pip, :pim]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 60
                }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC and per-energy-point exclusive signal MC
root_files = my_algorithm.execute_on(data_points + incMC_points + exMCs)