# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Six CM energy points: 4.600, 4.628, 4.641, 4.661, 4.682 and 4.698 GeV
data_4600 = DatasetManager.real_data.find("703_4600")   # 4.600 GeV
data_4628 = DatasetManager.real_data.find("706_4620")   # 4.628 GeV
data_4641 = DatasetManager.real_data.find("706_4640")   # 4.641 GeV
data_4661 = DatasetManager.real_data.find("706_4660")   # 4.661 GeV
data_4682 = DatasetManager.real_data.find("706_4680")   # 4.682 GeV
data_4698 = DatasetManager.real_data.find("706_4700")   # 4.698 GeV
data_points = [data_4600, data_4628, data_4641, data_4661, data_4682, data_4698]

# Corresponding inclusive MC samples
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4628 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4641 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4661 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4682 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4698 = DatasetManager.inclusive_mc.find("706_4700")
incMC_points = [incMC_4600, incMC_4628, incMC_4641, incMC_4661, incMC_4682, incMC_4698]

# Decay card (EvtGen format) for the signal process:
# psi(4260) -> Lambda_c+ anti-Lambda_c-  (charge conjugate included)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 p+ K+ K- PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K- K+ PHSP;
    Enddecay

    End
DECAYCARD

# Phase-space exclusive signal MC: 1,000,000 events per energy point
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_Lambdac_to_pKK"   # suffixed per energy point automatically
  config.events        = 1_000_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "LambdaCTopKK"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.600]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {          # charged track selection
                  cos_theta   0.93      # |cos(theta)| < 0.93
                  Vz          10.0      # |Vz| < 10 cm
                  Vr          1.0       # Vr < 1 cm
                  nChrp       ">=2"     # at least two positive tracks
                  nChrn       ">=1"     # at least one negative track
                }
               .pid(method: :probability) {                     # PID, probability method
                  prob_cut   0.0                                # zero probability cut
                  identify :proton, against: [:kaon, :pion]     # p / anti-p vs K, pi
                  identify :kaon,   against: [:pion, :proton]   # K+ / K- vs pi, p
                  identify :pion,   against: [:kaon, :proton]   # pi+ / pi- vs K, p
                  nprp  ">=1"                                   # at least one proton
                  nkp   ">=1"                                   # at least one K+
                  nkm   ">=1"                                   # at least one K-
                }
               # tighter proton vertex cut |Vr| < 0.2 cm
               .remove(:prp) { condition "Vr > 0.2" }
               # 4C kinematic fit of the p K+ K- system to the e+e- CM, chi2 < 200
               .kinematic_fit([:prp, :kp, :km]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }

# Generate the complete algorithm for the process defined in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute the same selection on real data, inclusive MC and exclusive signal MC
root_files = my_algorithm.execute_on(data_points + incMC_points + exMC_signal)