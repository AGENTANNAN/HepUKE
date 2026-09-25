### Dataset description ###
data_4918  = DatasetManager.real_data.find("707_4914")      # e+e- 4.918 GeV real data
incMC_4918 = DatasetManager.inclusive_mc.find("707_4914")   # matching inclusive MC
data_4951  = DatasetManager.real_data.find("707_4946")      # e+e- 4.9509 GeV real data (seeds the 4.9509 GeV MC)
incMC_4951 = DatasetManager.inclusive_mc.find("707_4946")   # 4.9509 GeV inclusive MC

# Signal decay card: e+e- -> Lambda_c+ Lambda_c(2595)- (+ c.c.)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c(2595)-   PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 p+ K- pi+   PHSP;
    Enddecay

    Decay anti-Lambda_c(2595)-
    1.0000 anti-Lambda_c- pi+ pi-   PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K+ pi-   PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC, one sample per energy point
exMCs_signal = DatasetManager.create_exclusive_mc_for([data_4918, data_4951]) do |config|
  config.sample_name   = "exmc_lamc_lamc2595"
  config.events        = 200000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name_4918 = "LamC2595_4918"
alg_4918 = Algorithm.new(alg_name_4918)
alg_4918.set_header(["#{alg_name_4918}Alg/#{alg_name_4918}.h"])
        .set_constant({"ECMS" => [:double, 4.918]})

alg_name_4951 = "LamC2595_4951"
alg_4951 = Algorithm.new(alg_name_4951)
alg_4951.set_header(["#{alg_name_4951}Alg/#{alg_name_4951}.h"])
        .set_constant({"ECMS" => [:double, 4.9509]})

# Common selection chain: charged tracks -> PID -> partial reconstruction of Lambda_c+
event_selection = Selection.new
    .select_track {
        cos_theta 0.93      # |cos(theta)| < 0.93
        Vz        10.0      # |Vz| < 10 cm
        Vr        1.0       # |Vxy| < 1 cm
        nTot      ">=3"     # at least three charged tracks
        nNet      "==0"     # net charge zero
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.6
        identify :prp, against: [:kaon, :pion]     # p+, veto K+ and pi+
        identify :km,  against: [:pion, :proton]   # K-, veto pi- and p-
        identify :pip, against: [:kaon, :proton]   # pi+, veto K+ and p+
        nprp ">=1"
        nkm  ">=1"
        npip ">=1"
    }
    # No kinematic fit: partial reconstruction combines p K- pi+ into the Lambda_c+
    # and infers the excited Lambda_c(2595)- from the recoil four-momentum.
    .partial_rec([1, 3, 4, 5]) do
        best_combination_by_mass :Lambda_c, 2.2865   # Lambda_c+ nominal mass
        require_recoil_mass 2.55, 2.75               # M_rec > 2.55 GeV/c2
    end

alg_4918
    .note(:lambda_c_mass_window, "M(p K- pi+) required in (2.27, 2.30) GeV/c2 to select the Lambda_c+ candidate")
    .note(:excited_lambda_c_isospin, "Lambda_c(2595)- -> anti-Lambda_c- pi0 pi0 mode present with isospin pi+pi-:pi0pi0 = 2:1 (varied 1.5:1-5:1 for systematics)")
    .with_decay_card(decay_card_signal)
    .apply(event_selection.dup)

alg_4951
    .note(:lambda_c_mass_window, "M(p K- pi+) required in (2.27, 2.30) GeV/c2 to select the Lambda_c+ candidate")
    .note(:excited_lambda_c_isospin, "Lambda_c(2595)- -> anti-Lambda_c- pi0 pi0 mode present with isospin pi+pi-:pi0pi0 = 2:1 (varied 1.5:1-5:1 for systematics)")
    .with_decay_card(decay_card_signal)
    .apply(event_selection.dup)

files_4918 = alg_4918.execute_on([data_4918, incMC_4918, exMCs_signal[0]])
files_4951 = alg_4951.execute_on([incMC_4951, exMCs_signal[1]])