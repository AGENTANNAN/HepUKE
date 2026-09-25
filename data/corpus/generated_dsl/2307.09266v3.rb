# ============================================================================
# BOSS (dataset + event selection) DSL for the branching-fraction measurement
# of Lambda_c+ -> p eta and Lambda_c+ -> p omega
# with the single-tag, direct-reconstruction method
# over the seven energy points 4.600 - 4.699 GeV
# ============================================================================

### Dataset description ###
# Seven real-data points of the 4.600-4.699 GeV scan (BOSS releases 703 / 706)
data_points = [
  DatasetManager.real_data.find("703_4600"),   # 4.5995 GeV
  DatasetManager.real_data.find("706_4610"),   # 4.6119 GeV
  DatasetManager.real_data.find("706_4620"),   # 4.6280 GeV
  DatasetManager.real_data.find("706_4640"),   # 4.6409 GeV
  DatasetManager.real_data.find("706_4660"),   # 4.6612 GeV
  DatasetManager.real_data.find("706_4680"),   # 4.6819 GeV
  DatasetManager.real_data.find("706_4700")    # 4.6988 GeV
]

# Inclusive MC at 4.680 GeV
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")

### Decay cards (EvtGen format) ###
# Mode I: Lambda_c+ -> p eta, eta -> gamma gamma
decay_card_mode1 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 p+ eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: Lambda_c+ -> p eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_mode2 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 p+ eta PHSP;
    Enddecay

    Decay eta
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode III: Lambda_c+ -> p omega, omega -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_mode3 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 p+ omega PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive signal MC ###
# 100k events per mode; the same signal is generated at every one of the seven
# energy points via create_exclusive_mc_for (per-point MC needed by the scan).
exMC_mode1 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_LcToPEta_etagg"
  config.events        = 100_000
  config.decay_card    = decay_card_mode1
  config.cross_section = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_LcToPEta_eta3pi"
  config.events        = 100_000
  config.decay_card    = decay_card_mode2
  config.cross_section = :default
end

exMC_mode3 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_LcToPOmega_omega3pi"
  config.events        = 100_000
  config.decay_card    = decay_card_mode3
  config.cross_section = :default
end

### Mode I: Lambda_c+ -> p eta, eta -> gamma gamma ###
alg_name_1 = "LcToPEtaGG"
alg_mode1 = Algorithm.new(alg_name_1)
alg_mode1.set_header(["#{alg_name_1}Alg/#{alg_name_1}.h"])
         .set_constant({"ECMS" => [:double, 4.680]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_mode1 = Selection.new
    .select_track {
        cos_theta 0.93   # |cos(theta)| < 0.93
        Vz        10.0   # |Vz| < 10 cm
        Vr        1.0    # Vr < 1 cm
        nChrp     "==1"  # exactly one positive track (the proton)
        nChrn     "==0"  # no negative track
        nNet      "==1"  # net charge +1
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    8.0      # at least 8 degrees from any charged track
        energyThreshold_b 0.025    # 25 MeV in the EMC barrel
        energyThreshold_e 0.050    # 50 MeV in the EMC endcap
        nGam              ">=2"    # at least two photons
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:pion, :kaon]   # separate p from pi/K
        nprp     "==1"
    }
    .for_each(:prp) {           # extra requirement: proton Vr < 0.2 cm
        where { vr > 0.2 }
        remove
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit of the two photons to the eta mass
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 20
        neta     "==1"          # keep exactly one eta
    }
    .kinematic_fit([:prp, :eta]) {              # nominal 4C fit to p eta
        nominal
        constrain_four_momentum
        chi2_cut 200
    }

alg_mode1.with_decay_card(decay_card_mode1).apply(sel_mode1)

### Mode II: Lambda_c+ -> p eta, eta -> pi+ pi- pi0, pi0 -> gamma gamma ###
alg_name_2 = "LcToPEta3pi"
alg_mode2 = Algorithm.new(alg_name_2)
alg_mode2.set_header(["#{alg_name_2}Alg/#{alg_name_2}.h"])
         .set_constant({"ECMS" => [:double, 4.680]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_mode2 = Selection.new
    .select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"   # p and pi+
        nChrn     "==1"   # pi-
        nNet      "==1"   # net charge +1
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    8.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=2"
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:pion, :kaon]    # p vs pi/K
        identify :pion,   against: [:proton, :kaon]  # pi vs p/K
        nprp     "==1"
        npip     "==1"
        npim     "==1"
    }
    .for_each(:prp) {           # extra requirement: proton Vr < 0.2 cm
        where { vr > 0.2 }
        remove
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 mass constraint, chi2 < 50
        invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)   # M(gamma gamma) window
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 50
        npi0     "==1"
    }
    .kalman_kinematic_fit([:pip, :pim, :pi0]) {  # combine pi+ pi- pi0 into an eta (mass window)
        invariant_mass_of(:pip, :pim, :pi0).within(0.536, 0.560)
        neta     "==1"
    }
    .kinematic_fit([:prp, :eta]) {               # nominal 4C fit to p eta
        nominal
        constrain_four_momentum
        chi2_cut 200
    }

alg_mode2.with_decay_card(decay_card_mode2).apply(sel_mode2)

### Mode III: Lambda_c+ -> p omega, omega -> pi+ pi- pi0, pi0 -> gamma gamma ###
alg_name_3 = "LcToPOmega3pi"
alg_mode3 = Algorithm.new(alg_name_3)
alg_mode3.set_header(["#{alg_name_3}Alg/#{alg_name_3}.h"])
         .set_constant({"ECMS" => [:double, 4.680]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_mode3 = Selection.new
    .select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"   # p and pi+
        nChrn     "==1"   # pi-
        nNet      "==1"   # net charge +1
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    8.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=2"
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:pion, :kaon]
        identify :pion,   against: [:proton, :kaon]
        nprp     "==1"
        npip     "==1"
        npim     "==1"
    }
    .for_each(:prp) {           # extra requirement: proton Vr < 0.2 cm
        where { vr > 0.2 }
        remove
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # pi0 mass constraint, chi2 < 50
        invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)   # M(gamma gamma) window
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 50
        npi0     "==1"
    }
    .kalman_kinematic_fit([:pip, :pim, :pi0]) {  # combine pi+ pi- pi0 into an omega (mass window)
        invariant_mass_of(:pip, :pim, :pi0).within(0.750, 0.810)
        nomega   "==1"
    }
    .kinematic_fit([:prp, :omega]) {             # nominal 4C fit to p omega
        nominal
        constrain_four_momentum
        chi2_cut 200
    }

# Dalitz R < 0.9 on the omega -> pi+ pi- pi0 combination has no dedicated DSL construct
alg_mode3.note(:dalitz_cut, "omega -> pi+ pi- pi0 candidates are required to satisfy the omega Dalitz variable R < 0.9; this per-candidate cut on the 3-pion combination precedes the 4C kinematic fit")

alg_mode3.with_decay_card(decay_card_mode3).apply(sel_mode3)

### Execute on real data, inclusive MC and the signal MC of each mode ###
root_files_mode1 = alg_mode1.execute_on(data_points + [incMC_4680] + exMC_mode1)
root_files_mode2 = alg_mode2.execute_on(data_points + [incMC_4680] + exMC_mode2)
root_files_mode3 = alg_mode3.execute_on(data_points + [incMC_4680] + exMC_mode3)