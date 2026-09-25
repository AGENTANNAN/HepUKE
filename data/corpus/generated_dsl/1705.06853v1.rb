# =============================================================================
# e+e- -> gamma eta_c(1S), 12 eta_c decay channels, at six c.m. energies
# (4.01, 4.23, 4.26, 4.36, 4.42, 4.60 GeV)
# =============================================================================

### -------------------- Datasets (six energy points) -------------------- ###
data_4010 = DatasetManager.real_data.find("703_4009")   # 4.01 GeV
data_4230 = DatasetManager.real_data.find("703_4230")   # 4.23 GeV
data_4260 = DatasetManager.real_data.find("703_4260")   # 4.26 GeV
data_4360 = DatasetManager.real_data.find("703_4360")   # 4.36 GeV
data_4420 = DatasetManager.real_data.find("703_4420")   # 4.42 GeV
data_4600 = DatasetManager.real_data.find("703_4600")   # 4.60 GeV
data_points = [data_4010, data_4230, data_4260, data_4360, data_4420, data_4600]

incMC_4010 = DatasetManager.inclusive_mc.find("703_4009")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360")
incMC_4420 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_points = [incMC_4010, incMC_4230, incMC_4260, incMC_4360, incMC_4420, incMC_4600]

### -------------------- Decay cards (12 eta_c channels) -------------------- ###

# Mode 1: pi+ pi- pi0 pi+ pi- pi0
dc_mode01 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi0 pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode 2: pi+ pi- pi0 pi0
dc_mode02 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi0 pi0 PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode 3: pi+ pi+ pi- pi- eta
dc_mode03 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi+ pi- pi- eta PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode 4: K+ K- pi+ pi- pi0
dc_mode04 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- pi+ pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode 5: pi+ pi- pi+ pi-
dc_mode05 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Mode 6: pi+ pi- pi+ pi- pi+ pi-
dc_mode06 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- pi+ pi- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Mode 7: pi+ pi- eta
dc_mode07 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 pi+ pi- eta PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode 8: K+- K_S pi-+ pi+ pi-
dc_mode08 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K_S0 pi- pi+ pi- PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Mode 9: K+- K_S pi-+
dc_mode09 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K_S0 pi- PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Mode 10: K+ K- pi0
dc_mode10 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- pi0 PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode 11: K+ K- pi+ pi-
dc_mode11 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Mode 12: K+ K- pi+ pi+ pi- pi-
dc_mode12 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma eta_c PHSP;
  Enddecay
  Decay eta_c
  1.0000 K+ K- pi+ pi+ pi- pi- PHSP;
  Enddecay
  End
DECAYCARD

### -------------------- Exclusive signal MC (100k events per mode/energy) -------------------- ###

exMCs_mode01 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etaC_gamma_2pi2pi02pi"
  config.events        = 100000
  config.decay_card    = dc_mode01
  config.cross_section = :default
end

exMCs_mode02 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etaC_gamma_2pi2pi0"
  config.events        = 100000
  config.decay_card    = dc_mode02
  config.cross_section = :default
end

exMCs_mode03 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etaC_gamma_2pi2pieta"
  config.events        = 100000
  config.decay_card    = dc_mode03
  config.cross_section = :default
end

exMCs_mode04 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etaC_gamma_KKpipipi0"
  config.events        = 100000
  config.decay_card    = dc_mode04
  config.cross_section = :default
end

exMCs_mode05 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etaC_gamma_4pi"
  config.events        = 100000
  config.decay_card    = dc_mode05
  config.cross_section = :default
end

exMCs_mode06 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etaC_gamma_6pi"
  config.events        = 100000
  config.decay_card    = dc_mode06
  config.cross_section = :default
end

exMCs_mode07 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etaC_gamma_pipieta"
  config.events        = 100000
  config.decay_card    = dc_mode07
  config.cross_section = :default
end

exMCs_mode08 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etaC_gamma_K_KS_3pi"
  config.events        = 100000
  config.decay_card    = dc_mode08
  config.cross_section = :default
end

exMCs_mode09 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etaC_gamma_K_KS_pi"
  config.events        = 100000
  config.decay_card    = dc_mode09
  config.cross_section = :default
end

exMCs_mode10 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etaC_gamma_KKpi0"
  config.events        = 100000
  config.decay_card    = dc_mode10
  config.cross_section = :default
end

exMCs_mode11 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etaC_gamma_KKpipi"
  config.events        = 100000
  config.decay_card    = dc_mode11
  config.cross_section = :default
end

exMCs_mode12 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "etaC_gamma_KK4pi"
  config.events        = 100000
  config.decay_card    = dc_mode12
  config.cross_section = :default
end

### -------------------- Event selection (BOSS), one algorithm per mode -------------------- ###

# --------------------------------------------------------------------------
# Mode 1 : eta_c -> pi+ pi- pi0 pi+ pi- pi0
# --------------------------------------------------------------------------
name1 = "EtaCGammaTo2pi2pi02pi"
alg_mode01 = Algorithm.new(name1)
alg_mode01.set_header(["#{name1}Alg/#{name1}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:background_veto, "the transition photon must not form a pi0 with any other EMC cluster")

sel_mode01 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=5"
  }
  .pid(method: :probability) {
    prob_cut 0.00001
    identify :pion, against: [:kaon]
    npip "==2"
    npim "==2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).within(0.107, 0.163)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim, :pi0, :pi0]) {
    nominal
    constrain_four_momentum     # 4C (+1C per pi0 from the mass-constrained pi0)
    chi2_cut 200
  }

alg_mode01.with_decay_card(dc_mode01).apply(sel_mode01)

# --------------------------------------------------------------------------
# Mode 2 : eta_c -> pi+ pi- pi0 pi0
# --------------------------------------------------------------------------
name2 = "EtaCGammaTo2pi2pi0"
alg_mode02 = Algorithm.new(name2)
alg_mode02.set_header(["#{name2}Alg/#{name2}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .note(:background_veto, "the transition photon must not form a pi0 with any other EMC cluster")

sel_mode02 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=5"
  }
  .pid(method: :probability) {
    prob_cut 0.00001
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).within(0.107, 0.163)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode02.with_decay_card(dc_mode02).apply(sel_mode02)

# --------------------------------------------------------------------------
# Mode 3 : eta_c -> pi+ pi+ pi- pi- eta
# --------------------------------------------------------------------------
name3 = "EtaCGammaTo2pi2pieta"
alg_mode03 = Algorithm.new(name3)
alg_mode03.set_header(["#{name3}Alg/#{name3}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .note(:background_veto, "the transition photon must not form a pi0 with any other EMC cluster")

sel_mode03 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"
  }
  .pid(method: :probability) {
    prob_cut 0.00001
    identify :pion, against: [:kaon]
    npip "==2"
    npim "==2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:gamma, :gamma).within(0.400, 0.700)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:gamma, :pip, :pip, :pim, :pim, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode03.with_decay_card(dc_mode03).apply(sel_mode03)

# --------------------------------------------------------------------------
# Mode 4 : eta_c -> K+ K- pi+ pi- pi0
# --------------------------------------------------------------------------
name4 = "EtaCGammaToKK2pipi0"
alg_mode04 = Algorithm.new(name4)
alg_mode04.set_header(["#{name4}Alg/#{name4}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .note(:background_veto, "the transition photon must not form a pi0 with any other EMC cluster")

sel_mode04 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"
  }
  .pid(method: :probability) {
    prob_cut 0.00001
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp  "==1"
    nkm  "==1"
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).within(0.107, 0.163)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :kp, :km, :pip, :pim, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode04.with_decay_card(dc_mode04).apply(sel_mode04)

# --------------------------------------------------------------------------
# Mode 5 : eta_c -> pi+ pi- pi+ pi-
# --------------------------------------------------------------------------
name5 = "EtaCGammaTo4pi"
alg_mode05 = Algorithm.new(name5)
alg_mode05.set_header(["#{name5}Alg/#{name5}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .note(:background_veto, "the transition photon must not form a pi0 with any other EMC cluster")

sel_mode05 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.00001
    identify :pion, against: [:kaon]
    npip "==2"
    npim "==2"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode05.with_decay_card(dc_mode05).apply(sel_mode05)

# --------------------------------------------------------------------------
# Mode 6 : eta_c -> pi+ pi- pi+ pi- pi+ pi-
# --------------------------------------------------------------------------
name6 = "EtaCGammaTo6pi"
alg_mode06 = Algorithm.new(name6)
alg_mode06.set_header(["#{name6}Alg/#{name6}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .note(:background_veto, "the transition photon must not form a pi0 with any other EMC cluster")

sel_mode06 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==3"
    nChrn     "==3"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.00001
    identify :pion, against: [:kaon]
    npip "==3"
    npim "==3"
  }
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode06.with_decay_card(dc_mode06).apply(sel_mode06)

# --------------------------------------------------------------------------
# Mode 7 : eta_c -> pi+ pi- eta
# --------------------------------------------------------------------------
name7 = "EtaCGammaTo2pieta"
alg_mode07 = Algorithm.new(name7)
alg_mode07.set_header(["#{name7}Alg/#{name7}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .note(:background_veto, "the transition photon must not form a pi0 with any other EMC cluster")
          .note(:background_veto, "in pi+ pi- eta the transition photon must be more than 17.5 deg away from all charged tracks")

sel_mode07 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"
  }
  .pid(method: :probability) {
    prob_cut 0.00001
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:gamma, :gamma).within(0.400, 0.700)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:gamma, :pip, :pim, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode07.with_decay_card(dc_mode07).apply(sel_mode07)

# --------------------------------------------------------------------------
# Mode 8 : eta_c -> K+- K_S pi-+ pi+ pi-
# --------------------------------------------------------------------------
name8 = "EtaCGammaToKKS3pi"
alg_mode08 = Algorithm.new(name8)
alg_mode08.set_header(["#{name8}Alg/#{name8}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .note(:background_veto, "the transition photon must not form a pi0 with any other EMC cluster")
          .note(:ks_daughter_vertex_exemption, "K_S daughters are exempt from the |Vz| and Vr requirements applied to prompt tracks; K_S is rebuilt from secondary_vertex_fit to pi+ pi- minimising the mass difference")

sel_mode08 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==3"
    nChrn     "==3"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.00001
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp "==1"
  }
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :kp, :pip, :pim, :pim, :K_S0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode08.with_decay_card(dc_mode08).apply(sel_mode08)

# --------------------------------------------------------------------------
# Mode 9 : eta_c -> K+- K_S pi-+
# --------------------------------------------------------------------------
name9 = "EtaCGammaToKKSp i".gsub(" ", "")
alg_mode09 = Algorithm.new(name9)
alg_mode09.set_header(["#{name9}Alg/#{name9}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .note(:background_veto, "the transition photon must not form a pi0 with any other EMC cluster")
          .note(:ks_daughter_vertex_exemption, "K_S daughters are exempt from the |Vz| and Vr requirements applied to prompt tracks; K_S is rebuilt from secondary_vertex_fit to pi+ pi- minimising the mass difference")

sel_mode09 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.00001
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp "==1"
  }
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :kp, :pim, :K_S0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode09.with_decay_card(dc_mode09).apply(sel_mode09)

# --------------------------------------------------------------------------
# Mode 10 : eta_c -> K+ K- pi0
# --------------------------------------------------------------------------
name10 = "EtaCGammaToKKpi0"
alg_mode10 = Algorithm.new(name10)
alg_mode10.set_header(["#{name10}Alg/#{name10}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .note(:background_veto, "the transition photon must not form a pi0 with any other EMC cluster")

sel_mode10 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"
  }
  .pid(method: :probability) {
    prob_cut 0.00001
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp "==1"
    nkm "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).within(0.107, 0.163)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :kp, :km, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode10.with_decay_card(dc_mode10).apply(sel_mode10)

# --------------------------------------------------------------------------
# Mode 11 : eta_c -> K+ K- pi+ pi-
# --------------------------------------------------------------------------
name11 = "EtaCGammaToKK2pi"
alg_mode11 = Algorithm.new(name11)
alg_mode11.set_header(["#{name11}Alg/#{name11}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .note(:background_veto, "the transition photon must not form a pi0 with any other EMC cluster")

sel_mode11 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.00001
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp  "==1"
    nkm  "==1"
    npip "==1"
    npim "==1"
  }
  .kinematic_fit([:gamma, :kp, :km, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode11.with_decay_card(dc_mode11).apply(sel_mode11)

# --------------------------------------------------------------------------
# Mode 12 : eta_c -> K+ K- pi+ pi+ pi- pi-
# --------------------------------------------------------------------------
name12 = "EtaCGammaToKK4pi"
alg_mode12 = Algorithm.new(name12)
alg_mode12.set_header(["#{name12}Alg/#{name12}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .note(:background_veto, "the transition photon must not form a pi0 with any other EMC cluster")

sel_mode12 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==3"
    nChrn     "==3"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.00001
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp  "==1"
    nkm  "==1"
    npip "==2"
    npim "==2"
  }
  .kinematic_fit([:gamma, :kp, :km, :pip, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_mode12.with_decay_card(dc_mode12).apply(sel_mode12)

### -------------------- Execute on data + inclusive MC + signal MC -------------------- ###
root_files = []
root_files += alg_mode01.execute_on(data_points + incMC_points + exMCs_mode01)
root_files += alg_mode02.execute_on(data_points + incMC_points + exMCs_mode02)
root_files += alg_mode03.execute_on(data_points + incMC_points + exMCs_mode03)
root_files += alg_mode04.execute_on(data_points + incMC_points + exMCs_mode04)
root_files += alg_mode05.execute_on(data_points + incMC_points + exMCs_mode05)
root_files += alg_mode06.execute_on(data_points + incMC_points + exMCs_mode06)
root_files += alg_mode07.execute_on(data_points + incMC_points + exMCs_mode07)
root_files += alg_mode08.execute_on(data_points + incMC_points + exMCs_mode08)
root_files += alg_mode09.execute_on(data_points + incMC_points + exMCs_mode09)
root_files += alg_mode10.execute_on(data_points + incMC_points + exMCs_mode10)
root_files += alg_mode11.execute_on(data_points + incMC_points + exMCs_mode11)
root_files += alg_mode12.execute_on(data_points + incMC_points + exMCs_mode12)