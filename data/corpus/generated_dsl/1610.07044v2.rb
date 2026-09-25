# Core DSL classes and dependencies will be loaded automatically at execution

### ==================== Dataset description ==================== ###
# 17 high-luminosity XYZ c.m. energy points spanning 3.8962 - 4.5995 GeV (L > 40 pb^-1)
data_names = ["703_3900", "703_4009", "703_4090", "703_4180", "703_4190",
              "703_4200", "703_4210", "703_4220", "703_4230", "703_4237",
              "703_4246", "703_4260", "703_4270", "703_4280", "703_4360",
              "703_4420", "703_4600"]
data_points = data_names.map { |n| DatasetManager.real_data.find(n) }

# Inclusive MC exists only for a nine-point subset of the scan
incMC_names = ["703_4180", "703_4190", "703_4200", "703_4210", "703_4220",
               "703_4237", "703_4246", "703_4260", "703_4270"]
incMCs = incMC_names.map { |n| DatasetManager.inclusive_mc.find(n) }

### ==================== Decay cards (EvtGen format) ==================== ###
# Every signal is produced through e+e- -> pi+ pi- h_c, h_c -> gamma eta_c.
# The top mother particle is psi(4260) (BESIII / KKMC convention); the actual
# beam energy of each scan point is taken from the associated real dataset.

# Mode 01: eta_c -> p pbar
decay_card_mode01_ppbar = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 anti-p- p+ PHSP;
    Enddecay

    End
DECAYCARD

# Mode 02: eta_c -> 2(pi+ pi-)
decay_card_mode02_2pip2pim = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 pi+ pi- pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode 03: eta_c -> 2(K+ K-)
decay_card_mode03_2K2K = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 K+ K- K+ K- PHSP;
    Enddecay

    End
DECAYCARD

# Mode 04: eta_c -> pi+ pi- K+ K-
decay_card_mode04_pipimKK = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 pi+ pi- K+ K- PHSP;
    Enddecay

    End
DECAYCARD

# Mode 05: eta_c -> p pbar pi+ pi-
decay_card_mode05_ppbarpipim = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 anti-p- p+ pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode 06: eta_c -> 3(pi+ pi-)
decay_card_mode06_3pip3pim = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 pi+ pi- pi+ pi- pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode 07: eta_c -> 2(pi+ pi-) K+ K-
decay_card_mode07_2pip2pimKK = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 pi+ pi- pi+ pi- K+ K- PHSP;
    Enddecay

    End
DECAYCARD

# Mode 08: eta_c -> K_S0 K+- pi-+ (both charge conjugates)
decay_card_mode08_KsKpi = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    0.5 K_S0 K+ pi- PHSP;
    0.5 K_S0 K- pi+ PHSP;
    Enddecay

    Decay K_S0
    1.0 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode 09: eta_c -> K_S0 K+- pi-+ pi+ pi- (both charge conjugates)
decay_card_mode09_KsKpipipi = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    0.5 K_S0 K+ pi- pi+ pi- PHSP;
    0.5 K_S0 K- pi+ pi- pi+ PHSP;
    Enddecay

    Decay K_S0
    1.0 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode 10: eta_c -> K+ K- pi0
decay_card_mode10_KKpi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 K+ K- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 11: eta_c -> p pbar pi0
decay_card_mode11_ppbarpi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 anti-p- p+ pi0 PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 12: eta_c -> K+ K- eta
decay_card_mode12_KKeta = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 K+ K- eta PHSP;
    Enddecay

    Decay eta
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 13: eta_c -> pi+ pi- eta
decay_card_mode13_pipieta = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 14: eta_c -> 2(pi+ pi-) eta
decay_card_mode14_2pip2pimeta = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 pi+ pi- pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 15: eta_c -> pi+ pi- pi0 pi0
decay_card_mode15_pipimpi0pi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 pi+ pi- pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 16: eta_c -> 2(pi+ pi- pi0)
decay_card_mode16_2pip2pim2pi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0 pi+ pi- pi0 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### ============ Exclusive MC: 40k events per mode at every energy point ============ ###
exMC_mode01 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_ppbar"
  config.events = 40000
  config.decay_card = decay_card_mode01_ppbar
  config.cross_section = :default
end

exMC_mode02 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_2pip2pim"
  config.events = 40000
  config.decay_card = decay_card_mode02_2pip2pim
  config.cross_section = :default
end

exMC_mode03 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_2K2K"
  config.events = 40000
  config.decay_card = decay_card_mode03_2K2K
  config.cross_section = :default
end

exMC_mode04 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_pipimKK"
  config.events = 40000
  config.decay_card = decay_card_mode04_pipimKK
  config.cross_section = :default
end

exMC_mode05 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_ppbarpipim"
  config.events = 40000
  config.decay_card = decay_card_mode05_ppbarpipim
  config.cross_section = :default
end

exMC_mode06 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_3pip3pim"
  config.events = 40000
  config.decay_card = decay_card_mode06_3pip3pim
  config.cross_section = :default
end

exMC_mode07 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_2pip2pimKK"
  config.events = 40000
  config.decay_card = decay_card_mode07_2pip2pimKK
  config.cross_section = :default
end

exMC_mode08 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_KsKpi"
  config.events = 40000
  config.decay_card = decay_card_mode08_KsKpi
  config.cross_section = :default
end

exMC_mode09 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_KsKpipipi"
  config.events = 40000
  config.decay_card = decay_card_mode09_KsKpipipi
  config.cross_section = :default
end

exMC_mode10 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_KKpi0"
  config.events = 40000
  config.decay_card = decay_card_mode10_KKpi0
  config.cross_section = :default
end

exMC_mode11 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_ppbarpi0"
  config.events = 40000
  config.decay_card = decay_card_mode11_ppbarpi0
  config.cross_section = :default
end

exMC_mode12 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_KKeta"
  config.events = 40000
  config.decay_card = decay_card_mode12_KKeta
  config.cross_section = :default
end

exMC_mode13 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_pipieta"
  config.events = 40000
  config.decay_card = decay_card_mode13_pipieta
  config.cross_section = :default
end

exMC_mode14 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_2pip2pimeta"
  config.events = 40000
  config.decay_card = decay_card_mode14_2pip2pimeta
  config.cross_section = :default
end

exMC_mode15 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_pipimpi0pi0"
  config.events = 40000
  config.decay_card = decay_card_mode15_pipimpi0pi0
  config.cross_section = :default
end

exMC_mode16 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_hc_etac_2pip2pim2pi0"
  config.events = 40000
  config.decay_card = decay_card_mode16_2pip2pim2pi0
  config.cross_section = :default
end

### ==================== Event selection (BOSS) ==================== ###
# One Algorithm per eta_c decay mode (Rule T1).  Every chain shares the same
# charged-track / photon quality cuts and the same 0.001 probability PID, and
# ends in the 4C kinematic fit of e+e- -> pi+ pi- gamma X_i with chi2 < 200
# (the smallest-chi2 combination is kept by default).  The eta_c mass window
# and the h_c signal window quoted in the analysis are applied to
# kinematic-fit-corrected quantities and therefore belong to the ROOT stage.

### ---------------- Mode 01: eta_c -> p pbar ---------------- ###
alg01 = Algorithm.new("HcEtacPPbar")
alg01.set_header(["HcEtacPPbarAlg/HcEtacPPbar.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel01 = Selection.new
sel01.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=2"      # production pi+ and p
        nChrn ">=2"      # production pi- and pbar
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=1"       # prompt photon from h_c -> gamma eta_c
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :pion, against: [:kaon, :proton]
        nprp ">=1"
        nprm ">=1"
        npip ">=1"
        npim ">=1"
      }
      .kinematic_fit([:pip, :pim, :gamma, :prp, :prm]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg01.with_decay_card(decay_card_mode01_ppbar).apply(sel01)
alg01.execute_on(data_points + incMCs + exMC_mode01)

### ---------------- Mode 02: eta_c -> 2(pi+ pi-) ---------------- ###
alg02 = Algorithm.new("HcEtac2PiPi")
alg02.set_header(["HcEtac2PiPiAlg/HcEtac2PiPi.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel02 = Selection.new
sel02.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=3"
        nChrn ">=3"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=1"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon, :proton]
        npip ">=2"
        npim ">=2"
      }
      .kinematic_fit([:pip, :pim, :gamma, :pip, :pim, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg02.with_decay_card(decay_card_mode02_2pip2pim).apply(sel02)
alg02.execute_on(data_points + incMCs + exMC_mode02)

### ---------------- Mode 03: eta_c -> 2(K+ K-) ---------------- ###
alg03 = Algorithm.new("HcEtac2KK")
alg03.set_header(["HcEtac2KKAlg/HcEtac2KK.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel03 = Selection.new
sel03.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=3"
        nChrn ">=3"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=1"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]
        identify :pion, against: [:kaon, :proton]   # production pi+ pi-
        nkp ">=2"
        nkm ">=2"
        npip ">=1"
        npim ">=1"
      }
      .kinematic_fit([:kp, :km, :kp, :km, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg03.with_decay_card(decay_card_mode03_2K2K).apply(sel03)
alg03.execute_on(data_points + incMCs + exMC_mode03)

### ---------------- Mode 04: eta_c -> pi+ pi- K+ K- ---------------- ###
alg04 = Algorithm.new("HcEtacPiPiKK")
alg04.set_header(["HcEtacPiPiKKAlg/HcEtacPiPiKK.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel04 = Selection.new
sel04.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=3"
        nChrn ">=3"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=1"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]
        identify :pion, against: [:kaon, :proton]
        nkp ">=1"
        nkm ">=1"
        npip ">=1"
        npim ">=1"
      }
      .kinematic_fit([:pip, :pim, :kp, :km, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg04.with_decay_card(decay_card_mode04_pipimKK).apply(sel04)
alg04.execute_on(data_points + incMCs + exMC_mode04)

### ---------------- Mode 05: eta_c -> p pbar pi+ pi- ---------------- ###
alg05 = Algorithm.new("HcEtacPPbarPiPi")
alg05.set_header(["HcEtacPPbarPiPiAlg/HcEtacPPbarPiPi.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel05 = Selection.new
sel05.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=3"
        nChrn ">=3"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=1"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :pion, against: [:kaon, :proton]
        nprp ">=1"
        nprm ">=1"
        npip ">=1"
        npim ">=1"
      }
      .kinematic_fit([:prp, :prm, :pip, :pim, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg05.with_decay_card(decay_card_mode05_ppbarpipim).apply(sel05)
alg05.execute_on(data_points + incMCs + exMC_mode05)

### ---------------- Mode 06: eta_c -> 3(pi+ pi-) ---------------- ###
alg06 = Algorithm.new("HcEtac3PiPi")
alg06.set_header(["HcEtac3PiPiAlg/HcEtac3PiPi.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel06 = Selection.new
sel06.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=4"
        nChrn ">=4"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=1"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon, :proton]
        npip ">=3"
        npim ">=3"
      }
      .kinematic_fit([:pip, :pim, :gamma, :pip, :pim, :pip, :pim, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg06.with_decay_card(decay_card_mode06_3pip3pim).apply(sel06)
alg06.execute_on(data_points + incMCs + exMC_mode06)

### ---------------- Mode 07: eta_c -> 2(pi+ pi-) K+ K- ---------------- ###
alg07 = Algorithm.new("HcEtac2PiPiKK")
alg07.set_header(["HcEtac2PiPiKKAlg/HcEtac2PiPiKK.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel07 = Selection.new
sel07.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=4"
        nChrn ">=4"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=1"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]
        identify :pion, against: [:kaon, :proton]
        nkp ">=1"
        nkm ">=1"
        npip ">=2"
        npim ">=2"
      }
      .kinematic_fit([:kp, :km, :pip, :pim, :pip, :pim, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg07.with_decay_card(decay_card_mode07_2pip2pimKK).apply(sel07)
alg07.execute_on(data_points + incMCs + exMC_mode07)

### ---------------- Mode 08: eta_c -> K_S0 K+- pi-+ ---------------- ###
alg08 = Algorithm.new("HcEtacKsKPi")
alg08.set_header(["HcEtacKsKPiAlg/HcEtacKsKPi.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel08 = Selection.new
sel08.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=3"      # production pi+, K_S0 pi+ and K+/pi+
        nChrn ">=3"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=1"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]
        identify :pion, against: [:kaon, :proton]
        nkp ">=1"
        npip ">=2"
        npim ">=2"
      }
      # K_S0 -> pi+ pi- candidates from a secondary-vertex fit
      .secondary_vertex_fit([:pip, :pim]) {
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
      .kinematic_fit([:K_S0, :kp, :pim, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg08.note(:charge_conjugate_mode,
           "eta_c -> K_S0 K+- pi-+ contains both charge conjugates (K_S0 K+ pi- and
            K_S0 K- pi+); the track selection and PID are charge symmetric, while
            the 4C participant list shown here corresponds to the K+ pi- assignment")
alg08.with_decay_card(decay_card_mode08_KsKpi).apply(sel08)
alg08.execute_on(data_points + incMCs + exMC_mode08)

### ---------------- Mode 09: eta_c -> K_S0 K+- pi-+ pi+ pi- ---------------- ###
alg09 = Algorithm.new("HcEtacKsKPiPiPi")
alg09.set_header(["HcEtacKsKPiPiPiAlg/HcEtacKsKPiPiPi.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel09 = Selection.new
sel09.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=4"
        nChrn ">=4"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=1"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]
        identify :pion, against: [:kaon, :proton]
        nkp ">=1"
        npip ">=3"
        npim ">=3"
      }
      # K_S0 -> pi+ pi- candidates from a secondary-vertex fit
      .secondary_vertex_fit([:pip, :pim]) {
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
      .kinematic_fit([:K_S0, :kp, :pim, :pip, :pim, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg09.note(:charge_conjugate_mode,
           "eta_c -> K_S0 K+- pi-+ pi+ pi- contains both charge conjugates
            (K_S0 K+ pi- pi+ pi- and K_S0 K- pi+ pi- pi+); the 4C participant list
            shown here corresponds to the K+ pi- assignment")
alg09.with_decay_card(decay_card_mode09_KsKpipipi).apply(sel09)
alg09.execute_on(data_points + incMCs + exMC_mode09)

### ---------------- Mode 10: eta_c -> K+ K- pi0 ---------------- ###
alg10 = Algorithm.new("HcEtacKKPi0")
alg10.set_header(["HcEtacKKPi0Alg/HcEtacKKPi0.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel10 = Selection.new
sel10.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=2"
        nChrn ">=2"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=3"       # 2 photons from pi0 + the prompt h_c photon
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]
        identify :pion, against: [:kaon, :proton]   # production pi+ pi-
        nkp ">=1"
        nkm ">=1"
        npip ">=1"
        npim ">=1"
      }
      # pi0 -> gamma gamma candidates from a 1C mass-constrained fit
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      }
      .kinematic_fit([:kp, :km, :pi0, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg10.with_decay_card(decay_card_mode10_KKpi0).apply(sel10)
alg10.execute_on(data_points + incMCs + exMC_mode10)

### ---------------- Mode 11: eta_c -> p pbar pi0 ---------------- ###
alg11 = Algorithm.new("HcEtacPPbarPi0")
alg11.set_header(["HcEtacPPbarPi0Alg/HcEtacPPbarPi0.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel11 = Selection.new
sel11.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=2"
        nChrn ">=2"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=3"       # 2 photons from pi0 + the prompt h_c photon
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :pion, against: [:kaon, :proton]
        nprp ">=1"
        nprm ">=1"
        npip ">=1"
        npim ">=1"
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      }
      .kinematic_fit([:prp, :prm, :pi0, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg11.with_decay_card(decay_card_mode11_ppbarpi0).apply(sel11)
alg11.execute_on(data_points + incMCs + exMC_mode11)

### ---------------- Mode 12: eta_c -> K+ K- eta ---------------- ###
alg12 = Algorithm.new("HcEtacKKEta")
alg12.set_header(["HcEtacKKEtaAlg/HcEtacKKEta.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel12 = Selection.new
sel12.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=2"
        nChrn ">=2"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=3"       # 2 photons from eta + the prompt h_c photon
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]
        identify :pion, against: [:kaon, :proton]
        nkp ">=1"
        nkm ">=1"
        npip ">=1"
        npim ">=1"
      }
      # eta -> gamma gamma candidates from a 1C mass-constrained fit
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=1"
      }
      .kinematic_fit([:kp, :km, :eta, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg12.with_decay_card(decay_card_mode12_KKeta).apply(sel12)
alg12.execute_on(data_points + incMCs + exMC_mode12)

### ---------------- Mode 13: eta_c -> pi+ pi- eta ---------------- ###
alg13 = Algorithm.new("HcEtacPiPiEta")
alg13.set_header(["HcEtacPiPiEtaAlg/HcEtacPiPiEta.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel13 = Selection.new
sel13.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=2"
        nChrn ">=2"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=3"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon, :proton]
        npip ">=1"
        npim ">=1"
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=1"
      }
      .kinematic_fit([:pip, :pim, :eta, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg13.with_decay_card(decay_card_mode13_pipieta).apply(sel13)
alg13.execute_on(data_points + incMCs + exMC_mode13)

### ---------------- Mode 14: eta_c -> 2(pi+ pi-) eta ---------------- ###
alg14 = Algorithm.new("HcEtac2PiPiEta")
alg14.set_header(["HcEtac2PiPiEtaAlg/HcEtac2PiPiEta.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel14 = Selection.new
sel14.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=3"
        nChrn ">=3"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=3"
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon, :proton]
        npip ">=2"
        npim ">=2"
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=1"
      }
      .kinematic_fit([:pip, :pim, :eta, :gamma, :pip, :pim, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg14.with_decay_card(decay_card_mode14_2pip2pimeta).apply(sel14)
alg14.execute_on(data_points + incMCs + exMC_mode14)

### ---------------- Mode 15: eta_c -> pi+ pi- pi0 pi0 ---------------- ###
alg15 = Algorithm.new("HcEtacPiPiPi0Pi0")
alg15.set_header(["HcEtacPiPiPi0Pi0Alg/HcEtacPiPiPi0Pi0.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel15 = Selection.new
sel15.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=2"
        nChrn ">=2"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=5"       # 4 photons from the two pi0 + the prompt h_c photon
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon, :proton]
        npip ">=1"
        npim ">=1"
      }
      # two pi0 -> gamma gamma candidates from the 1C mass-constrained fit
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=2"
      }
      .kinematic_fit([:pip, :pim, :pi0, :pi0, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg15.with_decay_card(decay_card_mode15_pipimpi0pi0).apply(sel15)
alg15.execute_on(data_points + incMCs + exMC_mode15)

### ---------------- Mode 16: eta_c -> 2(pi+ pi- pi0) ---------------- ###
alg16 = Algorithm.new("HcEtac2PiPiPi0")
alg16.set_header(["HcEtac2PiPiPi0Alg/HcEtac2PiPiPi0.h"])
     .set_constant({"ECMS" => [:double, 4.26]})
sel16 = Selection.new
sel16.select_track {
        cos_theta 0.93
        Vz 10.0
        Vr 1.0
        nChrp ">=3"
        nChrn ">=3"
        nNet "==0"
      }
      .select_photon {
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam ">=5"       # 4 photons from the two pi0 + the prompt h_c photon
      }
      .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon, :proton]
        npip ">=2"
        npim ">=2"
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=2"
      }
      .kinematic_fit([:pip, :pim, :pi0, :pip, :pim, :pi0, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }
alg16.with_decay_card(decay_card_mode16_2pip2pim2pi0).apply(sel16)
alg16.execute_on(data_points + incMCs + exMC_mode16)