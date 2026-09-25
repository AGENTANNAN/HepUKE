# Search for e+e- -> gamma eta_c(1S) at six c.m. energies between 4.01 and 4.60 GeV
# (arXiv:1705.06853v1)
#
# The analysis combines twelve eta_c(1S) decay channels; each channel is an
# independent Algorithm (Rule T1: different track/photon multiplicities and
# different kinematic-fit participant lists).
# The Born cross section is measured at six energy points, so ECMS is NOT set
# via set_constant (multi-energy analysis) and the signal MC is created with
# create_exclusive_mc_for over the six data sets.

### Dataset description ###
data_401 = DatasetManager.real_data.find("703_4009")   # 4.01 GeV  (4007.62 MeV)
data_423 = DatasetManager.real_data.find("703_4230")   # 4.23 GeV  (4226.26 MeV)
data_426 = DatasetManager.real_data.find("703_4260")   # 4.26 GeV  (4257.97 MeV)
data_436 = DatasetManager.real_data.find("703_4360")   # 4.36 GeV  (4358.26 MeV)
data_442 = DatasetManager.real_data.find("703_4420")   # 4.42 GeV  (4415.58 MeV)
data_460 = DatasetManager.real_data.find("703_4600")   # 4.60 GeV  (4599.53 MeV)

incMC_401 = DatasetManager.inclusive_mc.find("703_4009")
incMC_423 = DatasetManager.inclusive_mc.find("703_4230")
incMC_426 = DatasetManager.inclusive_mc.find("703_4260")
incMC_436 = DatasetManager.inclusive_mc.find("703_4360")
incMC_442 = DatasetManager.inclusive_mc.find("703_4420")
incMC_460 = DatasetManager.inclusive_mc.find("703_4600")

energy_points = [data_401, data_423, data_426, data_436, data_442, data_460]

### Decay cards — e+e- -> gamma eta_c(1S), eta_c -> X_i (EvtGen, PHSP) ###
decay_card_m1 = <<~DECAYCARD
    Decay psi(4260)
    1.000  gamma  eta_c   PHSP;
    Enddecay
    Decay eta_c
    1.000  pi+ pi- pi0 pi+ pi- pi0  PHSP;
    Enddecay
    Decay pi0
    1.000  gamma gamma  PHSP;
    Enddecay
    End
DECAYCARD

decay_card_m2 = <<~DECAYCARD
    Decay psi(4260)
    1.000  gamma  eta_c   PHSP;
    Enddecay
    Decay eta_c
    1.000  pi+ pi- pi0 pi0  PHSP;
    Enddecay
    Decay pi0
    1.000  gamma gamma  PHSP;
    Enddecay
    End
DECAYCARD

decay_card_m3 = <<~DECAYCARD
    Decay psi(4260)
    1.000  gamma  eta_c   PHSP;
    Enddecay
    Decay eta_c
    1.000  pi+ pi+ pi- pi- eta  PHSP;
    Enddecay
    Decay eta
    1.000  gamma gamma  PHSP;
    Enddecay
    End
DECAYCARD

decay_card_m4 = <<~DECAYCARD
    Decay psi(4260)
    1.000  gamma  eta_c   PHSP;
    Enddecay
    Decay eta_c
    1.000  K+ K- pi+ pi- pi0  PHSP;
    Enddecay
    Decay pi0
    1.000  gamma gamma  PHSP;
    Enddecay
    End
DECAYCARD

decay_card_m5 = <<~DECAYCARD
    Decay psi(4260)
    1.000  gamma  eta_c   PHSP;
    Enddecay
    Decay eta_c
    1.000  pi+ pi- pi+ pi-  PHSP;
    Enddecay
    End
DECAYCARD

decay_card_m6 = <<~DECAYCARD
    Decay psi(4260)
    1.000  gamma  eta_c   PHSP;
    Enddecay
    Decay eta_c
    1.000  pi+ pi- pi+ pi- pi+ pi-  PHSP;
    Enddecay
    End
DECAYCARD

decay_card_m7 = <<~DECAYCARD
    Decay psi(4260)
    1.000  gamma  eta_c   PHSP;
    Enddecay
    Decay eta_c
    1.000  pi+ pi- eta  PHSP;
    Enddecay
    Decay eta
    1.000  gamma gamma  PHSP;
    Enddecay
    End
DECAYCARD

decay_card_m8 = <<~DECAYCARD
    Decay psi(4260)
    1.000  gamma  eta_c   PHSP;
    Enddecay
    Decay eta_c
    1.000  K+ K_S0 pi- pi+ pi-  PHSP;
    Enddecay
    Decay K_S0
    1.000  pi+ pi-  PHSP;
    Enddecay
    End
DECAYCARD

decay_card_m9 = <<~DECAYCARD
    Decay psi(4260)
    1.000  gamma  eta_c   PHSP;
    Enddecay
    Decay eta_c
    1.000  K+ K_S0 pi-  PHSP;
    Enddecay
    Decay K_S0
    1.000  pi+ pi-  PHSP;
    Enddecay
    End
DECAYCARD

decay_card_m10 = <<~DECAYCARD
    Decay psi(4260)
    1.000  gamma  eta_c   PHSP;
    Enddecay
    Decay eta_c
    1.000  K+ K- pi0  PHSP;
    Enddecay
    Decay pi0
    1.000  gamma gamma  PHSP;
    Enddecay
    End
DECAYCARD

decay_card_m11 = <<~DECAYCARD
    Decay psi(4260)
    1.000  gamma  eta_c   PHSP;
    Enddecay
    Decay eta_c
    1.000  K+ K- pi+ pi-  PHSP;
    Enddecay
    End
DECAYCARD

decay_card_m12 = <<~DECAYCARD
    Decay psi(4260)
    1.000  gamma  eta_c   PHSP;
    Enddecay
    Decay eta_c
    1.000  K+ K- pi+ pi+ pi- pi-  PHSP;
    Enddecay
    End
DECAYCARD

### Exclusive signal MC — same cards at each of the six energy points ###
exMC_m1  = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gamma_etac_m1_2pipipim_pi0"
  config.events        = 100000
  config.decay_card    = decay_card_m1
  config.cross_section = :default
end

exMC_m2  = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gamma_etac_m2_pipim_pi0pi0"
  config.events        = 100000
  config.decay_card    = decay_card_m2
  config.cross_section = :default
end

exMC_m3  = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gamma_etac_m3_pipipipimim_eta"
  config.events        = 100000
  config.decay_card    = decay_card_m3
  config.cross_section = :default
end

exMC_m4  = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gamma_etac_m4_kkpipim_pi0"
  config.events        = 100000
  config.decay_card    = decay_card_m4
  config.cross_section = :default
end

exMC_m5  = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gamma_etac_m5_2pipi"
  config.events        = 100000
  config.decay_card    = decay_card_m5
  config.cross_section = :default
end

exMC_m6  = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gamma_etac_m6_3pipi"
  config.events        = 100000
  config.decay_card    = decay_card_m6
  config.cross_section = :default
end

exMC_m7  = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gamma_etac_m7_pipi_eta"
  config.events        = 100000
  config.decay_card    = decay_card_m7
  config.cross_section = :default
end

exMC_m8  = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gamma_etac_m8_ks_k_pipim"
  config.events        = 100000
  config.decay_card    = decay_card_m8
  config.cross_section = :default
end

exMC_m9  = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gamma_etac_m9_ks_kpi"
  config.events        = 100000
  config.decay_card    = decay_card_m9
  config.cross_section = :default
end

exMC_m10 = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gamma_etac_m10_kk_pi0"
  config.events        = 100000
  config.decay_card    = decay_card_m10
  config.cross_section = :default
end

exMC_m11 = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gamma_etac_m11_kkpipi"
  config.events        = 100000
  config.decay_card    = decay_card_m11
  config.cross_section = :default
end

exMC_m12 = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_gamma_etac_m12_kk2pipi"
  config.events        = 100000
  config.decay_card    = decay_card_m12
  config.cross_section = :default
end

signal_mc = exMC_m1 + exMC_m2 + exMC_m3 + exMC_m4 + exMC_m5 + exMC_m6 +
            exMC_m7 + exMC_m8 + exMC_m9 + exMC_m10 + exMC_m11 + exMC_m12

### Event selection (BOSS) ###

# ---------------------------------------------------------------------------
# Mode 1: eta_c -> pi+ pi- pi0 pi+ pi- pi0   (4 charged tracks + 2 pi0)
# ---------------------------------------------------------------------------
alg_m1 = Algorithm.new("GamEtacM1")
alg_m1.set_header(["GamEtacM1Alg/GamEtacM1.h"])
sel_m1 = Selection.new
sel_m1.select_track do
        cos_theta 0.93    # |cos(theta)| < 0.93
        Vz        10.0    # |Vz| < 10 cm
        Vr        1.0     # Vr < 1 cm
        nChrp     "==2"
        nChrn     "==2"
        nNet      "==0"
      end
      .select_photon do
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        angle_to_track    10.0
        nGam              ">=5"   # 1 transition photon + 4 photons from the two pi0
      end
      .pid(method: :probability) do
        prob_cut 0.00001   # P(pi) > 1e-5
        identify :pion,  against: [:kaon]   # P_pi > 1e-5
        identify :kaon,  against: [:pion]   # P_K > 1e-5 and P_K > P_pi
      end
      .kalman_kinematic_fit([:gamma, :gamma]) do
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=2"
      end
      .kinematic_fit([:gamma, :pip, :pip, :pim, :pim, :pi0, :pi0]) do
        nominal
        constrain_four_momentum   # 4C + 1C per pi0 candidate
        chi2_cut 200
      end
alg_m1.with_decay_card(decay_card_m1).apply(sel_m1)

# ---------------------------------------------------------------------------
# Mode 2: eta_c -> pi+ pi- pi0 pi0   (2 charged tracks + 2 pi0)
# ---------------------------------------------------------------------------
alg_m2 = Algorithm.new("GamEtacM2")
alg_m2.set_header(["GamEtacM2Alg/GamEtacM2.h"])
sel_m2 = Selection.new
sel_m2.select_track do
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==1"
        nChrn     "==1"
        nNet      "==0"
      end
      .select_photon do
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        angle_to_track    10.0
        nGam              ">=5"
      end
      .pid(method: :probability) do
        prob_cut 0.00001
        identify :pion, against: [:kaon]
        identify :kaon, against: [:pion]
      end
      .kalman_kinematic_fit([:gamma, :gamma]) do
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=2"
      end
      .kinematic_fit([:gamma, :pip, :pim, :pi0, :pi0]) do
        nominal
        constrain_four_momentum
        chi2_cut 200
      end
alg_m2.with_decay_card(decay_card_m2).apply(sel_m2)

# ---------------------------------------------------------------------------
# Mode 3: eta_c -> pi+ pi+ pi- pi- eta   (4 charged tracks + eta)
# ---------------------------------------------------------------------------
alg_m3 = Algorithm.new("GamEtacM3")
alg_m3.set_header(["GamEtacM3Alg/GamEtacM3.h"])
sel_m3 = Selection.new
sel_m3.select_track do
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"
        nChrn     "==2"
        nNet      "==0"
      end
      .select_photon do
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        angle_to_track    10.0
        nGam              ">=3"
      end
      .pid(method: :probability) do
        prob_cut 0.00001
        identify :pion, against: [:kaon]
        identify :kaon, against: [:pion]
      end
      .kalman_kinematic_fit([:gamma, :gamma]) do
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=1"
      end
      .kinematic_fit([:gamma, :pip, :pip, :pim, :pim, :eta]) do
        nominal
        constrain_four_momentum
        chi2_cut 200
      end
alg_m3.with_decay_card(decay_card_m3).apply(sel_m3)

# ---------------------------------------------------------------------------
# Mode 4: eta_c -> K+ K- pi+ pi- pi0   (4 charged tracks + pi0)
# ---------------------------------------------------------------------------
alg_m4 = Algorithm.new("GamEtacM4")
alg_m4.set_header(["GamEtacM4Alg/GamEtacM4.h"])
sel_m4 = Selection.new
sel_m4.select_track do
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"
        nChrn     "==2"
        nNet      "==0"
      end
      .select_photon do
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        angle_to_track    10.0
        nGam              ">=3"
      end
      .pid(method: :probability) do
        prob_cut 0.00001
        identify :pion, against: [:kaon]
        identify :kaon, against: [:pion]
      end
      .kalman_kinematic_fit([:gamma, :gamma]) do
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      end
      .kinematic_fit([:gamma, :kp, :km, :pip, :pim, :pi0]) do
        nominal
        constrain_four_momentum
        chi2_cut 200
      end
alg_m4.with_decay_card(decay_card_m4).apply(sel_m4)

# ---------------------------------------------------------------------------
# Mode 5: eta_c -> pi+ pi- pi+ pi-   (4 charged tracks)
# ---------------------------------------------------------------------------
alg_m5 = Algorithm.new("GamEtacM5")
alg_m5.set_header(["GamEtacM5Alg/GamEtacM5.h"])
sel_m5 = Selection.new
sel_m5.select_track do
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"
        nChrn     "==2"
        nNet      "==0"
      end
      .select_photon do
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        angle_to_track    10.0
        nGam              ">=1"
      end
      .pid(method: :probability) do
        prob_cut 0.00001
        identify :pion, against: [:kaon]
        identify :kaon, against: [:pion]
      end
      .kinematic_fit([:gamma, :pip, :pip, :pim, :pim]) do
        nominal
        constrain_four_momentum
        chi2_cut 200
      end
alg_m5.with_decay_card(decay_card_m5).apply(sel_m5)

# ---------------------------------------------------------------------------
# Mode 6: eta_c -> pi+ pi- pi+ pi- pi+ pi-   (6 charged tracks)
# ---------------------------------------------------------------------------
alg_m6 = Algorithm.new("GamEtacM6")
alg_m6.set_header(["GamEtacM6Alg/GamEtacM6.h"])
sel_m6 = Selection.new
sel_m6.select_track do
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==3"
        nChrn     "==3"
        nNet      "==0"
      end
      .select_photon do
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        angle_to_track    10.0
        nGam              ">=1"
      end
      .pid(method: :probability) do
        prob_cut 0.00001
        identify :pion, against: [:kaon]
        identify :kaon, against: [:pion]
      end
      .kinematic_fit([:gamma, :pip, :pip, :pip, :pim, :pim, :pim]) do
        nominal
        constrain_four_momentum
        chi2_cut 200
      end
alg_m6.with_decay_card(decay_card_m6).apply(sel_m6)

# ---------------------------------------------------------------------------
# Mode 7: eta_c -> pi+ pi- eta   (2 charged tracks + eta)
# ---------------------------------------------------------------------------
alg_m7 = Algorithm.new("GamEtacM7")
alg_m7.set_header(["GamEtacM7Alg/GamEtacM7.h"])
sel_m7 = Selection.new
sel_m7.select_track do
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==1"
        nChrn     "==1"
        nNet      "==0"
      end
      .select_photon do
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        angle_to_track    10.0
        nGam              ">=3"
      end
      .pid(method: :probability) do
        prob_cut 0.00001
        identify :pion, against: [:kaon]
        identify :kaon, against: [:pion]
      end
      .kalman_kinematic_fit([:gamma, :gamma]) do
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=1"
      end
      .kinematic_fit([:gamma, :pip, :pim, :eta]) do
        nominal
        constrain_four_momentum
        chi2_cut 200
      end
alg_m7.with_decay_card(decay_card_m7).apply(sel_m7)

# ---------------------------------------------------------------------------
# Mode 8: eta_c -> K+- K_S pi-+ pi+ pi-   (K_S from a secondary vertex)
# ---------------------------------------------------------------------------
alg_m8 = Algorithm.new("GamEtacM8")
alg_m8.set_header(["GamEtacM8Alg/GamEtacM8.h"])
sel_m8 = Selection.new
sel_m8.select_track do
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==3"
        nChrn     "==3"
        nNet      "==0"
      end
      .select_photon do
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        angle_to_track    10.0
        nGam              ">=1"
      end
      .pid(method: :probability) do
        prob_cut 0.00001
        identify :pion, against: [:kaon]
        identify :kaon, against: [:pion]
      end
      .secondary_vertex_fit([:pip, :pim]) do
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      end
      .kinematic_fit([:gamma, :kp, :K_S0, :pim, :pip, :pim]) do
        nominal
        constrain_four_momentum
        chi2_cut 200
      end
alg_m8.with_decay_card(decay_card_m8).apply(sel_m8)

# ---------------------------------------------------------------------------
# Mode 9: eta_c -> K+- K_S pi-+   (K_S from a secondary vertex)
# ---------------------------------------------------------------------------
alg_m9 = Algorithm.new("GamEtacM9")
alg_m9.set_header(["GamEtacM9Alg/GamEtacM9.h"])
sel_m9 = Selection.new
sel_m9.select_track do
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"
        nChrn     "==2"
        nNet      "==0"
      end
      .select_photon do
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        angle_to_track    10.0
        nGam              ">=1"
      end
      .pid(method: :probability) do
        prob_cut 0.00001
        identify :pion, against: [:kaon]
        identify :kaon, against: [:pion]
      end
      .secondary_vertex_fit([:pip, :pim]) do
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      end
      .kinematic_fit([:gamma, :kp, :K_S0, :pim]) do
        nominal
        constrain_four_momentum
        chi2_cut 200
      end
alg_m9.with_decay_card(decay_card_m9).apply(sel_m9)

# ---------------------------------------------------------------------------
# Mode 10: eta_c -> K+ K- pi0   (2 charged tracks + pi0)
# ---------------------------------------------------------------------------
alg_m10 = Algorithm.new("GamEtacM10")
alg_m10.set_header(["GamEtacM10Alg/GamEtacM10.h"])
sel_m10 = Selection.new
sel_m10.select_track do
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     "==1"
         nChrn     "==1"
         nNet      "==0"
       end
       .select_photon do
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=3"
       end
       .pid(method: :probability) do
         prob_cut 0.00001
         identify :pion, against: [:kaon]
         identify :kaon, against: [:pion]
       end
       .kalman_kinematic_fit([:gamma, :gamma]) do
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
         chi2_cut 25
         npi0 ">=1"
       end
       .kinematic_fit([:gamma, :kp, :km, :pi0]) do
         nominal
         constrain_four_momentum
         chi2_cut 200
       end
alg_m10.with_decay_card(decay_card_m10).apply(sel_m10)

# ---------------------------------------------------------------------------
# Mode 11: eta_c -> K+ K- pi+ pi-   (4 charged tracks)
# ---------------------------------------------------------------------------
alg_m11 = Algorithm.new("GamEtacM11")
alg_m11.set_header(["GamEtacM11Alg/GamEtacM11.h"])
sel_m11 = Selection.new
sel_m11.select_track do
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     "==2"
         nChrn     "==2"
         nNet      "==0"
       end
       .select_photon do
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=1"
       end
       .pid(method: :probability) do
         prob_cut 0.00001
         identify :pion, against: [:kaon]
         identify :kaon, against: [:pion]
       end
       .kinematic_fit([:gamma, :kp, :km, :pip, :pim]) do
         nominal
         constrain_four_momentum
         chi2_cut 200
       end
alg_m11.with_decay_card(decay_card_m11).apply(sel_m11)

# ---------------------------------------------------------------------------
# Mode 12: eta_c -> K+ K- pi+ pi+ pi- pi-   (6 charged tracks)
# ---------------------------------------------------------------------------
alg_m12 = Algorithm.new("GamEtacM12")
alg_m12.set_header(["GamEtacM12Alg/GamEtacM12.h"])
sel_m12 = Selection.new
sel_m12.select_track do
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     "==3"
         nChrn     "==3"
         nNet      "==0"
       end
       .select_photon do
         tdc_emc_start     0
         tdc_emc_end       14
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         nGam              ">=1"
       end
       .pid(method: :probability) do
         prob_cut 0.00001
         identify :pion, against: [:kaon]
         identify :kaon, against: [:pion]
       end
       .kinematic_fit([:gamma, :kp, :km, :pip, :pip, :pim, :pim]) do
         nominal
         constrain_four_momentum
         chi2_cut 200
       end
alg_m12.with_decay_card(decay_card_m12).apply(sel_m12)

### BOSS-side procedures that have no formal DSL construct ###
etac_algorithms = [alg_m1, alg_m2, alg_m3, alg_m4, alg_m5, alg_m6,
                   alg_m7, alg_m8, alg_m9, alg_m10, alg_m11, alg_m12]

etac_algorithms.each do |alg|
  alg.note(:background_veto,
           "the candidate transition photon must not be pairable with any other " \
           "energy cluster in the event to form a pi0 (pi0 veto against " \
           "pi0 -> gamma gamma fakes)")
  alg.note(:pid_correction_method,
           "pion PID requires P_pi > 1e-5; kaon PID requires P_K > 1e-5 and P_K > P_pi, " \
           "using combined dE/dx (MDC) and TOF information")
  alg.note(:track_vertex_exception,
           "pions originating from K_S decays are exempt from the |Vz| < 10 cm and " \
           "Vr < 1 cm interaction-point requirements; the K_S daughters are selected " \
           "through the secondary vertex fit")
  alg.note(:mass_window,
           "pi0 candidates: 0.107 < M(gamma gamma) < 0.163 GeV/c^2; " \
           "eta candidates: 0.400 < M(gamma gamma) < 0.700 GeV/c^2; " \
           "K_S candidates: 0.471 < M(pi+ pi-) < 0.524 GeV/c^2")
  alg.note(:chi2_optimisation,
           "the chi2/dof requirement of the 4C+1C kinematic fit is optimised " \
           "separately for each decay channel in the range 3.0-5.2; only the " \
           "combination with the best chi2/dof is retained (loose chi2_cut 200 " \
           "applied in BOSS, tight cut applied in ROOT)")
end

etac_algorithms.each do |alg|
  alg.note(:untagged_isr_selection,
           "the eta_c is reconstructed in the recoil-mass distribution of the " \
           "transition photon; the ISR photon from e+e- -> gamma_ISR J/psi is a " \
           "peaking background parameterised by a double Gaussian, and the " \
           "continuum qqbar background by a second order polynomial")
end

alg_m7.note(:background_veto,
            "in the pi+ pi- eta channel the candidate transition photon must be " \
            "separated from clusters formed by charged tracks by more than 17.5 degrees")

### Execution ###
datasets_all = energy_points + [incMC_401, incMC_423, incMC_426,
                                incMC_436, incMC_442, incMC_460] + signal_mc

root_files_m1  = alg_m1.execute_on(datasets_all)
root_files_m2  = alg_m2.execute_on(datasets_all)
root_files_m3  = alg_m3.execute_on(datasets_all)
root_files_m4  = alg_m4.execute_on(datasets_all)
root_files_m5  = alg_m5.execute_on(datasets_all)
root_files_m6  = alg_m6.execute_on(datasets_all)
root_files_m7  = alg_m7.execute_on(datasets_all)
root_files_m8  = alg_m8.execute_on(datasets_all)
root_files_m9  = alg_m9.execute_on(datasets_all)
root_files_m10 = alg_m10.execute_on(datasets_all)
root_files_m11 = alg_m11.execute_on(datasets_all)
root_files_m12 = alg_m12.execute_on(datasets_all)
