# Paper: 2210.13058v2
# Title: e+e- -> gamma phi J/psi via phi chi_c1,c2 at 4.600-4.951 GeV
# Energy: 4.600, 4.612, 4.628, 4.641, 4.661, 4.682, 4.699, 4.740, 4.750, 4.781, 4.843, 4.918, 4.951 GeV (13 points)
# Final state: gamma phi J/psi, with phi -> K+K- or KS0KL0, J/psi -> l+l-
# Three reconstruction modes: 3-track K+K-, 4-track K+K-, KS0KL0

### Dataset preparation ###
# 13 energy points: 703(4600) + 706(4610,4620,4640,4660,4680,4700) + 707(4740,4750,4780,4840,4914,4946)
data_703_4600 = DatasetManager.load_real_data.find("703_4600")
data_706_4610 = DatasetManager.load_real_data.find("706_4610")
data_706_4620 = DatasetManager.load_real_data.find("706_4620")
data_706_4640 = DatasetManager.load_real_data.find("706_4640")
data_706_4660 = DatasetManager.load_real_data.find("706_4660")
data_706_4680 = DatasetManager.load_real_data.find("706_4680")
data_706_4700 = DatasetManager.load_real_data.find("706_4700")
data_707_4740 = DatasetManager.load_real_data.find("707_4740")
data_707_4750 = DatasetManager.load_real_data.find("707_4750")
data_707_4780 = DatasetManager.load_real_data.find("707_4780")
data_707_4840 = DatasetManager.load_real_data.find("707_4840")
data_707_4914 = DatasetManager.load_real_data.find("707_4914")
data_707_4946 = DatasetManager.load_real_data.find("707_4946")

all_data = [
  data_703_4600,
  data_706_4610, data_706_4620, data_706_4640,
  data_706_4660, data_706_4680, data_706_4700,
  data_707_4740, data_707_4750, data_707_4780,
  data_707_4840, data_707_4914, data_707_4946
]

incMC_703_4600 = DatasetManager.load_inclusive_mc.find("703_4600")
incMC_706_4610 = DatasetManager.load_inclusive_mc.find("706_4610")
incMC_706_4620 = DatasetManager.load_inclusive_mc.find("706_4620")
incMC_706_4640 = DatasetManager.load_inclusive_mc.find("706_4640")
incMC_706_4660 = DatasetManager.load_inclusive_mc.find("706_4660")
incMC_706_4680 = DatasetManager.load_inclusive_mc.find("706_4680")
incMC_706_4700 = DatasetManager.load_inclusive_mc.find("706_4700")
incMC_707_4740 = DatasetManager.load_inclusive_mc.find("707_4740")
incMC_707_4750 = DatasetManager.load_inclusive_mc.find("707_4750")
incMC_707_4780 = DatasetManager.load_inclusive_mc.find("707_4780")
incMC_707_4840 = DatasetManager.load_inclusive_mc.find("707_4840")
incMC_707_4914 = DatasetManager.load_inclusive_mc.find("707_4914")
incMC_707_4946 = DatasetManager.load_inclusive_mc.find("707_4946")

all_incMC = [
  incMC_703_4600,
  incMC_706_4610, incMC_706_4620, incMC_706_4640,
  incMC_706_4660, incMC_706_4680, incMC_706_4700,
  incMC_707_4740, incMC_707_4750, incMC_707_4780,
  incMC_707_4840, incMC_707_4914, incMC_707_4946
]

### Decay cards ###
# Card for phi -> K+K- modes
decay_card_KK = <<~DECAYCARD
    Decay psi(4260)
    1.0000  phi  gamma  J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000  e+  e-               PHOTOS VLL;
    Enddecay

    Decay phi
    1.0000  K+  K-               PHSP;
    Enddecay
End
DECAYCARD

# Card for phi -> KS0KL0 mode
decay_card_KSKL = <<~DECAYCARD
    Decay psi(4260)
    1.0000  phi  gamma  J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000  e+  e-               PHOTOS VLL;
    Enddecay

    Decay phi
    1.0000  K_S0  K_L0           PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-             PHSP;
    Enddecay
End
DECAYCARD

exMC_KK = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "phiKK_gammaJpsi_signal"
  config.events         = 100000
  config.decay_card     = decay_card_KK
  config.cross_section  = :default
end

exMC_KSKL = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "phiKSKL_gammaJpsi_signal"
  config.events         = 100000
  config.decay_card     = decay_card_KSKL
  config.cross_section  = :default
end

# Flatten exclusive MC arrays
all_exMC = []
exMC_KK.each { |m| all_exMC << m }
exMC_KSKL.each { |m| all_exMC << m }

### Mode 1: phi -> K+K- with 3 tracks (one kaon missing, 1C fit) ###
alg_KK_3trk = Algorithm.new("GammaPhiJpsi_KK_3trk")
alg_KK_3trk.set_header(["GammaPhiJpsi_KK_3trkAlg/GammaPhiJpsi_KK_3trk.h"])
            .set_alias({"std::vector<double>" => "Vdouble"})

sel_KK_3trk = Selection.new
sel_KK_3trk.select_track {
                  cos_theta 0.93
                  Vz   10.0
                  Vr   1.0
                  nTot "==3"    # 2 leptons + 1 kaon
                }
                .select_photon {
                  tdc_emc_start 0
                  tdc_emc_end   700
                  angle_to_track 10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                }
                .pid(method: :probability) {
                  prob_cut 0.001
                  # High-momentum lepton ID: P > 0.95 GeV/c
                  # Electron if E_EMC > 1.0 GeV, otherwise muon
                  identify_high_momentum_leptons(
                    treat_as_lepton_if_momentum_above: 0.95,
                    treat_as_electron_if_energy_above: 1.0
                  )
                  # Low-momentum kaon ID: L(K) > L(pi), L(K) > 0
                  identify :kaon, against: [:pion]
                  nkp ">=1"
                }
                # 1C kinematic fit: missing kaon mass constrained to nominal K- mass
                .kinematic_fit([:kp, :lp, :lm]) {
                  nominal
                  miss_track_of(:km)
                  chi2_cut 200
                }

alg_KK_3trk.with_decay_card(decay_card_KK).apply(sel_KK_3trk)

### Mode 2: phi -> K+K- with 4 tracks (both kaons detected, 4C fit) ###
alg_KK_4trk = Algorithm.new("GammaPhiJpsi_KK_4trk")
alg_KK_4trk.set_header(["GammaPhiJpsi_KK_4trkAlg/GammaPhiJpsi_KK_4trk.h"])
            .set_alias({"std::vector<double>" => "Vdouble"})

sel_KK_4trk = Selection.new
sel_KK_4trk.select_track {
                  cos_theta 0.93
                  Vz   10.0
                  Vr   1.0
                  nTot "==4"    # 2 leptons + 2 kaons
                }
                .select_photon {
                  tdc_emc_start 0
                  tdc_emc_end   700
                  angle_to_track 10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                }
                .pid(method: :probability) {
                  prob_cut 0.001
                  identify_high_momentum_leptons(
                    treat_as_lepton_if_momentum_above: 0.95,
                    treat_as_electron_if_energy_above: 1.0
                  )
                  identify :kaon, against: [:pion]
                  nkp "==1"
                  nkm "==1"
                }
                # 4C kinematic fit: total four-momentum constrained to CMS
                .kinematic_fit([:kp, :km, :lp, :lm]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }

alg_KK_4trk.with_decay_card(decay_card_KK).apply(sel_KK_4trk)

### Mode 3: phi -> KS0KL0 (1C fit with missing KL0) ###
alg_KSKL = Algorithm.new("GammaPhiJpsi_KSKL")
alg_KSKL.set_header(["GammaPhiJpsi_KSKLAlg/GammaPhiJpsi_KSKL.h"])
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_KSKL = Selection.new
sel_KSKL.select_track {
                  cos_theta 0.93
                  Vz   10.0
                  Vr   1.0
                  nTot ">=4"    # 2 leptons + 2 pions from KS0
                }
                .select_photon {
                  tdc_emc_start 0
                  tdc_emc_end   700
                  angle_to_track 10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                }
                .pid(method: :probability) {
                  prob_cut 0.001
                  identify_high_momentum_leptons(
                    treat_as_lepton_if_momentum_above: 0.95,
                    treat_as_electron_if_energy_above: 1.0
                  )
                }
                .secondary_vertex_fit([:pip, :pim]) {
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
                # 1C kinematic fit: missing KL0 mass constrained
                .kinematic_fit([:K_S0, :lp, :lm]) {
                  nominal
                  miss_track_of(:K_L0)
                  chi2_cut 200
                }

alg_KSKL.with_decay_card(decay_card_KSKL).apply(sel_KSKL)

### Common notes ###
# J/psi mass window applied in ROOT
alg_KK_3trk.note(:jpsi_mass_window,
  "M(l+l-) within (m_J/psi - 3*sigma, m_J/psi + 3*sigma) for signal region; sidebands for background estimation; applied in ROOT")
alg_KK_4trk.note(:jpsi_mass_window,
  "M(l+l-) within (m_J/psi - 3*sigma, m_J/psi + 3*sigma) for signal region; applied in ROOT")
alg_KSKL.note(:jpsi_mass_window,
  "M(l+l-) within (m_J/psi - 3*sigma, m_J/psi + 3*sigma) for signal region; applied in ROOT")

# KS0 mass window (for mode 3)
alg_KSKL.note(:ks0_mass_window,
  "M(pi+pi-) within 3*sigma of KS0 mass; decay length > 2*sigma; applied at vertex fit level")

# Chi2 selection for 3-track: the 1C fit with missing kaon selects the
# combination with smallest chi2_1C
alg_KK_3trk.note(:chi2_selection,
  "Combination with smallest chi2_1C selected; additional muon hit depth > 40cm requirement for muon ID; applied in ROOT")

# The three alg outputs are merged in the ROOT analysis to extract
# chi_c1 and chi_c2 signals via fit to M(phi J/psi)
alg_KK_3trk.note(:merge_with_other_modes,
  "Results merged with KK_4trk and KSKL modes. chi_c1/chi_c2 yields extracted from M(phi J/psi) fit in ROOT analysis.")
alg_KK_4trk.note(:merge_with_other_modes,
  "Results merged with KK_3trk and KSKL modes.")
alg_KSKL.note(:merge_with_other_modes,
  "Results merged with KK_3trk and KK_4trk modes.")

# chi_c1 vs chi_c2: in the signal, both chi_c1 and chi_c2 contribute
# and are separated in the M(phi J/psi) fit in ROOT
alg_KK_3trk.note(:chic1_chic2_separation,
  "chi_c1 and chi_c2 separated by fitting M(phi J/psi) spectrum; applied in ROOT")
alg_KK_4trk.note(:chic1_chic2_separation,
  "chi_c1 and chi_c2 separated by fitting M(phi J/psi) spectrum; applied in ROOT")
alg_KSKL.note(:chic1_chic2_separation,
  "chi_c1 and chi_c2 separated by fitting M(phi J/psi) spectrum; applied in ROOT")

# Execute
all_datasets = all_data + all_incMC + all_exMC
root_files_3trk = alg_KK_3trk.execute_on(all_datasets)
root_files_4trk = alg_KK_4trk.execute_on(all_datasets)
root_files_KSKL = alg_KSKL.execute_on(all_datasets)