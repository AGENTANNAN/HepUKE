# Paper: 2211.08561v5
# Title: Observation of Y(4230) and Y(4710) in e+e- -> KS0 KS0 J/psi
# Energy: 4.128 to 4.950 GeV (36 energy points)
# Final state: KS0 KS0 J/psi, J/psi -> l+l- (e+e- or mu+mu-)
# Two methods: two-KS0 (full, 4C fit) & one-KS0 (partial, 1C fit)

### Dataset preparation ###
# 36 energy points across BOSS 703, 705, 706, 707
data_703_4180 = DatasetManager.load_real_data.find("703_4180")
data_703_4190 = DatasetManager.load_real_data.find("703_4190")
data_703_4200 = DatasetManager.load_real_data.find("703_4200")
data_703_4210 = DatasetManager.load_real_data.find("703_4210")
data_703_4220 = DatasetManager.load_real_data.find("703_4220")
data_703_4230 = DatasetManager.load_real_data.find("703_4230")
data_703_4237 = DatasetManager.load_real_data.find("703_4237")
data_703_4246 = DatasetManager.load_real_data.find("703_4246")
data_703_4260 = DatasetManager.load_real_data.find("703_4260")
data_703_4270 = DatasetManager.load_real_data.find("703_4270")
data_703_4280 = DatasetManager.load_real_data.find("703_4280")
data_703_4360 = DatasetManager.load_real_data.find("703_4360")
data_703_4420 = DatasetManager.load_real_data.find("703_4420")
data_703_4470 = DatasetManager.load_real_data.find("703_4470")
data_703_4530 = DatasetManager.load_real_data.find("703_4530")
data_703_4575 = DatasetManager.load_real_data.find("703_4575")
data_703_4600 = DatasetManager.load_real_data.find("703_4600")
data_705_4130 = DatasetManager.load_real_data.find("705_4130")
data_705_4160 = DatasetManager.load_real_data.find("705_4160")
data_705_4290 = DatasetManager.load_real_data.find("705_4290")
data_705_4315 = DatasetManager.load_real_data.find("705_4315")
data_705_4340 = DatasetManager.load_real_data.find("705_4340")
data_705_4380 = DatasetManager.load_real_data.find("705_4380")
data_705_4400 = DatasetManager.load_real_data.find("705_4400")
data_705_4440 = DatasetManager.load_real_data.find("705_4440")
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
  data_703_4180, data_703_4190, data_703_4200, data_703_4210,
  data_703_4220, data_703_4230, data_703_4237, data_703_4246,
  data_703_4260, data_703_4270, data_703_4280, data_703_4360,
  data_703_4420, data_703_4470, data_703_4530, data_703_4575,
  data_703_4600,
  data_705_4130, data_705_4160, data_705_4290, data_705_4315,
  data_705_4340, data_705_4380, data_705_4400, data_705_4440,
  data_706_4610, data_706_4620, data_706_4640, data_706_4660,
  data_706_4680, data_706_4700,
  data_707_4740, data_707_4750, data_707_4780, data_707_4840,
  data_707_4914, data_707_4946
]

incMC_703_4180 = DatasetManager.load_inclusive_mc.find("703_4180")
incMC_703_4190 = DatasetManager.load_inclusive_mc.find("703_4190")
incMC_703_4200 = DatasetManager.load_inclusive_mc.find("703_4200")
incMC_703_4210 = DatasetManager.load_inclusive_mc.find("703_4210")
incMC_703_4220 = DatasetManager.load_inclusive_mc.find("703_4220")
incMC_703_4230 = DatasetManager.load_inclusive_mc.find("703_4230")
incMC_703_4237 = DatasetManager.load_inclusive_mc.find("703_4237")
incMC_703_4246 = DatasetManager.load_inclusive_mc.find("703_4246")
incMC_703_4260 = DatasetManager.load_inclusive_mc.find("703_4260")
incMC_703_4270 = DatasetManager.load_inclusive_mc.find("703_4270")
incMC_703_4280 = DatasetManager.load_inclusive_mc.find("703_4280")
incMC_703_4360 = DatasetManager.load_inclusive_mc.find("703_4360")
incMC_703_4420 = DatasetManager.load_inclusive_mc.find("703_4420")
incMC_703_4600 = DatasetManager.load_inclusive_mc.find("703_4600")
incMC_705_4130 = DatasetManager.load_inclusive_mc.find("705_4130")
incMC_705_4160 = DatasetManager.load_inclusive_mc.find("705_4160")
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
  incMC_703_4180, incMC_703_4190, incMC_703_4200, incMC_703_4210,
  incMC_703_4220, incMC_703_4230, incMC_703_4237, incMC_703_4246,
  incMC_703_4260, incMC_703_4270, incMC_703_4280, incMC_703_4360,
  incMC_703_4420, incMC_703_4600,
  incMC_705_4130, incMC_705_4160,
  incMC_706_4610, incMC_706_4620, incMC_706_4640, incMC_706_4660,
  incMC_706_4680, incMC_706_4700,
  incMC_707_4740, incMC_707_4750, incMC_707_4780, incMC_707_4840,
  incMC_707_4914, incMC_707_4946
]

# Decay card: e+e- -> KS0 KS0 J/psi (continuum + resonance production)
# Uses KKMC generator with psi(4260) as top mother
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.0000  K_S0  K_S0  J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000  e+  e-               PHOTOS VLL;
    Enddecay

    Decay K_S0
    1.0000  pi+  pi-             PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "KSKSJpsi_signal"
  config.events         = 100000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
# Two-KS reconstruction method (full 6-track reconstruction, 4C kinematic fit)
alg_2KS = Algorithm.new("KSKSJpsi_2KS")
alg_2KS.set_header(["KSKSJpsi_2KSAlg/KSKSJpsi_2KS.h"])
        .set_alias({"std::vector<double>" => "Vdouble"})

selection_2KS = Selection.new
selection_2KS.select_track {
                  cos_theta 0.93
                  Vz   10.0
                  Vr   1.0
                  nTot ">=6"
                }
                .select_photon {
                  tdc_emc_start 0
                  tdc_emc_end   700
                  angle_to_track 10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                }
                # High-momentum lepton ID:
                # P > 0.95 GeV/c; electron if E_EMC > 0.95 GeV, else muon
                .pid(method: :probability) {
                  prob_cut 0.001
                  identify_high_momentum_leptons(
                    treat_as_lepton_if_momentum_above: 0.95,
                    treat_as_electron_if_energy_above: 0.95
                  )
                }
                # KS0 from secondary vertex fit of oppositely charged pions
                # Pions: P < 0.95 GeV/c (non-lepton tracks), |Vz| < 20 cm
                .secondary_vertex_fit([:pip, :pim]) {
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
                .secondary_vertex_fit([:pip, :pim]) {
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
                # 4C kinematic fit: four-momentum constrained to CM
                .kinematic_fit([:K_S0, :K_S0, :lp, :lm]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }

alg_2KS.with_decay_card(decay_card).apply(selection_2KS)

# One-KS reconstruction method (4-5 track, 1C fit with missing KS0)
# Missing KS0 decays to pi0pi0 (unobserved)
alg_1KS = Algorithm.new("KSKSJpsi_1KS")
alg_1KS.set_header(["KSKSJpsi_1KSAlg/KSKSJpsi_1KS.h"])
        .set_alias({"std::vector<double>" => "Vdouble"})

selection_1KS = Selection.new
selection_1KS.select_track {
                  cos_theta 0.93
                  Vz   10.0
                  Vr   1.0
                  nTot ">=4"
                }
                .select_photon {
                  tdc_emc_start 0
                  tdc_emc_end   700
                  angle_to_track 10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam ">=2"
                }
                .pid(method: :probability) {
                  prob_cut 0.001
                  identify_high_momentum_leptons(
                    treat_as_lepton_if_momentum_above: 0.95,
                    treat_as_electron_if_energy_above: 0.95
                  )
                }
                .secondary_vertex_fit([:pip, :pim]) {
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
                # 1C kinematic fit: missing mass constrained to KS0 mass
                .kinematic_fit([:K_S0, :lp, :lm]) {
                  nominal
                  miss_track_of(:K_S0)
                  chi2_cut 20
                }

alg_1KS.with_decay_card(decay_card).apply(selection_1KS)

# J/psi mass window: (m_J/psi - 3*sigma, m_J/psi + 3*sigma) applied in ROOT
alg_2KS.note(:jpsi_mass_window,
  "M(l+l-) in (m_J/psi - 3*sigma, m_J/psi + 3*sigma) for signal region; sidebands (m-13*sigma,m-7*sigma) and (m+7*sigma,m+13*sigma); applied in ROOT")
alg_1KS.note(:jpsi_mass_window,
  "M(l+l-) in (m_J/psi - 3*sigma, m_J/psi + 3*sigma) for signal region; applied in ROOT")

# KS0 mass window: (m_KS - 3*sigma, m_KS + 3*sigma) applied in ROOT
alg_2KS.note(:ks0_mass_window,
  "M(pi+pi-) in (m_KS0 - 3*sigma_KS0, m_KS0 + 3*sigma_KS0); applied in ROOT")
alg_1KS.note(:ks0_mass_window,
  "M(pi+pi-) in (m_KS0 - 3*sigma_KS0, m_KS0 + 3*sigma_KS0); applied in ROOT")

# KS0 decay length > 2*sigma requirement
alg_2KS.note(:ks0_decay_length, "Decay length > 2*vertex_resolution; applied in ROOT")
alg_1KS.note(:ks0_decay_length, "Decay length > 2*vertex_resolution; applied in ROOT")

# Note about the FOM optimization
alg_2KS.note(:chi2_optimization,
  "chi2 cuts (200 for 4C, 20 for 1C) optimized via FOM = S/sqrt(S+B); applied at execution")

all_datasets = all_data + all_incMC + exMC_signal
root_files_2KS = alg_2KS.execute_on(all_datasets)
root_files_1KS = alg_1KS.execute_on(all_datasets)