# Paper: 2310.07277v1 - J/psi -> D meson weak decays (5 modes)
# J/psi -> D0 pi0 / D0 eta / D0 rho0 / D- pi+ / D- rho+
# D meson tagged via semileptonic: D0->K+ e- nu_e_bar, D-->K_S0 e- nu_e_bar
# Missing neutrino handled via U_miss = E_miss - c|p_miss|

DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# J/psi -> D0 pi0 + c.c., with D0 -> K+ e- nu_e_bar, pi0 -> gamma gamma
decay_card_d0pi0 = <<~DECAYCARD
  Decay J/psi
  1.0000 anti-D0 pi0 PHSP;
  Enddecay
  Decay anti-D0
  1.0000 K+ e- anti-nu_e PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# J/psi -> D0 eta + c.c., D0 -> K+ e- nu_e_bar, eta -> gamma gamma
decay_card_d0eta = <<~DECAYCARD
  Decay J/psi
  1.0000 anti-D0 eta PHSP;
  Enddecay
  Decay anti-D0
  1.0000 K+ e- anti-nu_e PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# J/psi -> D0 rho0 + c.c., D0 -> K+ e- nu_e_bar, rho0 -> pi+ pi-
decay_card_d0rho0 = <<~DECAYCARD
  Decay J/psi
  1.0000 anti-D0 rho0 PHSP;
  Enddecay
  Decay anti-D0
  1.0000 K+ e- anti-nu_e PHSP;
  Enddecay
  Decay rho0
  1.0000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# J/psi -> D- pi+ + c.c., D- -> K_S0 e- nu_e_bar, K_S0 -> pi+ pi-
decay_card_dmpi = <<~DECAYCARD
  Decay J/psi
  1.0000 D+ pi- PHSP;
  Enddecay
  Decay D+
  1.0000 anti-K_S0 e+ nu_e PHSP;
  Enddecay
  Decay anti-K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# J/psi -> D- rho+ + c.c., D- -> K_S0 e- nu_e_bar, rho+ -> pi+ pi0
decay_card_dmrho = <<~DECAYCARD
  Decay J/psi
  1.0000 D+ rho- PHSP;
  Enddecay
  Decay D+
  1.0000 anti-K_S0 e+ nu_e PHSP;
  Enddecay
  Decay anti-K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay
  Decay rho-
  1.0000 pi- pi0 PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Signal MC for all 5 modes
signal_mc_d0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_d0pi0"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_d0pi0
  config.cross_section   = :default
end

signal_mc_d0eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_d0eta"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_d0eta
  config.cross_section   = :default
end

signal_mc_d0rho0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_d0rho0"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_d0rho0
  config.cross_section   = :default
end

signal_mc_dmpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_dmpi"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_dmpi
  config.cross_section   = :default
end

signal_mc_dmrho = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_dmrho"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_dmrho
  config.cross_section   = :default
end

# ============================================================================
# Mode 1: J/psi -> D0 pi0 + c.c., D0 -> K+ e- nu_e_bar
# ============================================================================
algorithm_d0pi0 = Algorithm.new("D0Pi0", "00-00-01")

selection_d0pi0 = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nTot        "==2"      # K+ and e- (2 tracks for D0 semileptonic tag)
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion]
    identify :electron, against: [:pion, :kaon]
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=1"
  end
  .partial_miss([:nu_e_bar]) do
    best_combination_by_mass
  end

algorithm_d0pi0
  .set_header(["DWeakDecayAlg/DWeakDecayAlg.h"])
  .set_constant(ECMS: 3.097)
  .note(:mode, "J/psi -> D0 pi0 + c.c. D0 -> K+ e- anti-nu_e. pi0 -> gamma gamma with 1C kinematic fit chi2 < 20. M(gamma gamma) in [0.115, 0.150] GeV/c2.")
  .note(:electron_pid, "Electron: CL_e > CL_pi, CL_e > CL_K, CL_e > 0.001. Additionally E_e/p_e > 0.8 required.")
  .note(:kaon_pid, "Kaon: CL_K > CL_pi.")
  .note(:U_miss, "U_miss = E_miss - c|p_miss| used to identify missing neutrino. Signal region: (-0.083, 0.119) GeV.")
  .note(:recoil_mass, "Recoil mass against pi0 required in (1.80, 1.95) GeV/c2 (D meson mass window).")
  .note(:upper_limit, "No significant signal observed. B(J/psi -> D0 pi0 + c.c.) < 4.7 x 10^-7 at 90% CL.")
  .with_decay_card(decay_card_d0pi0)
  .apply(selection_d0pi0)

algorithm_d0pi0.execute_on([jpsi_data, jpsi_incMC, signal_mc_d0pi0])

# ============================================================================
# Mode 2: J/psi -> D0 eta + c.c., D0 -> K+ e- nu_e_bar
# ============================================================================
algorithm_d0eta = Algorithm.new("D0Eta", "00-00-01")

selection_d0eta = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nTot        "==2"      # K+ and e-
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion]
    identify :electron, against: [:pion, :kaon]
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 20
    neta ">=1"
  end
  .partial_miss([:nu_e_bar]) do
    best_combination_by_mass
  end

algorithm_d0eta
  .set_header(["DWeakDecayAlg/DWeakDecayAlg.h"])
  .set_constant(ECMS: 3.097)
  .note(:mode, "J/psi -> D0 eta + c.c. D0 -> K+ e- anti-nu_e. eta -> gamma gamma with 1C kinematic fit chi2 < 20. M(gamma gamma) in [0.50, 0.57] GeV/c2.")
  .note(:U_miss, "U_miss signal region: (-0.050, 0.060) GeV.")
  .note(:recoil_mass, "Recoil mass against eta in (1.80, 1.95) GeV/c2.")
  .note(:upper_limit, "B(J/psi -> D0 eta + c.c.) < 6.8 x 10^-7 at 90% CL.")
  .with_decay_card(decay_card_d0eta)
  .apply(selection_d0eta)

algorithm_d0eta.execute_on([jpsi_data, jpsi_incMC, signal_mc_d0eta])

# ============================================================================
# Mode 3: J/psi -> D0 rho0 + c.c., D0 -> K+ e- nu_e_bar, rho0 -> pi+ pi-
# ============================================================================
algorithm_d0rho0 = Algorithm.new("D0Rho0", "00-00-01")

selection_d0rho0 = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nTot        "==4"      # K+, e-, pi+, pi-
    nNet        "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion]
    identify :electron, against: [:pion, :kaon]
    identify :pion, against: [:kaon, :proton]
  end
  .partial_miss([:nu_e_bar]) do
    best_combination_by_mass
  end

algorithm_d0rho0
  .set_header(["DWeakDecayAlg/DWeakDecayAlg.h"])
  .set_constant(ECMS: 3.097)
  .note(:mode, "J/psi -> D0 rho0 + c.c. D0 -> K+ e- anti-nu_e. rho0 -> pi+ pi-. rho mass window: M(pi+ pi-) in [0.62, 0.95] GeV/c2.")
  .note(:U_miss, "U_miss signal region: (-0.040, 0.050) GeV.")
  .note(:recoil_mass, "Recoil mass against pi+ pi- system in (1.80, 1.95) GeV/c2.")
  .note(:upper_limit, "B(J/psi -> D0 rho0 + c.c.) < 5.2 x 10^-7 at 90% CL.")
  .with_decay_card(decay_card_d0rho0)
  .apply(selection_d0rho0)

algorithm_d0rho0.execute_on([jpsi_data, jpsi_incMC, signal_mc_d0rho0])

# ============================================================================
# Mode 4: J/psi -> D- pi+ + c.c., D- -> K_S0 e- nu_e_bar, K_S0 -> pi+ pi-
# ============================================================================
algorithm_dmpi = Algorithm.new("DmPi", "00-00-01")

selection_dmpi = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nTot        "==4"      # pi+ (from D- side), e+, pi+ pi- (from K_S0)
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :electron, against: [:pion, :kaon]
    identify :pion, against: [:kaon, :proton]
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .partial_miss([:nu_e_bar]) do
    best_combination_by_mass
  end

algorithm_dmpi
  .set_header(["DWeakDecayAlg/DWeakDecayAlg.h"])
  .set_constant(ECMS: 3.097)
  .note(:mode, "J/psi -> D- pi+ + c.c. D- -> K_S0 e- anti-nu_e. K_S0 -> pi+ pi- with secondary vertex fit, |M(pi+ pi-) - m_K_S0| < 12 MeV/c2. Decay length > 2 sigma from IP. Best K_S0 by min chi2.")
  .note(:K_S0_vertex, "K_S0: |Vz| < 20 cm for daughter tracks. Secondary vertex fit. Invariant mass within 12 MeV/c2 of nominal K_S0 mass. Decay length > 2x vertex resolution.")
  .note(:U_miss, "U_miss signal region: (-0.037, 0.040) GeV.")
  .note(:recoil_mass, "Recoil mass against pi+ in (1.80, 1.95) GeV/c2.")
  .note(:upper_limit, "B(J/psi -> D- pi+ + c.c.) < 7.0 x 10^-8 at 90% CL.")
  .with_decay_card(decay_card_dmpi)
  .apply(selection_dmpi)

algorithm_dmpi.execute_on([jpsi_data, jpsi_incMC, signal_mc_dmpi])

# ============================================================================
# Mode 5: J/psi -> D- rho+ + c.c., D- -> K_S0 e- nu_e_bar, rho+ -> pi+ pi0
# ============================================================================
algorithm_dmrho = Algorithm.new("DmRho", "00-00-01")

selection_dmrho = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nTot        "==4"      # pi+ (from rho+ side), e+, pi+ pi- (from K_S0)
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :electron, against: [:pion, :kaon]
    identify :pion, against: [:kaon, :proton]
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=1"
  end
  .partial_miss([:nu_e_bar]) do
    best_combination_by_mass
  end

algorithm_dmrho
  .set_header(["DWeakDecayAlg/DWeakDecayAlg.h"])
  .set_constant(ECMS: 3.097)
  .note(:mode, "J/psi -> D- rho+ + c.c. D- -> K_S0 e- anti-nu_e. rho+ -> pi+ pi0. K_S0 -> pi+ pi- with secondary vertex fit. pi0 from 1C kinematic fit. rho mass window: M(pi+ pi0) in [0.62, 0.95] GeV/c2.")
  .note(:U_miss, "U_miss signal region: (-0.058, 0.074) GeV.")
  .note(:recoil_mass, "Recoil mass against pi+ pi0 system in (1.80, 1.95) GeV/c2.")
  .note(:upper_limit, "B(J/psi -> D- rho+ + c.c.) < 6.0 x 10^-7 at 90% CL.")
  .with_decay_card(decay_card_dmrho)
  .apply(selection_dmrho)

algorithm_dmrho.execute_on([jpsi_data, jpsi_incMC, signal_mc_dmrho])