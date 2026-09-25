# Paper 2310.07277v1: Search for J/psi weak decays to D meson + light hadron
# BESIII, (10087±44)×10^6 J/psi events
# Five signal modes with semileptonic D-meson tagging (missing neutrino)

jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ==============================
# Mode 1: J/psi -> anti-D0 pi0 (with anti-D0 -> K+ e- nu_e_bar)
# ==============================
decay_card_d0pi0 = <<~DECAYCARD
    Decay J/psi
    1.000 anti-D0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ e- anti-nu_e PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_d0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_antiD0pi0_KeNu_mc"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_d0pi0
  config.cross_section = :default
end

alg_d0pi0 = Algorithm.new("Jpsi2AntiD0Pi0")
alg_d0pi0.set_header(["Jpsi2AntiD0Pi0Alg/Jpsi2AntiD0Pi0.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })

sel_d0pi0 = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot "==2"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.8
    identify :kaon, against: [:pion]
    nkp "==1"
    nkm "==0"
    nlp "==0"
    nlm "==1"
    nep "==0"
    nem "==1"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=1"
  end
  .kinematic_fit([:kp, :em, :pi0]) do
    nominal
    miss_track_of(:nu_e)
    chi2_cut 200
  end

alg_d0pi0.note(:Um_iss_cut, "Um_iss = E_miss - c|p_miss| in (-0.083, 0.119) GeV; identifies the missing neutrino")
  .note(:D0_recoil_mass, "Recoil mass of pi0 in (1.80, 1.95) GeV/c^2 (anti-D0 window)")
  .note(:electron_pid, "CL_e > CL_pi, CL_e > CL_K, CL_e > 0.001; E/p > 0.8")
  .note(:charge_conjugate, "J/psi -> D0 pi0 with D0 -> K- e+ nu_e also reconstructed")
  .with_decay_card(decay_card_d0pi0)
  .apply(sel_d0pi0)

alg_d0pi0.execute_on([jpsi_data, jpsi_incMC, exMC_d0pi0])

# ==============================
# Mode 2: J/psi -> anti-D0 eta (with anti-D0 -> K+ e- nu_e_bar)
# ==============================
decay_card_d0eta = <<~DECAYCARD
    Decay J/psi
    1.000 anti-D0 eta PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ e- anti-nu_e PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_d0eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_antiD0eta_KeNu_mc"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_d0eta
  config.cross_section = :default
end

alg_d0eta = Algorithm.new("Jpsi2AntiD0Eta")
alg_d0eta.set_header(["Jpsi2AntiD0EtaAlg/Jpsi2AntiD0Eta.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })

sel_d0eta = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot "==2"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.8
    identify :kaon, against: [:pion]
    nkp "==1"
    nkm "==0"
    nlp "==0"
    nlm "==1"
    nep "==0"
    nem "==1"
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 20
    neta ">=1"
  end
  .kinematic_fit([:kp, :em, :eta]) do
    nominal
    miss_track_of(:nu_e)
    chi2_cut 200
  end

alg_d0eta.note(:Um_iss_cut, "Um_iss in (-0.050, 0.060) GeV")
  .note(:D0_recoil_mass, "Recoil mass of eta in (1.80, 1.95) GeV/c^2 (anti-D0 window)")
  .note(:electron_pid, "CL_e > CL_pi, CL_e > CL_K, CL_e > 0.001; E/p > 0.8")
  .note(:eta_mass_window, "0.50 < M(gamma gamma) < 0.57 GeV/c^2")
  .with_decay_card(decay_card_d0eta)
  .apply(sel_d0eta)

alg_d0eta.execute_on([jpsi_data, jpsi_incMC, exMC_d0eta])

# ==============================
# Mode 3: J/psi -> anti-D0 rho0 (anti-D0 -> K+ e- nu_e_bar; rho0 -> pi+ pi-)
# ==============================
decay_card_d0rho0 = <<~DECAYCARD
    Decay J/psi
    1.000 anti-D0 rho0 PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ e- anti-nu_e PHSP;
    Enddecay

    Decay rho0
    1.000 pi+ pi- VSS;
    Enddecay

    End
DECAYCARD

exMC_d0rho0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_antiD0rho0_KeNu_mc"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_d0rho0
  config.cross_section = :default
end

alg_d0rho0 = Algorithm.new("Jpsi2AntiD0Rho0")
alg_d0rho0.set_header(["Jpsi2AntiD0Rho0Alg/Jpsi2AntiD0Rho0.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })

sel_d0rho0 = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot "==4"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.8
    identify :kaon, against: [:pion]
    nkp "==1"
    nkm "==0"
    nlp "==0"
    nlm "==1"
    nep "==0"
    nem "==1"
  end
  .remove([:kp <= :chrgp])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==1"
    npim "==1"
  end
  .kinematic_fit([:kp, :em, :pip, :pim]) do
    nominal
    miss_track_of(:nu_e)
    chi2_cut 200
  end

alg_d0rho0.note(:Um_iss_cut, "Um_iss in (-0.040, 0.050) GeV")
  .note(:D0_recoil_mass, "Recoil mass of rho0 (pi+pi-) in (1.80, 1.95) GeV/c^2")
  .note(:rho0_mass_window, "0.62 < M(pi+pi-) < 0.95 GeV/c^2")
  .with_decay_card(decay_card_d0rho0)
  .apply(sel_d0rho0)

alg_d0rho0.execute_on([jpsi_data, jpsi_incMC, exMC_d0rho0])

# ==============================
# Mode 4: J/psi -> D- pi+ (D- -> K_S0 e- nu_e_bar; K_S0 -> pi+ pi-)
# ==============================
decay_card_dmpi = <<~DECAYCARD
    Decay J/psi
    1.000 D- pi+ PHSP;
    Enddecay

    Decay D-
    1.000 K_S0 e- anti-nu_e PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_dmpi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_Dmpi_KS0eNu_mc"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_dmpi
  config.cross_section = :default
end

alg_dmpi = Algorithm.new("Jpsi2DmPiP")
alg_dmpi.set_header(["Jpsi2DmPiPAlg/Jpsi2DmPiP.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })

sel_dmpi = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot "==4"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.8
    nlp "==0"
    nlm "==1"
    nep "==0"
    nem "==1"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=3"
    npim ">=1"
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:K_S0, :em, :pip]) do
    nominal
    miss_track_of(:nu_e)
    chi2_cut 200
  end

alg_dmpi.note(:Um_iss_cut, "Um_iss in (-0.037, 0.040) GeV")
  .note(:Dm_recoil_mass, "Recoil mass of pi+ in (1.80, 1.95) GeV/c^2 (D- window)")
  .note(:KS0_reconstruction, "|M(pi+pi-) - m_K_S0| < 12 MeV/c^2; decay length > 2*resolution; secondary vertex chi2 chosen")
  .with_decay_card(decay_card_dmpi)
  .apply(sel_dmpi)

alg_dmpi.execute_on([jpsi_data, jpsi_incMC, exMC_dmpi])

# ==============================
# Mode 5: J/psi -> D- rho+ (D- -> K_S0 e- nu_e_bar; rho+ -> pi+ pi0)
# ==============================
decay_card_dmrho = <<~DECAYCARD
    Decay J/psi
    1.000 D- rho+ PHSP;
    Enddecay

    Decay D-
    1.000 K_S0 e- anti-nu_e PHSP;
    Enddecay

    Decay rho+
    1.000 pi+ pi0 VSS;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_dmrho = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_DmRhop_KS0eNu_mc"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_dmrho
  config.cross_section = :default
end

alg_dmrho = Algorithm.new("Jpsi2DmRhoP")
alg_dmrho.set_header(["Jpsi2DmRhoPAlg/Jpsi2DmRhoP.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })

sel_dmrho = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot "==4"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.8
    nlp "==0"
    nlm "==1"
    nep "==0"
    nem "==1"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=3"
    npim ">=1"
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
  .kinematic_fit([:K_S0, :em, :pip, :pi0]) do
    nominal
    miss_track_of(:nu_e)
    chi2_cut 200
  end

alg_dmrho.note(:Um_iss_cut, "Um_iss in (-0.058, 0.074) GeV")
  .note(:Dm_recoil_mass, "Recoil mass of rho+ (pi+ pi0) in (1.80, 1.95) GeV/c^2")
  .note(:rho_plus_mass_window, "0.62 < M(pi+ pi0) < 0.95 GeV/c^2")
  .with_decay_card(decay_card_dmrho)
  .apply(sel_dmrho)

alg_dmrho.execute_on([jpsi_data, jpsi_incMC, exMC_dmrho])