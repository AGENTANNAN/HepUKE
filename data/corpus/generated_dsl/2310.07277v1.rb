# =====================================================================
# BESIII: search for J/psi weak decays J/psi -> D + light hadron
#   Modes : anti-D0 pi0, anti-D0 eta, anti-D0 rho0, D- pi+, D- rho+
#   anti-D0 -> K+ e- anti-nu_e ;  D- -> K_S0 e- anti-nu_e
#   (charge conjugates included in the analysis)
# =====================================================================

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

### Decay cards (EvtGen format) ###
# Mode I : J/psi -> anti-D0 pi0
decay_card_d0bar_pi0 = <<~DECAYCARD
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

# Mode II : J/psi -> anti-D0 eta
decay_card_d0bar_eta = <<~DECAYCARD
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

# Mode III : J/psi -> anti-D0 rho0
decay_card_d0bar_rho0 = <<~DECAYCARD
    Decay J/psi
    1.0000 anti-D0 rho0 PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ e- anti-nu_e PHSP;
    Enddecay

    Decay rho0
    1.0000 pi+ pi- VSS;
    Enddecay

    End
DECAYCARD

# Mode IV : J/psi -> D- pi+
decay_card_dm_pip = <<~DECAYCARD
    Decay J/psi
    1.0000 D- pi+ PHSP;
    Enddecay

    Decay D-
    1.0000 K_S0 e- anti-nu_e PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode V : J/psi -> D- rho+
decay_card_dm_rhop = <<~DECAYCARD
    Decay J/psi
    1.0000 D- rho+ PHSP;
    Enddecay

    Decay D-
    1.0000 K_S0 e- anti-nu_e PHSP;
    Enddecay

    Decay rho+
    1.0000 pi+ pi0 VSS;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples (100k events each) ###
exMC_d0bar_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_weak_d0bar_pi0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_d0bar_pi0
  config.cross_section   = :default
end

exMC_d0bar_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_weak_d0bar_eta"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_d0bar_eta
  config.cross_section   = :default
end

exMC_d0bar_rho0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_weak_d0bar_rho0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_d0bar_rho0
  config.cross_section   = :default
end

exMC_dm_pip = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_weak_dm_pip"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_dm_pip
  config.cross_section   = :default
end

exMC_dm_rhop = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_weak_dm_rhop"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_dm_rhop
  config.cross_section   = :default
end

# =====================================================================
### Mode I : J/psi -> anti-D0 pi0, anti-D0 -> K+ e- anti-nu_e ###
# =====================================================================
alg_d0bar_pi0 = Algorithm.new("JpsiD0barPi0")
alg_d0bar_pi0.set_header(["JpsiD0barPi0Alg/JpsiD0barPi0.h"])
              .set_constant({"ECMS" => [:double, 3.097]})
              .set_alias({"std::vector<double>" => "Vdouble"})

sel_d0bar_pi0 = Selection.new
sel_d0bar_pi0
  .select_track {                 # common track selection + 2-track topology
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        10.0                # |Vz| < 10 cm
    Vr        1.0                 # Vr < 1 cm
    nChrp     "==1"               # one positive track (K+)
    nChrn     "==1"               # one negative track (e-)
    nNet      "==0"               # net charge zero
  }
  .select_photon {                # common photon selection; >=2 photons for pi0 -> gamma gamma
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025       # >25 MeV in barrel
    energyThreshold_e 0.050       # >50 MeV in endcap
    nGam              ">=2"
  }
  .pid(method: :probability) {    # probability PID, CL > 0.001
    prob_cut 0.001
    # tracks with p > 1.0 treated as leptons; electron if EMC energy > 0.8 GeV, else muon
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.8
    identify :kaon, against: [:pion]
    nkp "==1"                     # one K+
    nlm "==1"                     # one e-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # gamma gamma mass-constrained fit to pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=1"
  }
  .kinematic_fit([:kp, :lm, :pi0]) {          # nominal fit with a missing neutrino
    nominal
    constrain_four_momentum
    miss_track_of :nu_e
    chi2_cut 200
  }

alg_d0bar_pi0.with_decay_card(decay_card_d0bar_pi0).apply(sel_d0bar_pi0)

# =====================================================================
### Mode II : J/psi -> anti-D0 eta, anti-D0 -> K+ e- anti-nu_e ###
# =====================================================================
alg_d0bar_eta = Algorithm.new("JpsiD0barEta")
alg_d0bar_eta.set_header(["JpsiD0barEtaAlg/JpsiD0barEta.h"])
              .set_constant({"ECMS" => [:double, 3.097]})
              .set_alias({"std::vector<double>" => "Vdouble"})

sel_d0bar_eta = Selection.new
sel_d0bar_eta
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
    nGam              ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.8
    identify :kaon, against: [:pion]
    nkp "==1"
    nlm "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # gamma gamma mass-constrained fit to eta
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 20
    neta ">=1"
  }
  .kinematic_fit([:kp, :lm, :eta]) {          # nominal fit with a missing neutrino
    nominal
    constrain_four_momentum
    miss_track_of :nu_e
    chi2_cut 200
  }

alg_d0bar_eta.with_decay_card(decay_card_d0bar_eta).apply(sel_d0bar_eta)

# =====================================================================
### Mode III : J/psi -> anti-D0 rho0, anti-D0 -> K+ e- anti-nu_e ###
# =====================================================================
alg_d0bar_rho0 = Algorithm.new("JpsiD0barRho0")
alg_d0bar_rho0.set_header(["JpsiD0barRho0Alg/JpsiD0barRho0.h"])
               .set_constant({"ECMS" => [:double, 3.097]})
               .set_alias({"std::vector<double>" => "Vdouble"})

sel_d0bar_rho0 = Selection.new
sel_d0bar_rho0
  .select_track {                 # four tracks, net charge zero
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"               # K+, pi+
    nChrn     "==2"               # e-, pi-
    nNet      "==0"
  }
  .select_photon {                # common photon quality cuts
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.8
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp  "==1"                    # one K+
    nlm  "==1"                    # one e-
    npip "==1"                    # one pi+
    npim "==1"                    # one pi-
  }
  .kinematic_fit([:kp, :lm, :pip, :pim]) {    # nominal fit with a missing neutrino
    nominal
    constrain_four_momentum
    miss_track_of :nu_e
    chi2_cut 200
  }

alg_d0bar_rho0.with_decay_card(decay_card_d0bar_rho0).apply(sel_d0bar_rho0)

# =====================================================================
### Mode IV : J/psi -> D- pi+, D- -> K_S0 e- anti-nu_e ###
# =====================================================================
alg_dm_pip = Algorithm.new("JpsiDmPip")
alg_dm_pip.set_header(["JpsiDmPipAlg/JpsiDmPip.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_dm_pip = Selection.new
sel_dm_pip
  .select_track {                 # K_S0(pi+pi-) + e- + direct pi+
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"               # pi+ from K_S0 and the direct pi+
    nChrn     ">=2"               # e- and pi- from K_S0
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.8
    identify :pion, against: [:kaon]
    nlm  "==1"                    # one e-
    npip ">=1"
    npim ">=1"
  }
  .secondary_vertex_fit([:pip, :pim]) {       # reconstruct K_S0 from pi+pi- vertex
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:K_S0, :lm, :pip]) {        # nominal fit with a missing neutrino
    nominal
    constrain_four_momentum
    miss_track_of :nu_e
    chi2_cut 200
  }

# K_S0 mass window and decay-length requirement have no dedicated DSL primitive
alg_dm_pip.note(:ks0_vertex_selection,
  "K_S0 candidates from the pi+pi- secondary vertex are required to satisfy " \
  "|M(pi+pi-) - m_K_S0| < 12 MeV/c^2 and a decay length > 2 sigma, on top of the " \
  "mass-difference-minimising secondary vertex fit; efficiency from these vertex cuts " \
  "evaluated on signal exclusive MC")

alg_dm_pip.with_decay_card(decay_card_dm_pip).apply(sel_dm_pip)

# =====================================================================
### Mode V : J/psi -> D- rho+, D- -> K_S0 e- anti-nu_e, rho+ -> pi+ pi0 ###
# =====================================================================
alg_dm_rhop = Algorithm.new("JpsiDmRhoP")
alg_dm_rhop.set_header(["JpsiDmRhoPAlg/JpsiDmRhoP.h"])
           .set_constant({"ECMS" => [:double, 3.097]})
           .set_alias({"std::vector<double>" => "Vdouble"})

sel_dm_rhop = Selection.new
sel_dm_rhop
  .select_track {                 # same four-track topology as Mode IV
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"     # >=2 photons for pi0 -> gamma gamma
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.8
    identify :pion, against: [:kaon]
    nlm  "==1"
    npip ">=1"
    npim ">=1"
  }
  .secondary_vertex_fit([:pip, :pim]) {       # reconstruct K_S0 from pi+pi- vertex
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # gamma gamma mass-constrained fit to pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 20
    npi0 ">=1"
  }
  .kinematic_fit([:K_S0, :lm, :pip, :pi0]) {  # nominal fit with a missing neutrino
    nominal
    constrain_four_momentum
    miss_track_of :nu_e
    chi2_cut 200
  }

alg_dm_rhop.note(:ks0_vertex_selection,
  "K_S0 candidates from the pi+pi- secondary vertex are required to satisfy " \
  "|M(pi+pi-) - m_K_S0| < 12 MeV/c^2 and a decay length > 2 sigma, on top of the " \
  "mass-difference-minimising secondary vertex fit; efficiency from these vertex cuts " \
  "evaluated on signal exclusive MC")

alg_dm_rhop.with_decay_card(decay_card_dm_rhop).apply(sel_dm_rhop)

# =====================================================================
### Execute all five algorithms on data, inclusive MC, and their signal MC ###
# =====================================================================
root_files_d0bar_pi0  = alg_d0bar_pi0.execute_on([jpsi_data, jpsi_incMC, exMC_d0bar_pi0])
root_files_d0bar_eta  = alg_d0bar_eta.execute_on([jpsi_data, jpsi_incMC, exMC_d0bar_eta])
root_files_d0bar_rho0 = alg_d0bar_rho0.execute_on([jpsi_data, jpsi_incMC, exMC_d0bar_rho0])
root_files_dm_pip     = alg_dm_pip.execute_on([jpsi_data, jpsi_incMC, exMC_dm_pip])
root_files_dm_rhop    = alg_dm_rhop.execute_on([jpsi_data, jpsi_incMC, exMC_dm_rhop])