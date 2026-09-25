# Core DSL classes are loaded automatically at execution time.
#
# Analysis: five semileptonic("tagged")D-meson modes at the J/psi (3.097 GeV)
#   J/psi -> D0 pi0, D0 eta, D0 rho0, D- pi+, D- rho+   (charge conjugates implied)
#   D0 -> K+ e- anti-nu_e   (semileptonic D0)
#   D- -> K_S0 e- anti-nu_e (semileptonic D-)
#   pi0/eta -> gamma gamma, rho0 -> pi+ pi-, rho+ -> pi+ pi0, K_S0 -> pi+ pi-

### Datasets ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # matching inclusive MC

### Decay cards (EvtGen format) ###
decay_card_D0pi0 = <<~DECAYCARD
    Decay J/psi
    1.0000 D0 pi0 PHSP;
    Enddecay

    Decay D0
    1.0000 K+ e- anti-nu_e PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_D0eta = <<~DECAYCARD
    Decay J/psi
    1.0000 D0 eta PHSP;
    Enddecay

    Decay D0
    1.0000 K+ e- anti-nu_e PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_D0rho0 = <<~DECAYCARD
    Decay J/psi
    1.0000 D0 rho0 PHSP;
    Enddecay

    Decay D0
    1.0000 K+ e- anti-nu_e PHSP;
    Enddecay

    Decay rho0
    1.0000 pi+ pi- VSS;
    Enddecay

    End
DECAYCARD

decay_card_Dmpip = <<~DECAYCARD
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

decay_card_Dmrho = <<~DECAYCARD
    Decay J/psi
    1.0000 D- rho+ PHSP;
    Enddecay

    Decay D-
    1.0000 K_S0 e- anti-nu_e PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay rho+
    1.0000 pi+ pi0 VSS;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive signal MC — 500k events each, all at E_cms = 3.097 GeV ###
exmc_D0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_D0pi0_semilep"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_D0pi0
  config.cross_section   = :default
end

exmc_D0eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_D0eta_semilep"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_D0eta
  config.cross_section   = :default
end

exmc_D0rho0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_D0rho0_semilep"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_D0rho0
  config.cross_section   = :default
end

exmc_Dmpip = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_Dmpip_semilep"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_Dmpip
  config.cross_section   = :default
end

exmc_Dmrho = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_Dmrho_semilep"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_Dmrho
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ---- common selections -------------------------------------------------

# Two-track modes (D0 pi0, D0 eta): exactly one K+ and one e-, plus >= 2 photons
sel_2trk = Selection.new
  .select_track {
      cos_theta        0.93
      Vz               10.0
      Vr               1.0
      nChrp            "==1"
      nChrn            "==1"
      nNet             "==0"
  }
  .select_photon {
      tdc_emc_start    0
      tdc_emc_end      14
      angle_to_track   10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam             ">=2"
  }

# Four-track modes without photons from the signal (D0 rho0, D- pi+)
sel_4trk = Selection.new
  .select_track {
      cos_theta        0.93
      Vz               10.0
      Vr               1.0
      nChrp            "==2"
      nChrn            "==2"
      nNet             "==0"
  }
  .select_photon {
      tdc_emc_start    0
      tdc_emc_end      14
      angle_to_track   10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
  }

# Four-track mode with a pi0 from rho+ -> pi+ pi0 (D- rho+)
sel_4trk_ph = Selection.new
  .select_track {
      cos_theta        0.93
      Vz               10.0
      Vr               1.0
      nChrp            "==2"
      nChrn            "==2"
      nNet             "==0"
  }
  .select_photon {
      tdc_emc_start    0
      tdc_emc_end      14
      angle_to_track   10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam             ">=2"
  }

# ---- Mode 1: J/psi -> D0 pi0, D0 -> K+ e- anti-nu_e, pi0 -> gamma gamma ----
alg_D0pi0 = Algorithm.new("JpsiD0pi0Semilep")
alg_D0pi0.set_header(["JpsiD0pi0SemilepAlg/JpsiD0pi0Semilep.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_D0pi0 = sel_2trk.dup
  .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion]                 # K+ vs pion
      identify :pion, against: [:kaon, :proton]        # pion vs kaon and proton
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6   # e-/mu separation
      nkp "==1"
      nlm "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)   # pi0 mass window
      chi2_cut 20
      npi0 ">=1"
  }
  .kinematic_fit([:kp, :lm, :nu_e, :pi0]) {
      nominal
      miss_track_of :nu_e                 # undetected neutrino allowed missing in the fit
      constrain_four_momentum
      chi2_cut 200
  }

alg_D0pi0
  .note(:electron_ep_ratio_cut, "electron candidate required to satisfy E/p > 0.8; ep_ratio_of() not yet available in the DSL")
  .note(:missing_neutrino, "undetected neutrino handled by missing four-momentum U_miss = E_miss - c|p_miss| built from the recoil against the reconstructed system; best neutrino candidate chosen by mass")
  .note(:d_recoil_mass, "nominal D-meson recoil-mass window 1.80-1.95 GeV/c^2 imposed against the pi0 system to suppress background")
  .with_decay_card(decay_card_D0pi0)
  .apply(sel_D0pi0)

alg_D0pi0.execute_on([jpsi_data, jpsi_incMC, exmc_D0pi0])

# ---- Mode 2: J/psi -> D0 eta, D0 -> K+ e- anti-nu_e, eta -> gamma gamma ----
alg_D0eta = Algorithm.new("JpsiD0etaSemilep")
alg_D0eta.set_header(["JpsiD0etaSemilepAlg/JpsiD0etaSemilep.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_D0eta = sel_2trk.dup
  .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion]
      identify :pion, against: [:kaon, :proton]
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      nkp "==1"
      nlm "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      invariant_mass_of(:gamma, :gamma).within(0.50, 0.57)     # eta mass window
      chi2_cut 20
      neta ">=1"
  }
  .kinematic_fit([:kp, :lm, :nu_e, :eta]) {
      nominal
      miss_track_of :nu_e
      constrain_four_momentum
      chi2_cut 200
  }

alg_D0eta
  .note(:electron_ep_ratio_cut, "electron candidate required to satisfy E/p > 0.8; ep_ratio_of() not yet available in the DSL")
  .note(:missing_neutrino, "undetected neutrino handled by missing four-momentum U_miss = E_miss - c|p_miss|; best neutrino candidate chosen by mass")
  .note(:d_recoil_mass, "nominal D-meson recoil-mass window 1.80-1.95 GeV/c^2 imposed against the eta system to suppress background")
  .with_decay_card(decay_card_D0eta)
  .apply(sel_D0eta)

alg_D0eta.execute_on([jpsi_data, jpsi_incMC, exmc_D0eta])

# ---- Mode 3: J/psi -> D0 rho0, D0 -> K+ e- anti-nu_e, rho0 -> pi+ pi- ----
alg_D0rho0 = Algorithm.new("JpsiD0rho0Semilep")
alg_D0rho0.set_header(["JpsiD0rho0SemilepAlg/JpsiD0rho0Semilep.h"])
          .set_constant({"ECMS" => [:double, 3.097]})

sel_D0rho0 = sel_4trk.dup
  .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion]
      identify :pion, against: [:kaon, :proton]
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      nkp "==1"
      nlm "==1"
      npip "==1"
      npim "==1"
  }
  .kinematic_fit([:kp, :lm, :nu_e, :pip, :pim]) {
      nominal
      miss_track_of :nu_e
      constrain_four_momentum
      chi2_cut 200
  }

alg_D0rho0
  .note(:electron_ep_ratio_cut, "electron candidate required to satisfy E/p > 0.8; ep_ratio_of() not yet available in the DSL")
  .note(:missing_neutrino, "undetected neutrino handled by missing four-momentum U_miss = E_miss - c|p_miss|; best neutrino candidate chosen by mass")
  .note(:d_recoil_mass, "nominal D-meson recoil-mass window 1.80-1.95 GeV/c^2 imposed against the rho0 system to suppress background")
  .with_decay_card(decay_card_D0rho0)
  .apply(sel_D0rho0)

alg_D0rho0.execute_on([jpsi_data, jpsi_incMC, exmc_D0rho0])

# ---- Mode 4: J/psi -> D- pi+, D- -> K_S0 e- anti-nu_e, K_S0 -> pi+ pi- ----
alg_Dmpip = Algorithm.new("JpsiDmpipSemilep")
alg_Dmpip.set_header(["JpsiDmpipSemilepAlg/JpsiDmpipSemilep.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_Dmpip = sel_4trk.dup
  .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion]
      identify :pion, against: [:kaon, :proton]
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      nlm "==1"
      npip "==2"     # one K_S0 daughter pi+ and one signal pi+
      npim "==1"     # K_S0 daughter pi-
  }
  .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_verfit_chi2   # best chi2 combination
      remove_used_particle_from_candidate_list                  # remaining pi+ is the signal track
  }
  .kinematic_fit([:K_S0, :lm, :nu_e, :pip]) {
      nominal
      miss_track_of :nu_e
      constrain_four_momentum
      chi2_cut 200
  }

alg_Dmpip
  .note(:k_s0_selection, "K_S0 secondary-vertex requirements: |M(pi+pi-)-m_K_S0| < 12 MeV, decay length > 2 sigma from the IP, K_S0 daughter |Vz| < 20 cm, best chi2 combination")
  .note(:electron_ep_ratio_cut, "electron candidate required to satisfy E/p > 0.8; ep_ratio_of() not yet available in the DSL")
  .note(:missing_neutrino, "undetected neutrino handled by missing four-momentum U_miss = E_miss - c|p_miss|; best neutrino candidate chosen by mass")
  .note(:d_recoil_mass, "nominal D-meson recoil-mass window 1.80-1.95 GeV/c^2 imposed against the pi+ system to suppress background")
  .with_decay_card(decay_card_Dmpip)
  .apply(sel_Dmpip)

alg_Dmpip.execute_on([jpsi_data, jpsi_incMC, exmc_Dmpip])

# ---- Mode 5: J/psi -> D- rho+, D- -> K_S0 e- anti-nu_e, rho+ -> pi+ pi0 ----
alg_Dmrho = Algorithm.new("JpsiDmrhoSemilep")
alg_Dmrho.set_header(["JpsiDmrhoSemilepAlg/JpsiDmrhoSemilep.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_Dmrho = sel_4trk_ph.dup
  .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion]
      identify :pion, against: [:kaon, :proton]
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      nlm "==1"
      npip "==2"     # one K_S0 daughter pi+ and one rho+ daughter pi+
      npim "==1"     # K_S0 daughter pi-
  }
  .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_verfit_chi2
      remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)
      chi2_cut 20
      npi0 ">=1"
  }
  .kinematic_fit([:K_S0, :lm, :nu_e, :pip, :pi0]) {
      nominal
      miss_track_of :nu_e
      constrain_four_momentum
      chi2_cut 200
  }

alg_Dmrho
  .note(:k_s0_selection, "K_S0 secondary-vertex requirements: |M(pi+pi-)-m_K_S0| < 12 MeV, decay length > 2 sigma from the IP, K_S0 daughter |Vz| < 20 cm, best chi2 combination")
  .note(:electron_ep_ratio_cut, "electron candidate required to satisfy E/p > 0.8; ep_ratio_of() not yet available in the DSL")
  .note(:missing_neutrino, "undetected neutrino handled by missing four-momentum U_miss = E_miss - c|p_miss|; best neutrino candidate chosen by mass")
  .note(:d_recoil_mass, "nominal D-meson recoil-mass window 1.80-1.95 GeV/c^2 imposed against the pi+ pi0 (rho+) system to suppress background")
  .with_decay_card(decay_card_Dmrho)
  .apply(sel_Dmrho)

alg_Dmrho.execute_on([jpsi_data, jpsi_incMC, exmc_Dmrho])