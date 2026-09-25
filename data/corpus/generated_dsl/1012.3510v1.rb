# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")        # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")     # corresponding inclusive MC

# ---------------------------------------------------------------------------
# Decay cards
# ---------------------------------------------------------------------------
# Mode 1: J/psi -> gamma X(1835), X(1835) -> pi+ pi- eta', eta' -> gamma rho0, rho0 -> pi+ pi-
decay_card_mode1 = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma X(1835) PHSP;
    Enddecay

    Decay X(1835)
    1.0000 pi+ pi- eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma rho0 PHSP;
    Enddecay

    Decay rho0
    1.0000 pi+ pi- VSS;
    Enddecay

    End
DECAYCARD

# Mode 2: J/psi -> gamma X(1835), X(1835) -> pi+ pi- eta', eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_mode2 = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma X(1835) PHSP;
    Enddecay

    Decay X(1835)
    1.0000 pi+ pi- eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive MC — 500k events for each eta' decay mode
# ---------------------------------------------------------------------------
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaX1835_etap_gammarho"
  config.related_dataset = jpsi_data
  config.events         = 500000
  config.decay_card     = decay_card_mode1
  config.cross_section  = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gammaX1835_etap_pipieta"
  config.related_dataset = jpsi_data
  config.events         = 500000
  config.decay_card     = decay_card_mode2
  config.cross_section  = :default
end

### Event selection (BOSS) ###
# ===========================================================================
# Mode 1 : eta' -> gamma rho0, rho0 -> pi+ pi-
#          final state  gamma gamma pi+ pi- pi+ pi-
# ===========================================================================
alg_name_mode1 = "JpsiX1835EtaPToGammaRho"
alg_mode1 = Algorithm.new(alg_name_mode1)
alg_mode1.set_header(["#{alg_name_mode1}Alg/#{alg_name_mode1}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_mode1 = Selection.new
sel_mode1
  .select_track {                    # charged-track selection
    cos_theta 0.93                   # |cos(theta)| < 0.93
    Vz        20.0                   # |Vz| < 20 cm
    Vr        2.0                    # Vr < 2 cm
    nChrp "==2"                      # exactly two positive tracks
    nChrn "==2"                      # exactly two negative tracks
    nNet  "==0"                      # net charge zero
  }
  .select_photon {                   # photon selection
    tdc_emc_start 0                  # EMC TDC start
    tdc_emc_end   14                 # EMC TDC end
    angle_to_track 5.0               # angle to nearest charged track > 5 deg
    energyThreshold_b 0.1            # 100 MeV threshold in barrel
    energyThreshold_e 0.1            # 100 MeV threshold in endcap
    nGam ">=2"                       # mode 1: at least two photons
  }
  .pid(method: :probability) {       # PID with the probability method
    prob_cut 0.001                   # probability > 0.001
    identify :pion, against: [:kaon, :proton]   # separate pi from K and p
    npip ">=2"                       # at least two pi+
    npim ">=1"                       # at least one pi- (>=3 of 4 tracks are pions)
  }
  .kinematic_fit([:gamma, :gamma, :pip, :pip, :pim, :pim]) {   # 4C fit over full final state
    nominal                          # nominal fit: corrected four-momenta saved
    constrain_four_momentum          # 4C energy-momentum conservation
    # ---- background vetoes / windows ----
    invariant_mass_of(:gamma, :gamma).out_of(0.095, 0.175)      # |M(gg) - m_pi0| < 0.04
    invariant_mass_of(:gamma, :gamma).out_of(0.518, 0.578)      # |M(gg) - m_eta|  < 0.03
    invariant_mass_of(:gamma, :gamma).within(0.72, 0.82)        # 0.72 < M(gg) < 0.82
    invariant_mass_of(:gamma, :pip, :pim).out_of(0.541, 0.555)  # |M(g pi+ pi-) - m_eta| < 0.007
    invariant_mass_of(:pip, :pim).within(0.575, 0.975)          # rho window  |M(pi+pi-) - m_rho| < 0.2
    invariant_mass_of(:gamma, :pip, :pim).within(0.943, 0.973)  # eta' window |M(g pi+ pi-) - m_etap| < 0.015
    chi2_cut 40                      # chi2 < 40
  }

alg_mode1.with_decay_card(decay_card_mode1).apply(sel_mode1)

# ===========================================================================
# Mode 2 : eta' -> pi+ pi- eta, eta -> gamma gamma
#          final state  gamma gamma gamma pi+ pi- pi+ pi-
#          4C (four-momentum) + 1C (eta mass) = 5C fit
# ===========================================================================
alg_name_mode2 = "JpsiX1835EtaPToPiPiEta"
alg_mode2 = Algorithm.new(alg_name_mode2)
alg_mode2.set_header(["#{alg_name_mode2}Alg/#{alg_name_mode2}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

sel_mode2 = Selection.new
sel_mode2
  .select_track {                    # charged-track selection (identical to mode 1)
    cos_theta 0.93
    Vz        20.0
    Vr        2.0
    nChrp "==2"
    nChrn "==2"
    nNet  "==0"
  }
  .select_photon {                   # photon selection
    tdc_emc_start 0
    tdc_emc_end   14
    angle_to_track 5.0
    energyThreshold_b 0.1
    energyThreshold_e 0.1
    nGam ">=3"                       # mode 2: at least three photons
  }
  .pid(method: :probability) {       # PID with the probability method
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=2"                       # at least two pi+
    npim ">=1"                       # at least one pi-
  }
  # 1-C mass-constrained fit of two photons to the eta mass (reconstructs eta for subsequent 5C fit)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)   # M(gg) -> m_eta
    invariant_mass_of(:gamma, :gamma).out_of(0.095, 0.175)                 # veto pi0: |M(gg) - m_pi0| < 0.04
    chi2_cut 40                      # chi2 < 40
    neta ">=1"                       # at least one eta candidate
  }
  # main fit: 4C conservation + the eta mass constraint already imposed in the kalman step => 5C in total
  .kinematic_fit([:gamma, :pip, :pip, :pim, :pim, :eta]) {
    nominal                          # nominal fit
    constrain_four_momentum          # 4C energy-momentum conservation
    invariant_mass_of(:pip, :pim, :eta).within(0.948, 0.968)   # eta' window |M(pi+pi-eta) - m_etap| < 0.01
    chi2_cut 40                      # chi2 < 40
  }

alg_mode2
  .note(:background_veto, "|M(gamma gamma) - m_pi0| < 0.04 GeV/c^2 applied to ALL photon pairs " \
        "(the DSL form constrains the photon combination entering the fit); " \
        "|M(gamma gamma) - m_eta| < 0.03 GeV/c^2 used as eta preselection window")
  .with_decay_card(decay_card_mode2).apply(sel_mode2)

# ---------------------------------------------------------------------------
# Execute both algorithms on data, inclusive MC and the corresponding signal MC
# ---------------------------------------------------------------------------
root_files_mode1 = alg_mode1.execute_on([jpsi_data, jpsi_incMC, exMC_mode1])
root_files_mode2 = alg_mode2.execute_on([jpsi_data, jpsi_incMC, exMC_mode2])