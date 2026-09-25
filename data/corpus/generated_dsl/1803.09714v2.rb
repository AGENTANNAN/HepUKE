# ============================================================
# Datasets: psi(3686) at 3.686 GeV
# ============================================================
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(2S) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # corresponding inclusive MC

# ============================================================
# Decay cards (EvtGen format)
# ============================================================
# Mode I: psi(3686) -> eta' e+ e-, eta' -> gamma pi+ pi-
decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.0000 eta' e+ e-      PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi-   PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: psi(3686) -> eta' e+ e-, eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.0000 eta' e+ e-      PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- eta     PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma     PHSP;
    Enddecay

    End
DECAYCARD

# ============================================================
# Exclusive MC (100k events per eta' decay mode)
# ============================================================
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_etaprime_ee_modeI"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_etaprime_ee_modeII"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# ============================================================
# Mode I : psi(3686) -> eta' e+ e-, eta' -> gamma pi+ pi-
# ============================================================
alg_name_I = "EtaPrimeEeModeI"
alg_modeI  = Algorithm.new(alg_name_I)
alg_modeI.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
  .select_track {                 # 4 charged tracks: 2 positive, 2 negative, net charge 0
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        100.0               # |Vz| < 100 mm
    Vr        10.0                # Vr < 10 mm
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {                # at least one photon (Mode I)
    tdc_emc_start     0           # EMC TDC 0-700 ns
    tdc_emc_end       14
    angle_to_track    10.0        # angle to nearest charged track > 10 deg
    energyThreshold_b 0.025       # E > 25 MeV (barrel)
    energyThreshold_e 0.050       # E > 50 MeV (endcap)
    nGam              ">=1"
  }
  .pid(method: :probability) {    # PID: probability method, C.L. > 0.001
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6  # p>1.0 GeV -> lepton; EMC E>0.6 GeV -> e, else mu
    nlp ">=1"                     # at least one l+
    nlm ">=1"                     # at least one l-
  }
  .remove([:lp <= :chrgp, :lm <= :chrgn])    # leptons taken out of the charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})  # remaining charged tracks treated as pions
  .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {  # 4C fit to gamma pi+ pi- l+ l-
    nominal
    constrain_four_momentum
    chi2_cut 80
  }

alg_modeI
  .note(:background_veto, "e+e- pair vertex distance from IP delta_xy < 2 cm; " \
                          "cos(theta)(e+) < 0.8 and cos(theta)(e-) > -0.8; " \
                          "recoil mass of pi+pi- < 2.9 GeV; E/p > 0.8 on the " \
                          "higher-momentum e+- track -- BOSS-side background " \
                          "suppression cuts with no dedicated DSL method")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

root_files_modeI = alg_modeI.execute_on([psip_data, psip_incMC, exMC_modeI])

# ============================================================
# Mode II: psi(3686) -> eta' e+ e-, eta' -> pi+ pi- eta, eta -> gamma gamma
# ============================================================
alg_name_II = "EtaPrimeEeModeII"
alg_modeII  = Algorithm.new(alg_name_II)
alg_modeII.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
  .select_track {                 # same charged-track selection as Mode I
    cos_theta 0.93
    Vz        100.0
    Vr        10.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {                # at least two photons (eta -> gamma gamma)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"
  }
  .pid(method: :probability) {    # same PID as Mode I
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp ">=1"
    nlm ">=1"
  }
  .remove([:lp <= :chrgp, :lm <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {           # reconstruct eta -> gamma gamma
    invariant_mass_of(:gamma, :gamma).within(0.520, 0.575)         # gamma gamma mass window
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  .kinematic_fit([:pip, :pim, :lp, :lm, :eta]) {       # 4C fit to pi+ pi- l+ l- eta
    nominal
    constrain_four_momentum
    chi2_cut 80
  }

alg_modeII
  .note(:background_veto, "e+e- pair vertex distance from IP delta_xy < 2 cm; " \
                          "cos(theta)(e+) < 0.8 and cos(theta)(e-) > -0.8; " \
                          "recoil mass of pi+pi- < 3.2 GeV -- BOSS-side background " \
                          "suppression cuts with no dedicated DSL method")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

root_files_modeII = alg_modeII.execute_on([psip_data, psip_incMC, exMC_modeII])