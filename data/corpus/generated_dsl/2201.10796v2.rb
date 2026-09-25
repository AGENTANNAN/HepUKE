# ============================================================
# Datasets — J/psi(3097), 10.09 x 10^9 events
# ============================================================
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ============================================================
# Decay cards (all phase-space), one per eta' decay mode
# ============================================================
# Mode I: J/psi -> gamma pi+ pi- eta',  eta' -> gamma pi+ pi-
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma pi+ pi- eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: J/psi -> gamma pi+ pi- eta',  eta' -> pi+ pi- eta,  eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma pi+ pi- eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ============================================================
# Exclusive MC — 200k events per mode
# ============================================================
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_etap_gam_pip_pim"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_etap_pip_pim_eta"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# ============================================================
# Mode I : J/psi -> gamma pi+ pi- eta',  eta' -> gamma pi+ pi-
# Final state: gamma gamma pi+ pi- pi+ pi-
# ============================================================
alg_name_modeI = "JpsiEtapGamPipPim"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})
         .set_alias({"std::vector<double>" => "Vdouble"})

# Charged-track / photon / PID requirements are identical to Mode II
sel_modeI = Selection.new
  .select_track {
    cos_theta  0.93
    Vz         10.0          # |Vz| < 10 cm
    Vr         1.0           # Vr < 1 cm
    nChrp      ">=2"         # at least 2 positive tracks
    nChrn      ">=2"         # at least 2 negative tracks
    nNet       "==0"         # net charge zero (>=4 tracks in total)
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0   # > 10 deg from any charged track
    energyThreshold_b 0.1    # E > 100 MeV (barrel)
    energyThreshold_e 0.1    # E > 100 MeV (endcap)
    nGam              ">=2"  # Mode I needs >= 2 photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=1"
    npim ">=1"
  }
  .kinematic_fit([:gamma, :gamma, :pip, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum                                    # 4C fit to gamma gamma pi+ pi- pi+ pi-
    # select the gamma pi+ pi- combination closest to the eta' mass (|M - m_eta'| < 15 MeV)
    invariant_mass_of(:gamma, :pip, :pim).between(0.943, 0.973)
    # veto gamma-gamma pairs near pi0 (40 MeV), near eta (30 MeV), and in the omega region
    invariant_mass_of(:gamma, :gamma).out_of(0.095, 0.175)
    invariant_mass_of(:gamma, :gamma).out_of(0.518, 0.578)
    invariant_mass_of(:gamma, :gamma).out_of(0.720, 0.820)
    # veto gamma pi+ pi- combinations in the 400-563 MeV window
    invariant_mass_of(:gamma, :pip, :pim).out_of(0.400, 0.563)
    # veto gamma-gamma pairs near pi0 (15 MeV)
    invariant_mass_of(:gamma, :gamma).out_of(0.120, 0.150)
    chi2_cut 40
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

# ============================================================
# Mode II : J/psi -> gamma pi+ pi- eta', eta' -> pi+ pi- eta, eta -> gamma gamma
# Final state: gamma gamma gamma pi+ pi- pi+ pi-
# ============================================================
alg_name_modeII = "JpsiEtapPipPimEta"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:background_veto, "all gamma-gamma pairs with |M(gg) - m_pi0| < 40 MeV/c^2 " \
                "(eta-reconstruction stage) and within 15 MeV/c^2 of m_pi0 (eta' stage) are vetoed " \
                "to suppress pi0 -> gamma gamma contamination")

# Charged-track / photon / PID requirements identical to Mode I (photon count raised to >=3)
sel_modeII = Selection.new
  .select_track {
    cos_theta  0.93
    Vz         10.0
    Vr         1.0
    nChrp      ">=2"
    nChrn      ">=2"
    nNet       "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.1
    energyThreshold_e 0.1
    nGam              ">=3"  # Mode II needs >= 3 photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=1"
    npim ">=1"
  }
  # reconstruct eta from gamma gamma with a 1C Kalman fit to the nominal eta mass
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  }
  # nominal 5C fit (4C plus the eta mass carried through the constrained eta)
  .kinematic_fit([:gamma, :eta, :pip, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum
    # select the pi+ pi- eta combination closest to the eta' mass (|M - m_eta'| < 10 MeV)
    invariant_mass_of(:pip, :pim, :eta).between(0.948, 0.968)
    chi2_cut 40
  }
  # additional 4C fit without the eta mass constraint, for background suppression
  .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pip, :pim, :pim]) {
    constrain_four_momentum
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

# ============================================================
# Execute both algorithms
# ============================================================
alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])
alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])