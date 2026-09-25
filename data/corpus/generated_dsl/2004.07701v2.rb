### Dataset description ###
jpsi_data = DatasetManager.real_data.find("708_3097")        # J/psi real data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive MC
psip_data = DatasetManager.real_data.find("709_3686")        # psi(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # Corresponding inclusive MC

# Decay card: J/psi -> Sigma+ Sigma_bar- ; Sigma+ -> p pi0 ; Sigma_bar- -> anti-p pi0 ; pi0 -> gamma gamma
decay_card_jpsi = <<~DECAYCARD
    Decay J/psi
    1.0000 Sigma+ anti-Sigma-    PHSP;
    Enddecay

    Decay Sigma+
    1.0000 p+ pi0    PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000 anti-p- pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: psi(2S) -> Sigma+ Sigma_bar- with the same subsequent decays
decay_card_psip = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Sigma+ anti-Sigma-    PHSP;
    Enddecay

    Decay Sigma+
    1.0000 p+ pi0    PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000 anti-p- pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 100k events for the J/psi channel
exMC_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_jpsi_sigmasigmabar"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_jpsi
  config.cross_section = :default
end

# Exclusive MC: 100k events for the psi(2S) channel
exMC_psip = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_psip_sigmasigmabar"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card_psip
  config.cross_section = :default
end

### Event selection (BOSS) ###
# ------------------------------------------------------------------
# Analysis I: e+e- -> J/psi -> Sigma+ Sigma_bar-
# ------------------------------------------------------------------
alg_name_jpsi = "JpsiToSigmaSigma"
alg_jpsi = Algorithm.new(alg_name_jpsi)
alg_jpsi.set_header(["#{alg_name_jpsi}Alg/#{alg_name_jpsi}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})
        # Inexpressible: the p p_bar pi0 pi0 combination is chosen by minimizing the
        # mass-difference metric sqrt((M(p pi0)-m_Sigma+)^2 + (M(p_bar pi0)-m_Sigma_bar-)^2).
        .note(:best_combination, "best p p_bar pi0 pi0 combination chosen by minimizing sqrt((M(p pi0)-m_Sigma+)^2 + (M(p_bar pi0)-m_Sigma_bar-)^2)")

event_selection_jpsi = Selection.new
  .select_track {                 # exactly one positive and one negative charged track
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz 10.0                       # |Vz| < 10 cm
    Vr 2.0                        # Vr < 2 cm
    nChrp "==1"                   # exactly one positive track
    nChrn "==1"                   # exactly one negative track
    nNet "==0"                    # net charge zero
  }
  .select_photon {                # at least four good photons
    tdc_emc_start 0               # TDC start 0
    tdc_emc_end 14                # TDC end 14
    angle_to_track 10.0           # angle to nearest charged track > 10 degrees
    energyThreshold_b 0.025       # E > 25 MeV in EMC barrel
    energyThreshold_e 0.050       # E > 50 MeV in EMC endcap
    nGam ">=4"                    # at least four photons
  }
  .pid(method: :probability) {    # proton / anti-proton identification
    prob_cut 0.001                # PID probability > 0.001
    identify :proton, against: [:kaon, :pion]  # id p+ and anti-p- vs K and pi
    nprp "==1"                    # exactly one proton
    nprm "==1"                    # exactly one anti-proton
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])   # remove proton tracks from charged-track lists
  .kalman_kinematic_fit([:gamma, :gamma]) {    # first pi0 -> gamma gamma (1C Kalman fit)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"                    # at least one pi0 from this fit
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {    # second pi0 -> gamma gamma (1C Kalman fit)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"                    # at least one pi0 from this fit
  }
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) {   # 4C kinematic fit to p p_bar pi0 pi0
    nominal
    constrain_four_momentum
    chi2_cut 100
  }

alg_jpsi.with_decay_card(decay_card_jpsi).apply(event_selection_jpsi)
alg_jpsi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi])

# ------------------------------------------------------------------
# Analysis II: e+e- -> psi(2S) -> Sigma+ Sigma_bar-
# (additional veto against psi(2S) -> pi0 pi0 J/psi)
# ------------------------------------------------------------------
alg_name_psip = "PsipToSigmaSigma"
alg_psip = Algorithm.new(alg_name_psip)
alg_psip.set_header(["#{alg_name_psip}Alg/#{alg_name_psip}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .note(:best_combination, "best p p_bar pi0 pi0 combination chosen by minimizing sqrt((M(p pi0)-m_Sigma+)^2 + (M(p_bar pi0)-m_Sigma_bar-)^2)")
        .note(:background_veto, "veto |M(p p_bar) - 3.1 GeV| > 0.05 GeV to suppress psi(2S) -> pi0 pi0 J/psi")

event_selection_psip = Selection.new
  .select_track {                 # exactly one positive and one negative charged track
    cos_theta 0.93
    Vz 10.0
    Vr 2.0
    nChrp "==1"
    nChrn "==1"
    nNet "==0"
  }
  .select_photon {                # at least four good photons
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=4"
  }
  .pid(method: :probability) {    # proton / anti-proton identification
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"
    nprm "==1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])   # remove proton tracks from charged-track lists
  .kalman_kinematic_fit([:gamma, :gamma]) {    # first pi0 -> gamma gamma (1C Kalman fit)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {    # second pi0 -> gamma gamma (1C Kalman fit)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) {   # 4C kinematic fit to p p_bar pi0 pi0
    nominal
    constrain_four_momentum
    invariant_mass_of(:prp, :prm).out_of(3.05, 3.15)  # veto |M(p p_bar) - 3.1 GeV| > 0.05 GeV
    chi2_cut 100
  }

alg_psip.with_decay_card(decay_card_psip).apply(event_selection_psip)
alg_psip.execute_on([psip_data, psip_incMC, exMC_psip])