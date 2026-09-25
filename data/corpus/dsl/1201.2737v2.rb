### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi real data (225.2e6 events)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # J/psi inclusive MC (Lund-Charm, 2.25e8)

### Decay cards ###

# Mode I signal: J/psi -> gamma eta(1405), eta(1405) -> f0(980) pi0, f0(980) -> pi+ pi-
decay_card_ch = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta(1405)                                 PHSP;
    Enddecay

    Decay eta(1405)
    1.0000 f_0 pi0                                         PHSP;
    Enddecay

    Decay f_0
    1.0000 pi+ pi-                                         PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                                     PHSP;
    Enddecay

    End
DECAYCARD

# Mode II signal: J/psi -> gamma eta(1405), eta(1405) -> f0(980) pi0, f0(980) -> pi0 pi0
decay_card_ne = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta(1405)                                 PHSP;
    Enddecay

    Decay eta(1405)
    1.0000 f_0 pi0                                         PHSP;
    Enddecay

    Decay f_0
    1.0000 pi0 pi0                                         PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                                     PHSP;
    Enddecay

    End
DECAYCARD

# eta' decays for measurement (charged mode)
decay_card_etap_ch = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'                                      PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- pi0                                     PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                                     PHSP;
    Enddecay

    End
DECAYCARD

# eta' decays for measurement (neutral mode)
decay_card_etap_ne = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'                                      PHSP;
    Enddecay

    Decay eta'
    1.0000 pi0 pi0 pi0                                     PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                                     PHSP;
    Enddecay

    End
DECAYCARD

# Peaking background for charged eta' mode: J/psi -> gamma eta', eta' -> gamma rho0 (rho0 -> pi+ pi-)
decay_card_bg_etap_grho = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'                                      PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma rho0                                      PHSP;
    Enddecay

    Decay rho0
    1.0000 pi+ pi-                                         VSS;
    Enddecay

    End
DECAYCARD

# Peaking background for charged eta' mode: eta' -> gamma omega (omega -> pi+ pi- pi0)
decay_card_bg_etap_gomega = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'                                      PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma omega                                     PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0                                     OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                                     PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples ###
exMC_ch          = DatasetManager.create_exclusive_mc { |c| c.sample_name="jpsi_eta1405_f0pi0_ch";      c.related_dataset=jpsi_data; c.events=200000; c.decay_card=decay_card_ch;          c.cross_section=:default }
exMC_ne          = DatasetManager.create_exclusive_mc { |c| c.sample_name="jpsi_eta1405_f0pi0_ne";      c.related_dataset=jpsi_data; c.events=200000; c.decay_card=decay_card_ne;          c.cross_section=:default }
exMC_etap_ch     = DatasetManager.create_exclusive_mc { |c| c.sample_name="jpsi_gamma_etap_pipipi0";    c.related_dataset=jpsi_data; c.events=200000; c.decay_card=decay_card_etap_ch;     c.cross_section=:default }
exMC_etap_ne     = DatasetManager.create_exclusive_mc { |c| c.sample_name="jpsi_gamma_etap_3pi0";       c.related_dataset=jpsi_data; c.events=200000; c.decay_card=decay_card_etap_ne;     c.cross_section=:default }
exMC_bg_etap_grho   = DatasetManager.create_exclusive_mc { |c| c.sample_name="jpsi_gamma_etap_gamma_rho0"; c.related_dataset=jpsi_data; c.events=100000; c.decay_card=decay_card_bg_etap_grho;   c.cross_section=:default }
exMC_bg_etap_gomega = DatasetManager.create_exclusive_mc { |c| c.sample_name="jpsi_gamma_etap_gamma_omega"; c.related_dataset=jpsi_data; c.events=100000; c.decay_card=decay_card_bg_etap_gomega; c.cross_section=:default }


### Mode I: J/psi -> gamma pi+ pi- pi0 ###
alg_ch = Algorithm.new("JpsiGammaPiPiPi0")
alg_ch.set_header(["JpsiGammaPiPiPi0Alg/JpsiGammaPiPiPi0.h"])
      .set_constant({"ECMS" => [:double, 3.097]})

sel_ch = Selection.new
sel_ch.select_track {
        nChrp   "==1"
        nChrn   "==1"
        nNet    "==0"
        cos_theta 0.93
        Vz      20.0
        Vr      2.0
      }
      .select_photon {
        nGam                ">=3"
        energyThreshold_b   0.025
        energyThreshold_e   0.050
        angle_to_track      10.0        # >10 deg from any charged track
        tdc_emc_start       0
        tdc_emc_end         14
      }
      .pid(method: :probability) {
        prob_cut  0.001
        identify :pion, against: [:kaon, :proton]  # highest-CL hypothesis assignment
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      }
      .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 30                   # chi^2_4C < 30 for the 3-photon hypothesis
      }
      .kinematic_fit([:gamma, :gamma, :pip, :pim]) {   # competing 2-photon hypothesis
        constrain_four_momentum
      }
      .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim]) {  # competing 4-photon hypothesis
        constrain_four_momentum
      }

alg_ch.note(:pi0_mass_window,
            "select pi0 candidates via |M(gamma gamma) - m_pi0| < 0.015 GeV/c^2.")
      .note(:omega_veto,
            "reject events with |M(gamma pi0) - m_omega| < 0.05 GeV/c^2 to suppress J/psi -> omega pi+ pi- background.")

alg_ch.with_decay_card(decay_card_ch).apply(sel_ch)
alg_ch.execute_on([jpsi_data, jpsi_incMC, exMC_ch, exMC_etap_ch, exMC_bg_etap_grho, exMC_bg_etap_gomega])


### Mode II: J/psi -> gamma 3 pi0 ###
alg_ne = Algorithm.new("JpsiGamma3Pi0")
alg_ne.set_header(["JpsiGamma3Pi0Alg/JpsiGamma3Pi0.h"])
      .set_constant({"ECMS" => [:double, 3.097]})

sel_ne = Selection.new
sel_ne.select_track {
        nChrp   "==0"
        nChrn   "==0"
        nNet    "==0"
      }
      .select_photon {
        nGam                ">=7"     # >=7 and <9 good photons
        energyThreshold_b   0.025
        energyThreshold_e   0.050
        angle_to_track      10.0
        tdc_emc_start       0
        tdc_emc_end         14
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25                   # single-pi0 1-C fit chi^2 < 25
        npi0 ">=3"                    # at least three pi0 candidates
      }
      .kinematic_fit([:gamma, :pi0, :pi0, :pi0]) {
        nominal
        constrain_four_momentum       # 4C total energy-momentum + 3 pi0 mass constraints = 7C fit
        chi2_cut 60
      }

alg_ne.note(:pi0_decay_angle,
            "for each pi0 candidate, |cos(theta_decay)| < 0.95 where theta_decay is the polar angle of a daughter photon in the pi0 rest frame; suppresses wrong photon combinations.")
      .note(:nphoton_upper_limit,
            "require the total number of good photons to be less than 9 (7 <= nGam < 9).")
      .note(:omega_veto,
            "reject events with |M(gamma pi0) - m_omega| < 0.05 GeV/c^2 to suppress J/psi -> omega pi0 pi0 background.")

alg_ne.with_decay_card(decay_card_ne).apply(sel_ne)
alg_ne.execute_on([jpsi_data, jpsi_incMC, exMC_ne, exMC_etap_ne])
