### dataset description ###
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for Mode I: J/psi -> e+ e- eta, eta -> gamma gamma
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000 e+ e- eta    PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma   PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for Mode II: J/psi -> e+ e- eta, eta -> pi+ pi- pi0
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.0000 e+ e- eta    PHSP;
    Enddecay

    Decay eta
    1.0000 pi+ pi- pi0  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_jpsi_ee_eta_gg"
    config.related_dataset = jpsi_data
    config.events = 100000
    config.decay_card = decay_card_modeI
    config.cross_section = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_jpsi_ee_eta_3pi"
    config.related_dataset = jpsi_data
    config.events = 100000
    config.decay_card = decay_card_modeII
    config.cross_section = :default
end

### Mode I: J/psi -> e+ e- eta, eta -> gamma gamma ###
alg_modeI = Algorithm.new("JpsiEEtaGG")
alg_modeI.set_header(["JpsiEEtaGGAlg/JpsiEEtaGG.h"])
          .set_constant({"ECMS" => [:double, 3.097]})

sel_modeI = Selection.new
sel_modeI.select_track {
              cos_theta  0.93
              Vz         10.0
              Vr         1.0
              nChrp      "==1"
              nChrn      "==1"
              nNet       "==0"
          }
          .select_photon {
              tdc_emc_start     0
              tdc_emc_end       14
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              angle_to_track    10.0
              nGam              ">=2"
          }
          .pid(method: :probability) {
              identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                             treat_as_electron_if_energy_above: 0.6
              nlp "==1"; nlm "==1"
          }
          .kinematic_fit([:lp, :lm, :gamma, :gamma]) {
              nominal
              vertex_fit([0, 1])
              constrain_four_momentum
              chi2_cut 100
          }

alg_modeI
  .note(:photon_conversion_veto, "delta_xy < 2 cm requirement to suppress gamma conversion background from J/psi -> gamma eta; uses photon-conversion finder algorithm; removes ~98% of conversion events while retaining ~80% of signal")
  .note(:cos_theta_heli_cut, "|cos(theta_heli)| < 0.9 to suppress non-peaking QED background from e+e- -> e+e- gamma(gamma) and e+e- -> 3gamma processes in the eta -> gamma gamma decay mode")
  .note(:resonance_veto, "veto regions 0.65 < m(e+e-) < 0.90 GeV/c^2 and 0.96 < m(e+e-) < 1.08 GeV/c^2 to suppress J/psi -> V eta (V=rho,omega,phi) peaking background with V -> e+e- or V -> pi+pi- mis-ID")
  .note(:electron_momentum_cut, "p(e±) < 1.45 GeV/c for m(e+e-) > 0.5 GeV/c^2 to suppress radiative Bhabha background e+e- -> gamma e+e- in high m(e+e-) region")
  .note(:eta_mass_window, "eta candidate reconstructed from gamma gamma pair with m(gamma gamma) in [0.45, 0.65] GeV/c^2; best gamma-gamma pair chosen by minimizing chi2_4C; eta mass window applied to 4C-fit-corrected momenta")
  .note(:tff_model, "signal MC generated with TFF following modified multipole function (Eq.2) with Lambda=2.56 GeV/c^2; angular distribution parameterized as (1+alpha_theta*cos^2(theta_eta)) with alpha_theta=1.0 measured from data")
  .note(:dark_photon_search, "search for dark photon gamma' in J/psi -> gamma' eta, gamma' -> e+e- performed via unbinned ML fits to m(e+e-) distribution in steps of 2 MeV/c^2; 90% CL upper limits on coupling strength epsilon set in range 10^-2 - 10^-3 for 0.01 <= m_gamma' <= 2.4 GeV/c^2")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

### Mode II: J/psi -> e+ e- eta, eta -> pi+ pi- pi0 ###
alg_modeII = Algorithm.new("JpsiEEta3Pi")
alg_modeII.set_header(["JpsiEEta3PiAlg/JpsiEEta3Pi.h"])
           .set_constant({"ECMS" => [:double, 3.097]})

sel_modeII = Selection.new
sel_modeII.select_track {
               cos_theta  0.93
               Vz         10.0
               Vr         1.0
               nChrp      "==2"
               nChrn      "==2"
               nNet       "==0"
           }
           .select_photon {
               tdc_emc_start     0
               tdc_emc_end       14
               energyThreshold_b 0.025
               energyThreshold_e 0.050
               angle_to_track    10.0
               nGam              ">=2"
           }
           .pid(method: :probability) {
               identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                              treat_as_electron_if_energy_above: 0.6
               nlp "==1"; nlm "==1"
           }
           .remove([:lp <= :chrgp, :lm <= :chrgn])
           .assign({chrgp: :pip, chrgn: :pim})
           .kinematic_fit([:lp, :lm, :pip, :pim, :gamma, :gamma]) {
               nominal
               vertex_fit([0, 1, 2, 3])
               constrain_four_momentum
               chi2_cut 100
           }

alg_modeII
  .note(:photon_conversion_veto, "delta_xy < 2 cm requirement to suppress gamma conversion background from J/psi -> gamma eta; uses photon-conversion finder algorithm; removes ~98% of conversion events while retaining ~80% of signal")
  .note(:resonance_veto, "veto regions 0.65 < m(e+e-) < 0.90 GeV/c^2 and 0.96 < m(e+e-) < 1.08 GeV/c^2 to suppress J/psi -> V eta (V=rho,omega,phi) peaking background with V -> e+e- or V -> pi+pi- mis-ID")
  .note(:pi0_mass_window, "pi0 candidate reconstructed from gamma gamma pair with m(gamma gamma) in [0.08, 0.16] GeV/c^2; applied to 4C-fit-corrected momenta")
  .note(:eta_mass_window, "eta candidate reconstructed from pi+ pi- pi0 with m(pi+ pi- pi0) in [0.45, 0.65] GeV/c^2; best gamma-gamma pi0 candidate chosen by minimizing chi2_4C; applied to 4C-fit-corrected momenta")
  .note(:tff_model, "signal MC generated with TFF following modified multipole function (Eq.2) with Lambda=2.56 GeV/c^2; angular distribution parameterized as (1+alpha_theta*cos^2(theta_eta)) with alpha_theta=1.0 measured from data")
  .note(:dark_photon_search, "search for dark photon gamma' in J/psi -> gamma' eta, gamma' -> e+e- performed via unbinned ML fits to m(e+e-) distribution in steps of 2 MeV/c^2; 90% CL upper limits on coupling strength epsilon set in range 10^-2 - 10^-3 for 0.01 <= m_gamma' <= 2.4 GeV/c^2")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

alg_modeI.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])
alg_modeII.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])