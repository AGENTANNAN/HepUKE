# Observation of electromagnetic Dalitz decays J/psi -> P e+ e- (P = eta'/eta/pi0)
#   [arXiv:1403.7042]
#
# Five independent reconstruction modes, each with its own final state and
# kinematic-fit hypothesis, therefore one Algorithm per mode (Rule T1).

### Dataset description ###
jpsi_data     = DatasetManager.real_data.find("708_3097")   # (225.3 +/- 2.8) x 10^6 J/psi events
jpsi_incMC    = DatasetManager.inclusive_mc.find("708_3097")
cont_data_3773 = DatasetManager.real_data.find("712_3773")  # ~2.9 fb^-1 at 3.773 GeV, continuum study
cont_incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ---------------------------------------------------------------- decay cards
# Mode I: J/psi -> eta' e+ e-, eta' -> gamma pi+ pi-
decay_card_modeI = <<~DECAYCARD
  Decay J/psi
  1.0000 eta' e+ e-  PHSP;
  Enddecay

  Decay eta'
  1.0000 gamma pi+ pi-  PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: J/psi -> eta' e+ e-, eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
  Decay J/psi
  1.0000 eta' e+ e-  PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- eta  PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma  PHSP;
  Enddecay

  End
DECAYCARD

# Mode III: J/psi -> eta e+ e-, eta -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card_modeIII = <<~DECAYCARD
  Decay J/psi
  1.0000 eta e+ e-  PHSP;
  Enddecay

  Decay eta
  1.0000 pi+ pi- pi0  PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay

  End
DECAYCARD

# Mode IV: J/psi -> eta e+ e-, eta -> gamma gamma
decay_card_modeIV = <<~DECAYCARD
  Decay J/psi
  1.0000 eta e+ e-  PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma  PHSP;
  Enddecay

  End
DECAYCARD

# Mode V: J/psi -> pi0 e+ e-, pi0 -> gamma gamma
decay_card_modeV = <<~DECAYCARD
  Decay J/psi
  1.0000 pi0 e+ e-  PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay

  End
DECAYCARD

# Peaking background: J/psi -> P gamma with the photon converting into an
# e+ e- pair in the beam pipe / inner MDC wall (dominant peaking background).
decay_card_bkg_gammaconv = <<~DECAYCARD
  Decay J/psi
  1.0000 eta' gamma  PHSP;
  Enddecay

  Decay eta'
  1.0000 gamma pi+ pi-  PHSP;
  Enddecay

  End
DECAYCARD

exMC_modeI  = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "jpsi_etap_ee_etap2gpipi"; c.related_dataset = jpsi_data
  c.events = 100_000; c.decay_card = decay_card_modeI; c.cross_section = :default
end
exMC_modeII = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "jpsi_etap_ee_etap2pipieta"; c.related_dataset = jpsi_data
  c.events = 100_000; c.decay_card = decay_card_modeII; c.cross_section = :default
end
exMC_modeIII = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "jpsi_eta_ee_eta2pipipi0"; c.related_dataset = jpsi_data
  c.events = 100_000; c.decay_card = decay_card_modeIII; c.cross_section = :default
end
exMC_modeIV = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "jpsi_eta_ee_eta2gg"; c.related_dataset = jpsi_data
  c.events = 100_000; c.decay_card = decay_card_modeIV; c.cross_section = :default
end
exMC_modeV = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "jpsi_pi0_ee_pi02gg"; c.related_dataset = jpsi_data
  c.events = 100_000; c.decay_card = decay_card_modeV; c.cross_section = :default
end
exMC_gammaconv = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "jpsi_etapgamma_gammaconv"; c.related_dataset = jpsi_data
  c.events = 500_000; c.decay_card = decay_card_bkg_gammaconv; c.cross_section = :default
end

# Common notes shared by every mode
common_notes = lambda do |alg, mode_name|
  alg.note(:photon_conversion_veto,
           "Peaking background from J/psi -> P gamma with the photon converting in the " \
           "beam pipe or the inner MDC wall is suppressed by the photon-conversion finder: " \
           "the reconstructed conversion point must satisfy delta_xy = sqrt(Rx^2 + Ry^2) " \
           "< 2 cm. This retains ~80% of the signal and removes ~98% of the conversion " \
           "events. No DSL expression exists for the conversion-point finder.")
     .note(:cos_theta_decay_cut,
           "For the modes with eta -> gamma gamma and pi0 -> gamma gamma, the decay photon " \
           "angle in the eta/pi0 helicity frame is required to satisfy |cos(theta_decay)| " \
           "< 0.9, which rejects the peaking QED continuum (e+e- -> e+e-gamma(gamma), " \
           "e+e- -> 3gamma) background. This helicity-frame quantity is not expressible in " \
           "the BOSS selection.")
     .note(:eemass_window,
           "For J/psi -> pi0 e+e- the peaking background from J/psi -> pi+ pi- pi0 is " \
           "(pi+pi- misidentified as e+e-) rejected by requiring M(e+e-) <= 0.4 GeV/c^2, " \
           "retaining about 80% of the signal.")
     .note(:raweemass_veto,
           "For J/psi -> eta' e+e- (eta' -> gamma pi+ pi-) candidates with an invariant " \
           "gamma e+e- mass in [0.10, 0.16] GeV/c^2 are vetoed to reject " \
           "J/psi -> pi+ pi- pi0 (pi0 -> gamma e+e-) background.")
     .note(:eta_pi0_mass_window,
           "In the modes with eta' -> pi+ pi- eta and eta -> pi+ pi- pi0, photon pairs are " \
           "reconstructed as eta / pi0 candidates when the gamma-gamma invariant mass lies " \
           "in (0.480, 0.600) GeV/c^2 or (0.100, 0.160) GeV/c^2 respectively; the Kalman " \
           "mass-constrained fit used here replaces that window.")
     .note(:kinematic_fit_chi2,
           "The published analysis applies chi2_4C < 100. The BOSS first pass keeps the " \
           "loose chi2_cut 200; the tight cut is applied in the ROOT analysis.")
     .note(:unbinned_ml_fit,
           "#{mode_name}: the signal yield is extracted from an unbinned extended maximum " \
           "likelihood fit of the pseudoscalar mass spectrum, with the signal MC shape " \
           "convoluted with a Gaussian and the non-peaking background described by a " \
           "Chebychev polynomial; the expected gamma-conversion peaking background is " \
           "fixed to its normalised value and subtracted.")
end

### Event selection (BOSS) ###

# ------------------------------------------------------------------ Mode I
# J/psi -> eta' e+ e-, eta' -> gamma pi+ pi-  (final state e+e-gamma pi+pi-)
alg_name = "JpsiEtapEEtogpipi"
algI = Algorithm.new(alg_name)
algI.set_header(["#{alg_name}Alg/#{alg_name}.h"])
    .set_constant({"ECMS" => [:double, 3.097]})

selI = Selection.new
selI.select_track do
      cos_theta 0.93    # |cos(theta)| < 0.93
      Vz        10.0    # within 10 cm of the IP along the beam direction
      Vr        1.0     # within 1 cm of the IP in the plane perpendicular to the beam
      nChrp     "==2"   # tracks exactly matching the final state: e+ and pi+
      nChrn     "==2"   # e- and pi-
      nNet      "==0"   # net charge zero
    end
    .select_photon do
      tdc_emc_start     0      # cluster timing requirement (noise suppression)
      tdc_emc_end       14
      angle_to_track    10.0   # photon at least 10 degrees from any charged track
      energyThreshold_b 0.025  # E > 25 MeV in the barrel (|cos(theta)| < 0.80)
      energyThreshold_e 0.050  # E > 50 MeV in the end cap (0.86 < |cos(theta)| < 0.92)
      nGam              ">=1"  # at least one photon from eta' -> gamma pi+ pi-
    end
    .pid(method: :probability) do
      prob_cut 0.001
      # Prob_PID(e) > Prob_PID(pi) and Prob_PID(e) > Prob_PID(K)
      identify :electron, against: [:pion, :kaon]
      nep "==1"                   # one e+
      nem "==1"                   # one e-
    end
    .remove([:ep <= :chrgp, :em <= :chrgn])  # electrons no longer available as pions
    .assign({chrgp: :pip, chrgn: :pim})      # remaining tracks assumed pions, no PID
    .kinematic_fit([:ep, :em, :gamma, :pip, :pim]) do
      nominal
      vertex_fit([0, 1, 3, 4])  # vertex fit on the charged tracks (loose chi2 at IP)
      constrain_four_momentum   # 4C energy-momentum conserving fit
      chi2_cut 200
    end

common_notes.call(algI, "Mode I (eta' -> gamma pi+ pi-)")
algI.with_decay_card(decay_card_modeI).apply(selI)
algI.execute_on([jpsi_data, jpsi_incMC, cont_data_3773, cont_incMC_3773, exMC_modeI, exMC_gammaconv])

# ----------------------------------------------------------------- Mode II
# J/psi -> eta' e+ e-, eta' -> pi+ pi- eta, eta -> gamma gamma
alg_name = "JpsiEtapEEtopipieta"
algII = Algorithm.new(alg_name)
algII.set_header(["#{alg_name}Alg/#{alg_name}.h"])
     .set_constant({"ECMS" => [:double, 3.097]})

selII = Selection.new
selII.select_track do
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==2"   # e+ and pi+
       nChrn     "==2"   # e- and pi-
       nNet      "==0"
     end
     .select_photon do
       tdc_emc_start     0
       tdc_emc_end       14
       angle_to_track    10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam              ">=2"   # two photons from eta -> gamma gamma
     end
     .pid(method: :probability) do
       prob_cut 0.001
       identify :electron, against: [:pion, :kaon]
       nep "==1"
       nem "==1"
     end
     .remove([:ep <= :chrgp, :em <= :chrgn])
     .assign({chrgp: :pip, chrgn: :pim})
     .kalman_kinematic_fit([:gamma, :gamma]) do
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
       chi2_cut 25
       neta ">=1"
     end
     .kinematic_fit([:ep, :em, :pip, :pim, :eta]) do
       nominal
       vertex_fit([0, 1, 2, 3])
       constrain_four_momentum
       chi2_cut 200
     end

common_notes.call(algII, "Mode II (eta' -> pi+ pi- eta, eta -> gamma gamma)")
algII.with_decay_card(decay_card_modeII).apply(selII)
algII.execute_on([jpsi_data, jpsi_incMC, cont_data_3773, cont_incMC_3773, exMC_modeII, exMC_gammaconv])

# ---------------------------------------------------------------- Mode III
# J/psi -> eta e+ e-, eta -> pi+ pi- pi0, pi0 -> gamma gamma
alg_name = "JpsiEtaEEtopipipi0"
algIII = Algorithm.new(alg_name)
algIII.set_header(["#{alg_name}Alg/#{alg_name}.h"])
      .set_constant({"ECMS" => [:double, 3.097]})

selIII = Selection.new
selIII.select_track do
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"   # e+ and pi+
        nChrn     "==2"   # e- and pi-
        nNet      "==0"
      end
      .select_photon do
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=2"   # two photons from pi0 -> gamma gamma
      end
      .pid(method: :probability) do
        prob_cut 0.001
        identify :electron, against: [:pion, :kaon]
        nep "==1"
        nem "==1"
      end
      .remove([:ep <= :chrgp, :em <= :chrgn])
      .assign({chrgp: :pip, chrgn: :pim})
      .kalman_kinematic_fit([:gamma, :gamma]) do
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      end
      .kinematic_fit([:ep, :em, :pip, :pim, :pi0]) do
        nominal
        vertex_fit([0, 1, 2, 3])
        constrain_four_momentum
        chi2_cut 200
      end

common_notes.call(algIII, "Mode III (eta -> pi+ pi- pi0)")
algIII.with_decay_card(decay_card_modeIII).apply(selIII)
algIII.execute_on([jpsi_data, jpsi_incMC, cont_data_3773, cont_incMC_3773, exMC_modeIII, exMC_gammaconv])

# ----------------------------------------------------------------- Mode IV
# J/psi -> eta e+ e-, eta -> gamma gamma   (final state e+e-gamma gamma)
alg_name = "JpsiEtaEEtogg"
algIV = Algorithm.new(alg_name)
algIV.set_header(["#{alg_name}Alg/#{alg_name}.h"])
     .set_constant({"ECMS" => [:double, 3.097]})

selIV = Selection.new
selIV.select_track do
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==1"   # e+
       nChrn     "==1"   # e-
       nNet      "==0"
     end
     .select_photon do
       tdc_emc_start     0
       tdc_emc_end       14
       angle_to_track    10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       nGam              ">=2"   # two photons from eta -> gamma gamma
     end
     .pid(method: :probability) do
       prob_cut 0.001
       identify :electron, against: [:pion, :kaon]
       nep "==1"
       nem "==1"
     end
     .kalman_kinematic_fit([:gamma, :gamma]) do
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
       chi2_cut 25
       neta ">=1"
     end
     .kinematic_fit([:ep, :em, :eta]) do
       nominal
       vertex_fit([0, 1])
       constrain_four_momentum
       chi2_cut 200
     end

common_notes.call(algIV, "Mode IV (eta -> gamma gamma)")
algIV.with_decay_card(decay_card_modeIV).apply(selIV)
algIV.execute_on([jpsi_data, jpsi_incMC, cont_data_3773, cont_incMC_3773, exMC_modeIV, exMC_gammaconv])

# ------------------------------------------------------------------ Mode V
# J/psi -> pi0 e+ e-, pi0 -> gamma gamma   (final state e+e-gamma gamma)
alg_name = "JpsiPi0EEtogg"
algV = Algorithm.new(alg_name)
algV.set_header(["#{alg_name}Alg/#{alg_name}.h"])
    .set_constant({"ECMS" => [:double, 3.097]})

selV = Selection.new
selV.select_track do
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==1"   # e+
      nChrn     "==1"   # e-
      nNet      "==0"
    end
    .select_photon do
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=2"   # two photons from pi0 -> gamma gamma
    end
    .pid(method: :probability) do
      prob_cut 0.001
      identify :electron, against: [:pion, :kaon]
      nep "==1"
      nem "==1"
    end
    .kalman_kinematic_fit([:gamma, :gamma]) do
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"
    end
    .kinematic_fit([:ep, :em, :pi0]) do
      nominal
      vertex_fit([0, 1])
      constrain_four_momentum
      chi2_cut 200
    end

common_notes.call(algV, "Mode V (pi0 -> gamma gamma)")
algV.with_decay_card(decay_card_modeV).apply(selV)
algV.execute_on([jpsi_data, jpsi_incMC, cont_data_3773, cont_incMC_3773, exMC_modeV, exMC_gammaconv])

# Continuum QED background (e+e- -> e+e-gamma(gamma), e+e- -> 3gamma) was generated
# with the Babayaga generator and studied on the 3.773 GeV data sample; it has no
# EvtGen decay-card equivalent.
algV.note(:qed_continuum_background,
          "Non-peaking continuum background from e+e- -> e+e-gamma(gamma) and " \
          "e+e- -> 3gamma (one photon converting to e+e-) is modelled with the Babayaga " \
          "QED generator and validated against the 2.9 fb^-1 psi(3773) data sample; " \
          "this generator is not available in the EvtGen decay-card framework.")
