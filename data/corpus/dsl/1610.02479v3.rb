# Amplitude analysis of psi(3686) -> gamma chi_c1, chi_c1 -> eta pi+ pi-
# using (448.0 +/- 3.1) x 10^6 psi(3686) events collected with the BESIII detector.
# Three eta decay modes are used, covering 95% of the eta decays:
#   Mode I  : eta -> gamma gamma
#   Mode II : eta -> pi+ pi- pi0
#   Mode III: eta -> pi0 pi0 pi0
# The three modes have different final-state topologies and different selection
# chains, so (Rule T1) each mode gets its own Algorithm object and Selection.

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) dataset at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # 106 x 10^6 generic psi(2S) events

# Mode I: psi(3686) -> gamma chi_c1, chi_c1 -> eta pi+ pi-, eta -> gamma gamma
decay_card_eta_gg = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c1       P2GC1;
  Enddecay

  Decay chi_c1
  1.0000 eta pi+ pi-        PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma        PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: psi(3686) -> gamma chi_c1, chi_c1 -> eta pi+ pi-, eta -> pi+ pi- pi0
decay_card_eta_pipimpi0 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c1       P2GC1;
  Enddecay

  Decay chi_c1
  1.0000 eta pi+ pi-        PHSP;
  Enddecay

  Decay eta
  1.0000 pi+ pi- pi0        PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma        PHSP;
  Enddecay

  End
DECAYCARD

# Mode III: psi(3686) -> gamma chi_c1, chi_c1 -> eta pi+ pi-, eta -> pi0 pi0 pi0
decay_card_eta_3pi0 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c1       P2GC1;
  Enddecay

  Decay chi_c1
  1.0000 eta pi+ pi-        PHSP;
  Enddecay

  Decay eta
  1.0000 pi0 pi0 pi0        PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma        PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC: one sample per eta decay mode.  Sample sizes are chosen in
# proportion to the eta branching fractions (39.41%, 22.92%, 32.68%) so that
# the complete generated set can be used directly in the MC integration.
exMC_eta_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chic1_eta_gg"
  config.related_dataset = psip_data
  config.events          = 800000
  config.decay_card      = decay_card_eta_gg
  config.cross_section   = :default
end

exMC_eta_pipimpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chic1_eta_pipimpi0"
  config.related_dataset = psip_data
  config.events          = 460000
  config.decay_card      = decay_card_eta_pipimpi0
  config.cross_section   = :default
end

exMC_eta_3pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chic1_eta_3pi0"
  config.related_dataset = psip_data
  config.events          = 660000
  config.decay_card      = decay_card_eta_3pi0
  config.cross_section   = :default
end

### Event selection (BOSS): Mode I  eta -> gamma gamma ###
# Final state: 3 photons (one radiative photon + two from the eta) and one
# pi+ pi- pair.  5C fit = 4C energy-momentum conservation + M(gamma gamma) = m_eta.
alg_name_modeI = "PsiGammaChiC1EtaToGG"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

selection_modeI = Selection.new
selection_modeI
  .select_track {
     cos_theta 0.93
     Vz        20.0
     Vr        2.0
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
     nGam              ">=3"
  }
  # All charged tracks are assumed to be pions (kaon contamination is removed by
  # the kinematic-constraint selection below).
  .assign({chrgp: :pip, chrgn: :pim})
  # eta -> gamma gamma: the two photons are mass-constrained to m_eta
  # (chi2_gamma gamma < 15 as quoted in the paper).
  .kalman_kinematic_fit([:gamma, :gamma]) {
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
     chi2_cut 15
     neta ">=1"
  }
  # Nominal 5C fit under the psi(3686) -> gamma eta pi+ pi- hypothesis: 4C
  # energy-momentum conservation (the sum of the final-state momenta is
  # constrained to the initial psi(3686) momentum) plus the M(gamma gamma) =
  # m_eta constraint from the Kalman fit above.  The combination with the
  # smallest chi2_5C is retained by the fit itself.
  .kinematic_fit([:gamma, :eta, :pip, :pim]) {
     nominal
     constrain_four_momentum
     chi2_cut 200
  }

alg_modeI
  .note(:chi2_cut_published,
        "published selection requires chi2_5C < 40; kept loose (200) in BOSS and " \
        "applied at the optimal value in ROOT.")
  .note(:radiative_photon_energy,
        "chi_c1 candidates are selected from psi(3686) -> gamma chi_c1 by requiring " \
        "the energy of the radiative photon to satisfy 0.155 < E_gamma < 0.185 GeV; " \
        "applied on the fit-updated four-momenta in ROOT.")
  .note(:track_assumption,
        "all charged tracks are assumed to be pions; the kinematic-constraint " \
        "selection removes kaon and other charged-track contamination.")
  .note(:background_veto,
        "psi(3686) -> eta J/psi background suppressed by requiring the invariant mass " \
        "of the system recoiling against the eta (with respect to the psi(3686)) to be " \
        "separated by at least 20 MeV/c^2 from the J/psi mass.")
  .note(:background_veto_pi0,
        "pi0 contamination in the eta -> gamma gamma channel rejected by discarding " \
        "events in which any two-photon combination satisfies " \
        "0.110 < m(gamma gamma) < 0.155 GeV/c^2.")
  .note(:background_veto_gammagamma_jpsi,
        "psi(3686) -> gamma gamma J/psi production suppressed by vetoing events where a " \
        "two-photon combination not forming the eta has a total energy " \
        "0.52 < E(gamma gamma) < 0.60 GeV (E(gamma gamma) ~ 0.560 GeV for the doubly " \
        "radiative psi(3686) -> gamma chi_cJ, chi_cJ -> gamma J/psi transition).")
  .note(:sideband_subtraction,
        "background from the eta sidebands, 68 < |m(gamma gamma) - m_eta| < 113 MeV/c^2, " \
        "is subtracted in the ROOT amplitude analysis.")
  .with_decay_card(decay_card_eta_gg)
  .apply(selection_modeI)

### Event selection (BOSS): Mode II  eta -> pi+ pi- pi0 ###
# Final state: 3 photons (one radiative photon + two from the pi0) and two
# pi+ pi- pairs.  5C fit = 4C energy-momentum conservation + M(gamma gamma) = m_pi0.
alg_name_modeII = "PsiGammaChiC1EtaToPiPiPi0"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

selection_modeII = Selection.new
selection_modeII
  .select_track {
     cos_theta 0.93
     Vz        20.0
     Vr        2.0
     nChrp     "==2"
     nChrn     "==2"
     nNet      "==0"
  }
  .select_photon {
     tdc_emc_start     0
     tdc_emc_end       14
     angle_to_track    10.0
     energyThreshold_b 0.025
     energyThreshold_e 0.050
     nGam              ">=3"
  }
  .assign({chrgp: :pip, chrgn: :pim})
  # pi0 -> gamma gamma mass constraint (one photon pair consumed).
  .kalman_kinematic_fit([:gamma, :gamma]) {
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
     chi2_cut 25
     npi0 ">=1"
  }
  # Nominal 5C fit under the psi(3686) -> gamma eta pi+ pi- (eta -> pi+ pi- pi0)
  # hypothesis; the combination with the smallest chi2_5C is retained.
  .kinematic_fit([:gamma, :pi0, :pip, :pip, :pim, :pim]) {
     nominal
     constrain_four_momentum
     chi2_cut 200
  }

alg_modeII
  .note(:chi2_cut_published,
        "published selection requires chi2_5C < 40; kept loose (200) in BOSS and " \
        "applied at the optimal value in ROOT.")
  .note(:eta_mass_window,
        "the eta candidates in the three-pion decay are selected by requiring the " \
        "invariant mass of the three pions to satisfy " \
        "0.535 < m(3 pi) < 0.560 GeV/c^2; applied in ROOT.")
  .note(:radiative_photon_energy,
        "chi_c1 candidates are selected from psi(3686) -> gamma chi_c1 by requiring " \
        "the energy of the radiative photon to satisfy 0.155 < E_gamma < 0.185 GeV; " \
        "applied on the fit-updated four-momenta in ROOT.")
  .note(:track_assumption,
        "all charged tracks are assumed to be pions; the kinematic-constraint " \
        "selection removes kaon and other charged-track contamination.")
  .note(:background_veto,
        "psi(3686) -> eta J/psi background suppressed by requiring the invariant mass " \
        "of the system recoiling against the eta (with respect to the psi(3686)) to be " \
        "separated by at least 20 MeV/c^2 from the J/psi mass.")
  .note(:sideband_subtraction,
        "background from the eta sidebands, 37 < |m(3 pi) - m_eta| < 62 MeV/c^2, " \
        "is subtracted in the ROOT amplitude analysis.")
  .with_decay_card(decay_card_eta_pipimpi0)
  .apply(selection_modeII)

### Event selection (BOSS): Mode III  eta -> pi0 pi0 pi0 ###
# Final state: 7 photons (one radiative photon + six from the three pi0) and one
# pi+ pi- pair.  7C fit = 4C energy-momentum conservation + three
# M(gamma gamma) = m_pi0 constraints.
alg_name_modeIII = "PsiGammaChiC1EtaTo3Pi0"
alg_modeIII = Algorithm.new(alg_name_modeIII)
alg_modeIII.set_header(["#{alg_name_modeIII}Alg/#{alg_name_modeIII}.h"])
           .set_constant({"ECMS" => [:double, 3.686]})

selection_modeIII = Selection.new
selection_modeIII
  .select_track {
     cos_theta 0.93
     Vz        20.0
     Vr        2.0
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
     nGam              ">=7"
  }
  .assign({chrgp: :pip, chrgn: :pim})
  # Three pi0 -> gamma gamma mass constraints (three photon pairs consumed).
  .kalman_kinematic_fit([:gamma, :gamma]) {
     invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
     chi2_cut 25
     npi0 ">=3"
  }
  # Nominal 7C fit under the psi(3686) -> gamma eta pi+ pi- (eta -> 3 pi0)
  # hypothesis; the combination with the smallest chi2_7C is retained.
  .kinematic_fit([:gamma, :pi0, :pi0, :pi0, :pip, :pim]) {
     nominal
     constrain_four_momentum
     chi2_cut 200
  }

alg_modeIII
  .note(:chi2_cut_published,
        "published selection requires chi2_7C < 56; kept loose (200) in BOSS and " \
        "applied at the optimal value in ROOT.")
  .note(:eta_mass_window,
        "the eta candidates in the three-pion decay are selected by requiring the " \
        "invariant mass of the three pions to satisfy " \
        "0.535 < m(3 pi) < 0.560 GeV/c^2; the same criterion as in the " \
        "eta -> pi+ pi- pi0 channel is used; applied in ROOT.")
  .note(:radiative_photon_energy,
        "chi_c1 candidates are selected from psi(3686) -> gamma chi_c1 by requiring " \
        "the energy of the radiative photon to satisfy 0.155 < E_gamma < 0.185 GeV; " \
        "applied on the fit-updated four-momenta in ROOT.")
  .note(:track_assumption,
        "all charged tracks are assumed to be pions; the kinematic-constraint " \
        "selection removes kaon and other charged-track contamination.")
  .note(:background_veto,
        "psi(3686) -> eta J/psi background suppressed by requiring the invariant mass " \
        "of the system recoiling against the eta (with respect to the psi(3686)) to be " \
        "separated by at least 20 MeV/c^2 from the J/psi mass.")
  .note(:sideband_subtraction,
        "background from the eta sidebands, 37 < |m(3 pi) - m_eta| < 62 MeV/c^2, " \
        "is subtracted in the ROOT amplitude analysis.")
  .with_decay_card(decay_card_eta_3pi0)
  .apply(selection_modeIII)

### Execute the algorithms on data, inclusive MC and the per-mode exclusive MC ###
root_files_modeI   = alg_modeI.execute_on([psip_data, psip_incMC, exMC_eta_gg])
root_files_modeII  = alg_modeII.execute_on([psip_data, psip_incMC, exMC_eta_pipimpi0])
root_files_modeIII = alg_modeIII.execute_on([psip_data, psip_incMC, exMC_eta_3pi0])
