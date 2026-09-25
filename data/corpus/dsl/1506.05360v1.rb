# Analysis: Measurement of matrix elements for eta -> pi+pi-pi0 and eta/eta' -> pi0 pi0 pi0
# Dataset: 1.31e9 J/psi events at BESIII

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# --------------------------------------------------------------
# Decay cards
# --------------------------------------------------------------
decay_card_eta_pipipi0 = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta                        PHSP;
  Enddecay

  Decay eta
  1.000 pi+ pi- pi0                      ETA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma                      PHSP;
  Enddecay

  End
DECAYCARD

decay_card_eta_3pi0 = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta                        PHSP;
  Enddecay

  Decay eta
  1.000 pi0 pi0 pi0                      PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma                      PHSP;
  Enddecay

  End
DECAYCARD

decay_card_etap_3pi0 = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta'                       PHSP;
  Enddecay

  Decay eta'
  1.000 pi0 pi0 pi0                      PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma                      PHSP;
  Enddecay

  End
DECAYCARD

# --------------------------------------------------------------
# Exclusive MC samples
# --------------------------------------------------------------
exMC_eta_pipipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_eta_pipipi0"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_eta_pipipi0
  config.cross_section   = :default
end

exMC_eta_3pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_eta_3pi0"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_eta_3pi0
  config.cross_section   = :default
end

exMC_etap_3pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etap_3pi0"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_etap_3pi0
  config.cross_section   = :default
end

# ==============================================================
# Algorithm A: J/psi -> gamma eta, eta -> pi+ pi- pi0
# Topology: pi+ pi- gamma gamma gamma
# ==============================================================
alg_A_name = "EtaToPiPiPi0"
alg_A = Algorithm.new(alg_A_name)
alg_A.set_header(["#{alg_A_name}Alg/#{alg_A_name}.h"])
     .set_constant({"ECMS" => [:double, 3.097]})

sel_A = Selection.new
sel_A.select_track {
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
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=3"
      }
      # Nominal 6C kinematic fit:
      # E-p conservation + M(gamma gamma) = m(pi0) + M(pi+ pi- pi0) = m(eta)
      # (The most energetic photon is treated as the radiative photon; the two
      # remaining photons form the pi0.)
      .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        invariant_mass_of(:gamma, :gamma, :gamma, :pip, :pim).constrain_to_nominal_mass_of(:eta)
        chi2_cut 80
      }
      # Competing hypothesis: J/psi -> pi+ pi- gamma gamma gamma gamma (4C)
      .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim]) {
        constrain_four_momentum
      }
      # Competing hypothesis: J/psi -> pi+ pi- gamma gamma (4C)
      .kinematic_fit([:gamma, :gamma, :pip, :pim]) {
        constrain_four_momentum
      }

alg_A.note(:radiative_photon,
           "The most energetic photon in the event is treated as the radiative photon " \
           "from the J/psi decay; the two remaining photons are combined into the pi0.")
     .note(:background_veto,
           "Compare chi2 of signal 4C hypothesis (J/psi -> pi+ pi- 3 gamma) with 4C chi2 " \
           "of pi+ pi- 4 gamma and pi+ pi- 2 gamma background hypotheses; discard event " \
           "when the signal chi2 exceeds either background chi2.")

alg_A.with_decay_card(decay_card_eta_pipipi0).apply(sel_A)
alg_A.execute_on([jpsi_data, jpsi_incMC, exMC_eta_pipipi0])

# ==============================================================
# Algorithm B: J/psi -> gamma eta, eta -> pi0 pi0 pi0
# Topology: 7 photons, no charged tracks
# ==============================================================
alg_B_name = "EtaTo3Pi0"
alg_B = Algorithm.new(alg_B_name)
alg_B.set_header(["#{alg_B_name}Alg/#{alg_B_name}.h"])
     .set_constant({"ECMS" => [:double, 3.097]})

sel_B = Selection.new
sel_B.select_track {
        cos_theta  0.93
        Vz         10.0
        Vr         1.0
        nChrp      "==0"
        nChrn      "==0"
        nNet       "==0"
      }
      .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=7"
      }
      # Reconstruct pi0 from photon pairs (1C Kalman fit)
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0     ">=3"
      }
      # Nominal 7C kinematic fit: J/psi -> gamma pi0 pi0 pi0 with
      # energy-momentum conservation + 3 pi0 mass constraints,
      # plus additional eta mass constraint on the 3 pi0 system.
      .kinematic_fit([:gamma, :pi0, :pi0, :pi0]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:pi0, :pi0, :pi0).constrain_to_nominal_mass_of(:eta)
        chi2_cut 70
      }

alg_B.note(:radiative_photon,
           "The most energetic photon is taken as the radiative photon from the J/psi decay; " \
           "the remaining photons are paired into three pi0 candidates.")
     .note(:pi0_decay_angle,
           "|cos(theta_decay)| < 0.95, where theta_decay is the polar angle of a photon in the " \
           "gamma-gamma rest frame, to suppress pi0 mis-combination.")
     .note(:pi0_pi0_pi0_mass_constraint,
           "An additional kinematic fit constraining M(pi0 pi0 pi0) to the eta nominal mass " \
           "is applied to improve the Dalitz-plot Z-variable resolution.")

alg_B.with_decay_card(decay_card_eta_3pi0).apply(sel_B)
alg_B.execute_on([jpsi_data, jpsi_incMC, exMC_eta_3pi0])

# ==============================================================
# Algorithm C: J/psi -> gamma eta', eta' -> pi0 pi0 pi0
# Topology: 7 photons, no charged tracks
# ==============================================================
alg_C_name = "EtaPrimeTo3Pi0"
alg_C = Algorithm.new(alg_C_name)
alg_C.set_header(["#{alg_C_name}Alg/#{alg_C_name}.h"])
     .set_constant({"ECMS" => [:double, 3.097]})

sel_C = Selection.new
sel_C.select_track {
        cos_theta  0.93
        Vz         10.0
        Vr         1.0
        nChrp      "==0"
        nChrn      "==0"
        nNet       "==0"
      }
      .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=7"
      }
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0     ">=3"
      }
      # Nominal 7C kinematic fit: J/psi -> gamma pi0 pi0 pi0 with
      # E-p conservation + 3 pi0 mass constraints, plus M(3 pi0) = m(eta').
      .kinematic_fit([:gamma, :pi0, :pi0, :pi0]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:pi0, :pi0, :pi0).constrain_to_nominal_mass_of(:etap)
        chi2_cut 70
      }
      # Competing hypothesis: J/psi -> gamma eta pi0 pi0 (7C fit) — used to
      # veto eta' -> eta pi0 pi0 background.
      .kinematic_fit([:gamma, :pi0, :pi0, :gamma, :gamma]) {
        constrain_four_momentum
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      }

alg_C.note(:radiative_photon,
           "The most energetic photon is taken as the radiative photon from the J/psi decay.")
     .note(:pi0_decay_angle,
           "|cos(theta_decay)| < 0.95 for each pi0 candidate to suppress pi0 mis-combination.")
     .note(:omega_veto,
           "Veto J/psi -> omega pi0 pi0 background by requiring |M(gamma pi0) - m_omega| >= 0.05 GeV/c^2 " \
           "for all gamma-pi0 combinations (gamma is the radiative photon).")
     .note(:eta_pair_veto,
           "Reject events in which any gamma gamma pair has |M(gg) - m_eta| < 0.03 GeV/c^2.")
     .note(:eta_pi0_pi0_veto,
           "Alternative 7C fit under J/psi -> gamma eta pi0 pi0 hypothesis; if its chi2 is smaller " \
           "than the signal-hypothesis chi2, the event is discarded.")
     .note(:pi0_pi0_pi0_mass_constraint,
           "Additional kinematic constraint M(pi0 pi0 pi0) = m(eta') is applied to improve the " \
           "resolution of the Dalitz-plot Z variable.")

alg_C.with_decay_card(decay_card_etap_3pi0).apply(sel_C)
alg_C.execute_on([jpsi_data, jpsi_incMC, exMC_etap_3pi0])
