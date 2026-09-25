# BESIII arXiv:1412.5258v1
# Study of J/psi -> eta phi pi+ pi- with eta -> gamma gamma and phi -> K+ K-
# (2.25 x 10^8 J/psi events at sqrt(s) = 3.097 GeV)
# Signal chain: J/psi -> eta Y(2175), Y(2175) -> phi f_0(980), f_0(980) -> pi+ pi-
# Also studied in the same final state: J/psi -> phi f_1(1285) / phi eta(1405) / phi X(1835) / phi X(1870)

### Dataset ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi data sample
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # inclusive J/psi MC (background estimate)

### Decay cards — all share the same final state K+ K- gamma gamma pi+ pi- ###
# Nominal signal: J/psi -> eta Y(2175), Y(2175) -> phi f_0(980), f_0(980) -> pi+ pi-
# l(eta-Y(2175)) = 1 (P-wave), l(phi-f_0(980)) = 0 (S-wave), l(pi+pi-) = 0 (S-wave);
# the f_0(980) shape follows the Flatte parameterisation (BESII parameters).
decay_card_y2175 = <<~DECAYCARD
  Decay J/psi
  1.000 eta Y(2175)                               PHSP;
  Enddecay

  Decay Y(2175)
  1.000 phi f_0(980)                              PHSP;
  Enddecay

  Decay f_0(980)
  1.000 pi+ pi-                                   PHSP;
  Enddecay

  Decay phi
  1.000 K+ K-                                     VSS;
  Enddecay

  Decay eta
  1.000 gamma gamma                               PHSP;
  Enddecay

  End
DECAYCARD

# J/psi -> phi f_1(1285), f_1(1285) -> eta pi+ pi- (same final state)
decay_card_f11285 = <<~DECAYCARD
  Decay J/psi
  1.000 phi f_1(1285)                             PHSP;
  Enddecay

  Decay f_1(1285)
  1.000 eta pi+ pi-                               PHSP;
  Enddecay

  Decay phi
  1.000 K+ K-                                     VSS;
  Enddecay

  Decay eta
  1.000 gamma gamma                               PHSP;
  Enddecay

  End
DECAYCARD

# J/psi -> phi eta(1405), eta(1405) -> eta pi+ pi- (same final state)
decay_card_eta1405 = <<~DECAYCARD
  Decay J/psi
  1.000 phi eta(1405)                             PHSP;
  Enddecay

  Decay eta(1405)
  1.000 eta pi+ pi-                               PHSP;
  Enddecay

  Decay phi
  1.000 K+ K-                                     VSS;
  Enddecay

  Decay eta
  1.000 gamma gamma                               PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive signal MC ###
exMC_y2175 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_etaY2175_phi_f0"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_y2175
  config.cross_section   = :default
end

exMC_f11285 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_phi_f1_1285"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_f11285
  config.cross_section   = :default
end

exMC_eta1405 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_phi_eta_1405"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_eta1405
  config.cross_section   = :default
end

### Algorithm: J/psi -> eta phi pi+ pi- (common selection for all resonances in this final state) ###
alg = Algorithm.new("JpsiEtaPhiPiPi")
alg.set_header(["JpsiEtaPhiPiPiAlg/JpsiEtaPhiPiPi.h"])
   .set_constant({"ECMS" => [:double, 3.097]})

selection = Selection.new
selection.select_track {
        cos_theta 0.93      # |cos(theta)| < 0.93 in the MDC
        Vz        20.0      # within +/-20 cm of the IP along the beam direction
        Vr        2.0       # within 2 cm in the plane perpendicular to the beam
        nChrp     "==2"     # K+ and pi+ (opposite charges required)
        nChrn     "==2"     # K- and pi-
        nNet      "==0"     # net charge zero
      }
      .select_photon {
        angle_to_track    10.0   # > 10 degrees from the nearest charged track (bremsstrahlung veto)
        energyThreshold_b 0.025  # > 25 MeV for barrel showers (|cos(theta)| < 0.80)
        energyThreshold_e 0.050  # > 50 MeV for end-cap showers (0.86 < |cos(theta)| < 0.92)
        tdc_emc_start     0      # EMC cluster timing, 0-700 ns
        tdc_emc_end       14
        nGam              ">=2"  # eta -> gamma gamma
      }
      .pid(method: :probability) {
        prob_cut 0.001
        # Each track is assigned to the hypothesis (pi, K, p) with the highest PID
        # confidence level; two kaons and two pions with opposite charges are required
        identify :kaon, against: [:pion, :proton]
        identify :pion, against: [:kaon, :proton]
        nkp  "==1"
        nkm  "==1"
        npip "==1"
        npim "==1"
      }
      # 4C kinematic fit to the J/psi -> K+ K- pi+ pi- gamma gamma hypothesis;
      # all photon-pair combinations are tried and the smallest chi2 is retained
      .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
        invariant_mass_of(:gamma, :gamma).within(0.528, 0.566)  # eta: |M(gg) - M_eta| < 0.019
        invariant_mass_of(:kp, :km).within(1.006, 1.032)        # phi: |M(KK) - M_phi| < 0.013
      }

alg.note(:pid_assignment,
  "The PID combines TOF and dE/dx confidence levels for the pi / K / p hypotheses and " \
  "assigns each track to the hypothesis with the highest confidence level, rather than " \
  "applying a fixed probability threshold.")
    .note(:kinematic_fit_combination,
  "The 4C fit is performed on the J/psi -> K+K- pi+pi- gamma gamma hypothesis; all " \
  "combinations of two photons are tried and the one with the smallest chi2_4C is retained.")
    .note(:mass_windows,
  "After the kinematic fit, the eta and phi signal regions are defined as " \
  "|M(gamma gamma) - M_eta| < 0.019 GeV/c^2 and |M(K+ K-) - M_phi| < 0.013 GeV/c^2 " \
  "(world-average eta and phi masses). The corresponding sideband regions are " \
  "0.480 < M(gamma gamma) < 0.499 GeV/c^2 or 0.577 < M(gamma gamma) < 0.596 GeV/c^2 for " \
  "the eta and 1.070 < M(K+ K-) < 1.096 GeV/c^2 for the phi; the f_0(980) window is " \
  "0.90 < M(pi+ pi-) < 1.05 GeV/c^2. These regions are used in the ROOT-level fit.")
    .note(:f0_980_lineshape,
  "The f_0(980) shape is parameterised with the Flatte formula using the parameters " \
  "measured by BESII; the alternative parameterisation of Zou and Bugg gives a 7.6% " \
  "difference in detection efficiency, which is taken as a systematic uncertainty.")
    .note(:signal_model,
  "The Y(2175) is generated with a P-wave eta-Y(2175) system and S-wave phi-f_0(980) and " \
  "pi+ pi- systems. For the J/psi -> phi f_1(1285) and phi eta(1405) signal MC samples the " \
  "angular distributions are also taken into account.")

alg.with_decay_card(decay_card_y2175).apply(selection)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_y2175, exMC_f11285, exMC_eta1405])
