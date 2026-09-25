# ============================================================================
# arXiv:1210.3746v1
# "Precision measurements of branching fractions for psi' -> pi0 J/psi and
#  eta J/psi"
# BESIII, 106.4e6 psi(2S) events at sqrt(s) = 3.686 GeV
#
# Signal chain : psi' -> pi0(eta) J/psi, pi0(eta) -> gamma gamma,
#                J/psi -> l+ l-  (l = e, mu)
# The pi0 and eta channels share the same final state and the same selection
# chain; two Algorithms are used because each carries its own decay card.
# The J/psi -> e+e- and mu+mu- sub-channels are covered by the same lepton PID.
# ============================================================================

### ---------------------------- Dataset preparation ---------------------------- ###
psip_data    = DatasetManager.real_data.find("709_3686")      # psi(2S) real data
psip_incMC   = DatasetManager.inclusive_mc.find("709_3686")   # Inclusive psi(2S) MC
cont_data    = DatasetManager.real_data.find("709_3650")      # 3.65 GeV continuum for QED background

# --- psi' -> pi0 J/psi ------------------------------------------------------
decay_card_pi0_ee = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi0 J/psi PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

decay_card_pi0_mumu = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi0 J/psi PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  1.0000 mu+ mu- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

exMC_pi0_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip2pi0Jpsi_ee"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_pi0_ee
  config.cross_section   = :default
end

exMC_pi0_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip2pi0Jpsi_mumu"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_pi0_mumu
  config.cross_section   = :default
end

# --- psi' -> eta J/psi ------------------------------------------------------
decay_card_eta_ee = <<~DECAYCARD
  Decay psi(2S)
  1.0000 eta J/psi PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

decay_card_eta_mumu = <<~DECAYCARD
  Decay psi(2S)
  1.0000 eta J/psi PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  1.0000 mu+ mu- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

exMC_eta_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip2etaJpsi_ee"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_eta_ee
  config.cross_section   = :default
end

exMC_eta_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip2etaJpsi_mumu"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_eta_mumu
  config.cross_section   = :default
end

### ============================================================================
### Mode I event selection : psi' -> pi0 J/psi -> gamma gamma l+ l-
### ============================================================================
alg_pi0 = Algorithm.new("Psippipi0Jpsi")
alg_pi0.set_header(["Psippipi0JpsiAlg/Psippipi0Jpsi.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

sel_pi0 = Selection.new
sel_pi0.select_track {                    # exactly two charged tracks, net charge zero
          cos_theta 0.93                  # |cos(theta)| < 0.93
          Vz        10.0                  # |Vz| < 10 cm along the beam direction
          Vr        1.0                   # Vr < 1 cm in the transverse plane
          nChrp     "==1"
          nChrn     "==1"
          nNet      "==0"
        }
       .select_photon {                   # two photons from the pi0
          tdc_emc_start     0
          tdc_emc_end       14            # EMC-cluster timing window
          angle_to_track    10.0          # photon >= 10 deg from the nearest charged track
          energyThreshold_b 0.025         # barrel  E > 25 MeV
          energyThreshold_e 0.050         # endcap  E > 50 MeV
          nGam              ">=2"
        }
       .pid(method: :probability) {       # J/psi -> l+ l- : lepton (e/mu) identification
          identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                         treat_as_electron_if_energy_above: 0.6
          nlp "==1"                       # one lepton of each charge gives the J/psi
          nlm "==1"
        }
       # 4C kinematic fit with the hypothesis psi' -> gamma gamma l+ l- constrained to
       # the initial e+ e- beam four-momentum. The paper applies chi2 < 100 (60) for the
       # pi0 e+e- (pi0 mu+mu-) final states; the tight value is optimised on S/sqrt(S+B)
       # in ROOT, so only the loose default cut is applied here.
       .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
          nominal
          constrain_four_momentum
          chi2_cut 200
        }

alg_pi0
  .note(:lep_pid_method, "leptons are identified by the E/p ratio rather than the
    standard dE/dx-TOF probability: muons require 0.08 < E/p < 0.22 and electrons
    require E/p > 0.8 (E = EMC deposited energy, p = track momentum)")
  .note(:mass_window, "M(l+l-) within (3.05, 3.15) GeV/c^2 (J/psi mass window) is
    applied in ROOT after the 4C kinematic fit")
  .note(:background_veto, "chi_cJ background suppressed by requiring
    M(gamma_h J/psi) outside (3.50, 3.57) GeV/c^2, where gamma_h is the more
    energetic photon; applied in ROOT after the nominal kinematic fit")
  .note(:efficiency_curve, "the angular distributions (1 + cos^2(theta) for the J/psi
    and for the leptons in the J/psi rest frame) are modelled in the event generator
    to determine the detection efficiency; efficiencies are 23.05% (gamma gamma e+e-)
    and 29.11% (gamma gamma mu+mu-)")
  .with_decay_card(decay_card_pi0_ee)
  .apply(sel_pi0)

### ============================================================================
### Mode II event selection : psi' -> eta J/psi -> gamma gamma l+ l-
### ============================================================================
alg_eta = Algorithm.new("PsipetaJpsi")
alg_eta.set_header(["PsipetaJpsiAlg/PsipetaJpsi.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

sel_eta = Selection.new
sel_eta.select_track {
          cos_theta 0.93
          Vz        10.0
          Vr        1.0
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
          nGam              ">=2"
        }
       .pid(method: :probability) {
          identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                         treat_as_electron_if_energy_above: 0.6
          nlp "==1"
          nlm "==1"
        }
       # 4C fit, psi' -> gamma gamma l+ l-. Paper cut: chi2 < 70 (50) for the
       # eta e+e- (eta mu+mu-) final states, optimised in ROOT and applied there.
       .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
          nominal
          constrain_four_momentum
          chi2_cut 200
        }

alg_eta
  .note(:lep_pid_method, "leptons are identified by the E/p ratio rather than the
    standard dE/dx-TOF probability: muons require 0.08 < E/p < 0.22 and electrons
    require E/p > 0.8 (E = EMC deposited energy, p = track momentum)")
  .note(:mass_window, "M(l+l-) within (3.05, 3.15) GeV/c^2 (J/psi mass window) is
    applied in ROOT after the 4C kinematic fit")
  .note(:background_veto, "chi_cJ background suppressed by requiring
    M(gamma_h J/psi) < 3.5 GeV/c^2, where gamma_h is the more energetic photon;
    applied in ROOT after the nominal kinematic fit. An additional smooth
    background from psi' -> pi0 pi0 J/psi lies inside the eta mass region")
  .note(:efficiency_curve, "the angular distributions (1 + cos^2(theta) for the J/psi
    and for the leptons in the J/psi rest frame) are modelled in the event generator
    to determine the detection efficiency; efficiencies are 35.41% (gamma gamma e+e-)
    and 46.28% (gamma gamma mu+mu-)")
  .with_decay_card(decay_card_eta_ee)
  .apply(sel_eta)

### --------------------------------- Execution --------------------------------- ###
# The QED continuum sample at 3.65 GeV is used for background studies; it is
# analysed with the same selection and therefore shares the algorithms' jobs.
alg_pi0.execute_on([psip_data, psip_incMC, cont_data,
                    exMC_pi0_ee, exMC_pi0_mumu])

alg_eta.execute_on([psip_data, psip_incMC, cont_data,
                    exMC_eta_ee, exMC_eta_mumu])
