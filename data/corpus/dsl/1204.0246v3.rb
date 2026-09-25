# =============================================================================
# BESIII: First measurement of the two-photon transition
#         psi(3686) -> gamma gamma J/psi , J/psi -> l+ l-  (l = e, mu)
# arXiv:1204.0246v3
# 106 million psi(3686) decays; 156.4 pb^-1 at sqrt(s) = 3.686 GeV
# Continuum data (42.6 pb^-1 at 3.65 GeV) used for background evaluation
# =============================================================================

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# ---- Decay cards -------------------------------------------------------------
# Signal: two-photon transition psi(3686) -> gamma gamma J/psi, J/psi -> e+ e-
decay_card_ee = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma gamma J/psi   PHSP;
  Enddecay

  Decay J/psi
  1.000 e+ e-               VLL;
  Enddecay

  End
DECAYCARD

# Signal: psi(3686) -> gamma gamma J/psi, J/psi -> mu+ mu-
decay_card_mumu = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma gamma J/psi   PHSP;
  Enddecay

  Decay J/psi
  1.000 mu+ mu-             VLL;
  Enddecay

  End
DECAYCARD

# Dominant background / also measured: cascade E1 transitions
# psi(3686) -> gamma chi_cJ, chi_cJ -> gamma J/psi (J = 0, 1, 2)
# The chi_cJ line shapes are simulated with E_gamma1^3 E_gamma2^3 weights to
# account for the double E1 transitions.
decay_card_cascade = <<~DECAYCARD
  Decay psi(2S)
  1.000 gamma chi_c0        PHSP;
  Enddecay

  Decay chi_c0
  1.000 gamma J/psi         PHSP;
  Enddecay

  Decay J/psi
  1.000 e+ e-               VLL;
  Enddecay

  Decay J/psi
  1.000 mu+ mu-             VLL;
  Enddecay

  End
DECAYCARD

# ---- Exclusive MC samples ----------------------------------------------------
exMC_ee = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psip_gammagamma_Jpsi_ee"
  c.related_dataset = psip_data
  c.events          = 200000
  c.decay_card      = decay_card_ee
  c.cross_section   = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psip_gammagamma_Jpsi_mumu"
  c.related_dataset = psip_data
  c.events          = 200000
  c.decay_card      = decay_card_mumu
  c.cross_section   = :default
end

exMC_cascade = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psip_gamma_chicJ_gamma_Jpsi"
  c.related_dataset = psip_data
  c.events          = 500000
  c.decay_card      = decay_card_cascade
  c.cross_section   = :default
end

# =============================================================================
# ALGORITHM 1: psi(3686) -> gamma gamma J/psi, J/psi -> e+ e-
# =============================================================================
alg_ee = Algorithm.new("PsipGammaGammaJpsiEE")
alg_ee.set_header(["PsipGammaGammaJpsiEEAlg/PsipGammaGammaJpsiEE.h"])
      .set_constant({"ECMS" => [:double, 3.686]})   # psi(3686) center-of-mass energy (GeV)

sel_ee = Selection.new
sel_ee.select_track {
         cos_theta  0.93   # |cos(theta)| < 0.93 with respect to the beam direction
         Vz         10.0   # track-vertex distance < 10 cm along the beam axis
         Vr         1.0    # track-vertex distance < 1 cm in the transverse plane
         nChrp      "==1"  # exactly two oppositely charged good tracks (the dilepton)
         nChrn      "==1"
         nNet       "==0"
       }
       .select_photon {
         # unmatched EMC showers: > 25 MeV barrel (|cos|<0.8), > 50 MeV endcap
         # (0.86<|cos|<0.92); 700 ns timing window suppresses noise
         tdc_emc_start     0
         tdc_emc_end       14
         angle_to_track    10.0   # rejects bremsstrahlung photons from either lepton
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         nGam              "==2"  # events are required to have only two photon candidates
       }
       .for_each(:charged) {
         # suppress non-J/psi decay leptons: each lepton must have p > 0.8 GeV/c
         where { p < 0.8 }
         remove
       }
       # Lepton identification with the EMC shower energy over MDC track momentum ratio
       # (E/p > 0.7 for an electron)
       .pid(method: :probability) {
         prob_cut 0.001
         identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                        treat_as_electron_if_energy_above: 0.6
         nlp "==1"
         nlm "==1"
       }
       # Kinematic fit: the vertexed dilepton is constrained to the nominal J/psi mass
       # (1C) and the resulting J/psi together with the two photon candidates to the
       # known initial psi(3686) four-momentum (4C). A vertex fit constrains the
       # production vertex and the dilepton tracks to a common vertex.
       .kinematic_fit([:ep, :em, :gamma, :gamma]) {
         nominal
         vertex_fit([0, 1])   # common vertex for the two lepton tracks
         invariant_mass_of(:ep, :em).constrain_to_nominal_mass_of(:jpsi)
         constrain_four_momentum
         chi2_cut 200        # chi2_KF/ndof < 12 applied at the ROOT level
       }

alg_ee
  .note(:vertex_fit_chi2,
        "the vertex fit (VF) constrains the run-by-run production vertex and the dilepton tracks to a " \
        "common vertex; only events with chi2_VF/ndof < 20 are accepted. The vertex constraint is applied " \
        "inside the nominal kinematic fit via vertex_fit([0,1]); the chi2/ndof quality requirement itself " \
        "is not expressible in the DSL.")
  .note(:kf_chi2,
        "the kinematic-fit quality chi2_KF/ndof < 12 is required; the corresponding chi2 value is stored " \
        "and the cut is applied at the ROOT level.")
  .note(:lepton_id_ep,
        "the electron is identified with the ratio of EMC shower energy to MDC track momentum, E/p > 0.7, " \
        "in addition to the high-momentum lepton PID; the E/p ratio cut is applied at the ROOT level.")
  .note(:jpsi_recoil_window,
        "J/psi candidates are identified by requiring the recoil mass of the two photons, " \
        "M_{gamma gamma - recoil}, to lie within (3.08, 3.14) GeV/c^2; this uses the kinematic-fit " \
        "corrected four-momenta and is applied at the ROOT level.")
  .note(:background_veto,
        "continuum background (Bhabha scattering, dimuon production and ISR J/psi) is excluded by " \
        "discarding events with M_{gamma_sm - recoil} > 3.6 GeV/c^2. Backgrounds from " \
        "psi(3686) -> pi0(eta) J/psi are suppressed by requiring the diphoton invariant mass " \
        "M_{gamma gamma} > 0.15 GeV/c^2 and the diphoton recoil momentum > 0.25 GeV/c.")
  .note(:spin_structure,
        "the signal MC momenta are generated according to the measured polarization structure of the " \
        "two-photon transition; the input spin structure contributes a 20% systematic uncertainty.")
  .note(:continuum_background,
        "non-psi(3686) backgrounds are estimated from the 42.6 pb^-1 continuum data sample at " \
        "sqrt(s) = 3.65 GeV, scaled by luminosity and the 1/s dependence of the cross sections.")
  .note(:signal_extraction,
        "the signal yield is extracted from a global unbinned maximum-likelihood fit to the " \
        "M_{gamma_sm - recoil} spectrum with the psi(3686) backgrounds fixed from MC simulation and " \
        "the signal shape smeared with an asymmetric Gaussian (ROOT level).")

alg_ee.with_decay_card(decay_card_ee).apply(sel_ee)
root_files_ee = alg_ee.execute_on([psip_data, psip_incMC, exMC_ee, exMC_cascade])

# =============================================================================
# ALGORITHM 2: psi(3686) -> gamma gamma J/psi, J/psi -> mu+ mu-
# =============================================================================
alg_mumu = Algorithm.new("PsipGammaGammaJpsiMuMu")
alg_mumu.set_header(["PsipGammaGammaJpsiMuMuAlg/PsipGammaGammaJpsiMuMu.h"])
        .set_constant({"ECMS" => [:double, 3.686]})

sel_mumu = Selection.new
sel_mumu.select_track {
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
           nGam              "==2"
         }
         .for_each(:charged) {
           where { p < 0.8 }
           remove
         }
         # E/p < 0.6 for a muon
         .pid(method: :probability) {
           prob_cut 0.001
           identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                          treat_as_electron_if_energy_above: 0.6
           nlp "==1"
           nlm "==1"
         }
         .kinematic_fit([:mup, :mum, :gamma, :gamma]) {
           nominal
           vertex_fit([0, 1])
           invariant_mass_of(:mup, :mum).constrain_to_nominal_mass_of(:jpsi)
           constrain_four_momentum
           chi2_cut 200
         }

alg_mumu
  .note(:vertex_fit_chi2,
        "the vertex fit constrains the production vertex and the dimuon tracks to a common vertex; only " \
        "events with chi2_VF/ndof < 20 are accepted (applied at the ROOT level).")
  .note(:kf_chi2,
        "chi2_KF/ndof < 12 is required; the chi2 is stored and the cut is applied at the ROOT level.")
  .note(:lepton_id_ep,
        "the muon is identified with E/p < 0.6, in addition to the high-momentum lepton PID and the " \
        "MUC information; the E/p cut is applied at the ROOT level.")
  .note(:jpsi_recoil_window,
        "M_{gamma gamma - recoil} within (3.08, 3.14) GeV/c^2 identifies the J/psi candidate (ROOT level).")
  .note(:background_veto,
        "events with M_{gamma_sm - recoil} > 3.6 GeV/c^2 are discarded to remove continuum background; " \
        "M_{gamma gamma} > 0.15 GeV/c^2 and the diphoton recoil momentum > 0.25 GeV/c suppress " \
        "psi(3686) -> pi0(eta) J/psi backgrounds.")
  .note(:spin_structure,
        "signal MC generated according to the measured polarization structure; 20% systematic uncertainty.")
  .note(:continuum_background,
        "continuum background estimated from the 42.6 pb^-1 sample at sqrt(s) = 3.65 GeV, scaled by " \
        "luminosity and the 1/s dependence of the cross sections.")
  .note(:signal_extraction,
        "global unbinned maximum-likelihood fit to M_{gamma_sm - recoil} (ROOT level).")

alg_mumu.with_decay_card(decay_card_mumu).apply(sel_mumu)
root_files_mumu = alg_mumu.execute_on([psip_data, psip_incMC, exMC_mumu, exMC_cascade])
