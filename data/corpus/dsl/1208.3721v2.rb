# =============================================================================
# BESIII paper arXiv:1208.3721v2
# Measurement of the branching fractions of chi_cJ -> p nbar pi- and
# chi_cJ -> p nbar pi- pi0 (J = 0, 1, 2), produced via psi(2S) -> gamma chi_cJ.
#
# Data: (1.06 +/- 0.04) x 10^8 psi(2S) events (156.4 pb^-1) at sqrt(s) = 3.686 GeV
#       + 42.6 pb^-1 continuum data at sqrt(s) = 3.65 GeV.
#
# The two final states (p nbar pi- and p nbar pi- pi0) require different photon
# multiplicities and kinematic-fit hypotheses, and the two charge-conjugate
# channels are measured separately in the paper (left/right columns of
# Tables I and II) with different PID confidence-level requirements; therefore
# four Algorithms are built:
#   psi(2S) -> gamma chi_cJ, chi_cJ -> p nbar pi-        (1 photon, 1 p, 1 pi-)
#   psi(2S) -> gamma chi_cJ, chi_cJ -> pbar n pi+        (charge conjugate)
#   psi(2S) -> gamma chi_cJ, chi_cJ -> p nbar pi- pi0    (3 photons)
#   psi(2S) -> gamma chi_cJ, chi_cJ -> pbar n pi+ pi0    (charge conjugate)
# Within one algorithm the chi_c0, chi_c1 and chi_c2 decays share the same final
# state and the same selection (only the decay-card mother differs).
# =============================================================================

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(2S) data, 1.06e8 events
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # psi(2S) inclusive MC
cont_data  = DatasetManager.real_data.find("709_3650")     # continuum data at 3.65 GeV
cont_incMC = DatasetManager.inclusive_mc.find("709_3650")  # continuum inclusive MC

# -----------------------------------------------------------------------------
# Decay cards (EvtGen format)
# -----------------------------------------------------------------------------
# Signal: psi(2S) -> gamma chi_cJ, chi_cJ -> p nbar pi-      (J = 0, 1, 2)
decay_card_pnbarpi_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0          PHSP;
    Enddecay

    Decay chi_c0
    1.0000 p+ anti-n0 pi-        PHSP;
    Enddecay

    End
DECAYCARD

decay_card_pnbarpi_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1          PHSP;
    Enddecay

    Decay chi_c1
    1.0000 p+ anti-n0 pi-        PHSP;
    Enddecay

    End
DECAYCARD

decay_card_pnbarpi_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2          PHSP;
    Enddecay

    Decay chi_c2
    1.0000 p+ anti-n0 pi-        PHSP;
    Enddecay

    End
DECAYCARD

# Signal (charge conjugate): psi(2S) -> gamma chi_cJ, chi_cJ -> pbar n pi+
decay_card_pbarnpi_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0          PHSP;
    Enddecay

    Decay chi_c0
    1.0000 anti-p- n0 pi+        PHSP;
    Enddecay

    End
DECAYCARD

decay_card_pbarnpi_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1          PHSP;
    Enddecay

    Decay chi_c1
    1.0000 anti-p- n0 pi+        PHSP;
    Enddecay

    End
DECAYCARD

decay_card_pbarnpi_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2          PHSP;
    Enddecay

    Decay chi_c2
    1.0000 anti-p- n0 pi+        PHSP;
    Enddecay

    End
DECAYCARD

# Signal: psi(2S) -> gamma chi_cJ, chi_cJ -> p nbar pi- pi0, pi0 -> gamma gamma
decay_card_pnbarpipi0_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0          PHSP;
    Enddecay

    Decay chi_c0
    1.0000 p+ anti-n0 pi- pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma           PHSP;
    Enddecay

    End
DECAYCARD

decay_card_pnbarpipi0_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1          PHSP;
    Enddecay

    Decay chi_c1
    1.0000 p+ anti-n0 pi- pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma           PHSP;
    Enddecay

    End
DECAYCARD

decay_card_pnbarpipi0_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2          PHSP;
    Enddecay

    Decay chi_c2
    1.0000 p+ anti-n0 pi- pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma           PHSP;
    Enddecay

    End
DECAYCARD

# Signal (charge conjugate): chi_cJ -> pbar n pi+ pi0
decay_card_pbarnpipi0_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0          PHSP;
    Enddecay

    Decay chi_c0
    1.0000 anti-p- n0 pi+ pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma           PHSP;
    Enddecay

    End
DECAYCARD

decay_card_pbarnpipi0_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1          PHSP;
    Enddecay

    Decay chi_c1
    1.0000 anti-p- n0 pi+ pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma           PHSP;
    Enddecay

    End
DECAYCARD

decay_card_pbarnpipi0_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2          PHSP;
    Enddecay

    Decay chi_c2
    1.0000 anti-p- n0 pi+ pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma           PHSP;
    Enddecay

    End
DECAYCARD

# -----------------------------------------------------------------------------
# Exclusive MC samples (one per chi_cJ state and per charge-conjugate channel)
# -----------------------------------------------------------------------------
exMC_pnbarpi_c0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic0_pnbarpi"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_pnbarpi_chic0
    config.cross_section   = :default
end

exMC_pnbarpi_c1 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic1_pnbarpi"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_pnbarpi_chic1
    config.cross_section   = :default
end

exMC_pnbarpi_c2 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic2_pnbarpi"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_pnbarpi_chic2
    config.cross_section   = :default
end

exMC_pbarnpi_c0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic0_pbarnpi"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_pbarnpi_chic0
    config.cross_section   = :default
end

exMC_pbarnpi_c1 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic1_pbarnpi"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_pbarnpi_chic1
    config.cross_section   = :default
end

exMC_pbarnpi_c2 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic2_pbarnpi"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_pbarnpi_chic2
    config.cross_section   = :default
end

exMC_pnbarpipi0_c0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic0_pnbarpipi0"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_pnbarpipi0_chic0
    config.cross_section   = :default
end

exMC_pnbarpipi0_c1 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic1_pnbarpipi0"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_pnbarpipi0_chic1
    config.cross_section   = :default
end

exMC_pnbarpipi0_c2 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic2_pnbarpipi0"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_pnbarpipi0_chic2
    config.cross_section   = :default
end

exMC_pbarnpipi0_c0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic0_pbarnpipi0"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_pbarnpipi0_chic0
    config.cross_section   = :default
end

exMC_pbarnpipi0_c1 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic1_pbarnpipi0"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_pbarnpipi0_chic1
    config.cross_section   = :default
end

exMC_pbarnpipi0_c2 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_gamma_chic2_pbarnpipi0"
    config.related_dataset = psip_data
    config.events          = 500_000
    config.decay_card      = decay_card_pbarnpipi0_chic2
    config.cross_section   = :default
end

# =============================================================================
# ALGORITHM 1: psi(2S) -> gamma chi_cJ, chi_cJ -> p nbar pi-
# =============================================================================
alg_name_pnbarpi = "PsipGammaChicJPNbarPi"
alg_pnbarpi = Algorithm.new(alg_name_pnbarpi)
alg_pnbarpi.set_header(["#{alg_name_pnbarpi}Alg/#{alg_name_pnbarpi}.h"])
           .set_constant({"ECMS" => [:double, 3.686]})

sel_pnbarpi = Selection.new
sel_pnbarpi.select_track {
             cos_theta 0.93    # |cos(theta)| < 0.93 for MDC tracks
             Vz        5.0     # within +-5 cm of the IP along the beam direction
             Vr        0.5     # within 0.5 cm of the IP in the plane perpendicular to the beam
             nChrp     "==1"   # exactly one proton
             nChrn     "==1"   # exactly one pi-
             nNet      "==0"
           }
           .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             angle_to_track    10.0    # > 10 deg from the nearest charged track
             energyThreshold_b 0.025   # barrel (|cos(theta)| < 0.80) : E > 25 MeV
             energyThreshold_e 0.050   # endcaps (0.86 < |cos(theta)| < 0.92) : E > 50 MeV
             nGam              ">=1"   # at least one radiative photon
           }
           .pid(method: :probability) {
             prob_cut 0.001
             identify :proton, against: [:kaon, :pion]   # p+ and anti-p- hypotheses
             identify :pion,   against: [:kaon, :proton]
             nprp "==1"        # exactly one proton candidate
             npim "==1"        # exactly one pi- candidate
           }
           # Photons from the interaction of the antiproton in the detector are
           # rejected: the photon must be more than 30 deg away from the
           # antiproton track (the > 20 deg cut against the antineutron and the
           # > 10 deg cut against the remaining particles are covered by the
           # baseline cuts above).
           .select_isolated_photon {
             angle_to_prm_track 30.0
             angle_to_prp_track 10.0
             nGam               ">=1"
           }
           # At least one photon with E > 80 MeV is required; the soft photons
           # are dropped from the candidate list before the kinematic fit.
           .for_each(:gamma) {
             where { energy < 0.08 }
             remove
           }
           # 1C kinematic fit under the psi(2S) -> gamma p nbar pi- hypothesis.
           # The antineutron is treated as the missing particle, its mass being
           # constrained to the nominal neutron mass (one constraint -> 1C).
           .kinematic_fit([:gamma, :prp, :pim]) {
             nominal
             miss_track_of :n_bar
             constrain_four_momentum
             chi2_cut 200   # loose BOSS cut; the paper's chi2_1C < 10 is applied in ROOT
           }

alg_pnbarpi
  .note(:photon_energy,
        "at least one photon with E > 80 MeV is required; implemented by removing all photon " \
        "candidates with E < 80 MeV before the kinematic fit.")
  .note(:antineutron_photon_angle,
        "alpha < 15 deg between the expected antineutron direction and the nearest photon is " \
        "required to purify the events containing an antineutron. The antineutron direction is " \
        "only defined from the missing momentum after the 1C fit, so this cut is applied at the " \
        "ROOT level.")
  .note(:lambda_veto,
        "M(p pi-) > 1.15 GeV/c^2 is required to remove background events containing a Lambda; " \
        "evaluated with the 1C-fit-corrected four-momenta and applied at the ROOT level.")
  .note(:transverse_momentum,
        "the transverse momentum of the proton is required to be > 0.3 GeV/c in order to reduce " \
        "the systematic uncertainty from the data-MC difference of the tracking efficiency at " \
        "low transverse momentum. Applied at the ROOT level.")
  .note(:intermediate_states,
        "structures are observed in the p pi-, nbar pi- and p nbar (threshold enhancement) " \
        "invariant mass spectra; they are taken into account in the signal generator (including " \
        "the polar angle distribution of the proton/neutron), and the chi_c0, chi_c1 and chi_c2 " \
        "MC samples are weighted with the amplitudes observed in data to obtain the detection " \
        "efficiency. The efficiency difference with respect to a phase-space sample re-weighted " \
        "in the two-dimensional (M(p pi-), M(nbar pi-)) plane is assigned as the systematic " \
        "uncertainty.")
  .note(:fit_details,
        "signal yields are extracted from unbinned maximum-likelihood fits to the p nbar pi- " \
        "invariant mass in 3.30-3.60 GeV/c^2 with the sum of three Breit-Wigner shapes " \
        "(natural widths fixed to the PDG values) convolved with a modified Gaussian resolution " \
        "plus a third-order Chebyshev background. Applied at the ROOT level.")
  .note(:background_estimation,
        "backgrounds are studied with the inclusive psi(2S) MC and the 3.65 GeV continuum data. " \
        "The dominant contributions are psi(2S) -> p nbar pi- pi0, psi(2S) -> p nbar pi-, " \
        "psi(2S) -> pi0 pi0 J/psi with J/psi -> p nbar pi-, and psi(2S) -> gamma chi_cJ with " \
        "chi_cJ -> p nbar pi- pi0, normalized with PDG branching fractions (our own measurement " \
        "for the last one). No peaking background is found in the signal region.")
  .note(:systematic_uncertainties,
        "tracking 2% (two tracks), PID 2% (3% for final states with anti-p- pi+), photon " \
        "detection 1% per photon, 1C kinematic fit 2.9% (2.7% for the charge conjugate channel), " \
        "alpha < 15 deg requirement 1.8%, intermediate states up to 8.2%, fitting procedure " \
        "(fit range, signal lineshape, resolution parameterization, resolution difference, " \
        "background shape), number of psi(2S) events 4%, and the PDG branching fractions of " \
        "psi(2S) -> gamma chi_cJ and pi0 -> gamma gamma.")

alg_pnbarpi.with_decay_card(decay_card_pnbarpi_chic0).apply(sel_pnbarpi)
alg_pnbarpi.execute_on([psip_data, psip_incMC, cont_data, cont_incMC,
                        exMC_pnbarpi_c0, exMC_pnbarpi_c1, exMC_pnbarpi_c2])

# =============================================================================
# ALGORITHM 2 (charge conjugate): psi(2S) -> gamma chi_cJ, chi_cJ -> pbar n pi+
# =============================================================================
alg_name_pbarnpi = "PsipGammaChicJPbarNPi"
alg_pbarnpi = Algorithm.new(alg_name_pbarnpi)
alg_pbarnpi.set_header(["#{alg_name_pbarnpi}Alg/#{alg_name_pbarnpi}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_pbarnpi = Selection.new
sel_pbarnpi.select_track {
             cos_theta 0.93
             Vz        5.0
             Vr        0.5
             nChrp     "==1"   # exactly one pi+
             nChrn     "==1"   # exactly one anti-p-
             nNet      "==0"
           }
           .select_photon {
             tdc_emc_start     0
             tdc_emc_end       14
             angle_to_track    10.0
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             nGam              ">=1"
           }
           .pid(method: :probability) {
             prob_cut 0.001
             identify :proton, against: [:kaon, :pion]   # anti-p- hypothesis
             identify :pion,   against: [:kaon, :proton] # pi+ hypothesis
             nprm "==1"        # exactly one antiproton candidate
             npip "==1"        # exactly one pi+ candidate
           }
           .select_isolated_photon {
             angle_to_prp_track 30.0
             angle_to_prm_track 10.0
             nGam               ">=1"
           }
           .for_each(:gamma) {
             where { energy < 0.08 }
             remove
           }
           # 1C kinematic fit under the psi(2S) -> gamma pbar n pi+ hypothesis,
           # the neutron being the missing particle (its mass is constrained to
           # the nominal neutron mass).
           .kinematic_fit([:gamma, :prm, :pip]) {
             nominal
             miss_track_of :n0
             constrain_four_momentum
             chi2_cut 200   # loose BOSS cut; the paper's chi2_1C < 10 is applied in ROOT
           }

alg_pbarnpi
  .note(:photon_energy,
        "at least one photon with E > 80 MeV is required; implemented by removing all photon " \
        "candidates with E < 80 MeV before the kinematic fit.")
  .note(:neutron_photon_angle,
        "alpha < 15 deg between the expected neutron direction and the nearest photon is required " \
        "to purify the events containing a neutron; applied at the ROOT level since the neutron " \
        "direction follows from the 1C fit.")
  .note(:lambda_veto,
        "M(pbar pi+) > 1.15 GeV/c^2 is required to remove anti-Lambda background; evaluated with " \
        "the 1C-fit-corrected four-momenta and applied at the ROOT level.")
  .note(:transverse_momentum,
        "the transverse momentum of the antiproton is required to be > 0.3 GeV/c to reduce the " \
        "systematic uncertainty from the data-MC tracking efficiency difference at low " \
        "transverse momentum. Applied at the ROOT level.")
  .note(:charge_conjugate_channel,
        "this is the charge conjugate of chi_cJ -> p nbar pi-; the selection is identical with " \
        "the roles of particle and antiparticle interchanged and is measured separately (right " \
        "column of Table I of the paper).")
  .note(:fit_details,
        "signal yields are extracted from unbinned maximum-likelihood fits to the pbar n pi+ " \
        "invariant mass in 3.30-3.60 GeV/c^2, as for the p nbar pi- channel. Applied at the " \
        "ROOT level.")
  .note(:systematic_uncertainties,
        "same sources as for the chi_cJ -> p nbar pi- channel; the PID uncertainty is 3% for " \
        "final states containing anti-p- pi+.")

alg_pbarnpi.with_decay_card(decay_card_pbarnpi_chic0).apply(sel_pbarnpi)
alg_pbarnpi.execute_on([psip_data, psip_incMC, cont_data, cont_incMC,
                        exMC_pbarnpi_c0, exMC_pbarnpi_c1, exMC_pbarnpi_c2])

# =============================================================================
# ALGORITHM 3: psi(2S) -> gamma chi_cJ, chi_cJ -> p nbar pi- pi0
# =============================================================================
alg_name_pnbarpipi0 = "PsipGammaChicJPNbarPiPi0"
alg_pnbarpipi0 = Algorithm.new(alg_name_pnbarpipi0)
alg_pnbarpipi0.set_header(["#{alg_name_pnbarpipi0}Alg/#{alg_name_pnbarpipi0}.h"])
              .set_constant({"ECMS" => [:double, 3.686]})

sel_pnbarpipi0 = Selection.new
sel_pnbarpipi0.select_track {
                cos_theta 0.93
                Vz        5.0
                Vr        0.5
                nChrp     "==1"   # exactly one proton
                nChrn     "==1"   # exactly one pi-
                nNet      "==0"
              }
              .select_photon {
                tdc_emc_start     0
                tdc_emc_end       14
                angle_to_track    10.0
                energyThreshold_b 0.025
                energyThreshold_e 0.050
                nGam              ">=3"   # three photons: pi0 pair + radiative photon
              }
              .pid(method: :probability) {
                prob_cut 0.001
                identify :proton, against: [:kaon, :pion]
                identify :pion,   against: [:kaon, :proton]
                nprp "==1"
                npim "==1"
              }
              # Same photon isolation criteria as for the p nbar pi- channel.
              .select_isolated_photon {
                angle_to_prm_track 30.0
                angle_to_prp_track 10.0
                nGam               ">=3"
              }
              .for_each(:gamma) {
                where { energy < 0.08 }
                remove
              }
              # Nominal 1C kinematic fit under the psi(2S) -> gamma gamma gamma
              # p nbar pi- hypothesis; the mass of the missing particle
              # (antineutron) is constrained to the nominal neutron mass, and all
              # possible three-photon combinations are tried, the smallest chi2
              # being retained by the fit.
              .kinematic_fit([:gamma, :gamma, :gamma, :prp, :pim]) {
                nominal
                miss_track_of :n_bar
                constrain_four_momentum
                chi2_cut 200   # loose BOSS cut; the paper's chi2 < 10 is applied in ROOT
              }
              # Competing-hypothesis fits (no chi2_cut, no nominal): the chi2
              # values are stored for the ROOT-level veto
              # chi2(gamma gamma gamma ...) < chi2(gamma gamma ...) and
              # chi2(gamma gamma gamma ...) < chi2(gamma gamma gamma gamma ...).
              .kinematic_fit([:gamma, :gamma, :prp, :pim]) {
                miss_track_of :n_bar
                constrain_four_momentum
              }
              .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :prp, :pim]) {
                miss_track_of :n_bar
                constrain_four_momentum
              }

alg_pnbarpipi0
  .note(:photon_energy,
        "the same photon isolation criteria as in the chi_cJ -> p nbar pi- channel are applied; " \
        "photon candidates with E < 80 MeV are removed before the kinematic fit.")
  .note(:pi0_selection,
        "pi0 candidates are formed from any pair of photon candidates kinematically fitted to the " \
        "pi0 mass with chi2 < 20, and at least one pi0 candidate is required. The 1C fit is " \
        "carried out over all possible three-photon combinations; when more than one pi0 " \
        "candidate can be formed from the three photons, the pair with invariant mass closest to " \
        "the pi0 mass is assigned to the pi0.")
  .note(:competing_hypotheses,
        "chi2(gamma gamma gamma p nbar pi-) < chi2(gamma gamma p nbar pi-) and " \
        "chi2(gamma gamma gamma p nbar pi-) < chi2(gamma gamma gamma gamma p nbar pi-) are " \
        "required to suppress the gamma gamma p nbar pi- and gamma gamma gamma gamma p nbar pi- " \
        "backgrounds; the chi2 values are stored by the two additional non-nominal kinematic " \
        "fits and the comparison is performed at the ROOT level.")
  .note(:background_veto,
        "the psi(2S) -> pi0 pi0 J/psi background is suppressed for events with at least four " \
        "photons: the pi0 pi0 combination minimising " \
        "Delta = sqrt((m(g1 g2) - m_pi0)^2 + (m(g3 g4) - m_pi0)^2) is selected and " \
        "|M_recoil(pi0 pi0) - M_J/psi| > 50 MeV/c^2 is required. Applied at the ROOT level.")
  .note(:antineutron_photon_angle,
        "alpha < 15 deg between the expected antineutron direction and the nearest photon is " \
        "required; applied at the ROOT level (the antineutron direction comes from the 1C fit).")
  .note(:transverse_momentum,
        "the transverse momentum of the proton is required to be > 0.3 GeV/c in order to reduce " \
        "the systematic uncertainty from the data-MC tracking efficiency difference at low " \
        "transverse momentum. Applied at the ROOT level.")
  .note(:intermediate_states,
        "no obvious N* intermediate state is observed in chi_cJ -> p nbar pi- pi0, but a " \
        "significant rho+- signal is seen in the pi+- pi0 invariant mass spectrum while the " \
        "efficiency MC is generated according to phase space. A chi_cJ -> p nbar rho- sample " \
        "generated with the correct rho+- -> pi+- pi0 angular distribution is used to estimate " \
        "the corresponding systematic uncertainty.")
  .note(:fit_details,
        "signal yields are extracted from unbinned maximum-likelihood fits to the p nbar pi- pi0 " \
        "invariant mass in 3.30-3.64 GeV/c^2 with the sum of three Breit-Wigner shapes " \
        "(natural widths fixed to the PDG values) convolved with a modified Gaussian resolution " \
        "plus a third-order Chebyshev background. Applied at the ROOT level.")
  .note(:background_estimation,
        "the dominant backgrounds, studied with the inclusive psi(2S) MC and the 3.65 GeV " \
        "continuum data, are psi(2S) -> p nbar pi- pi0 (about 25%), psi(2S) -> pi0 pi0 J/psi " \
        "with J/psi -> p nbar pi- (10%), and psi(2S) -> gamma chi_cJ with chi_cJ -> p nbar pi- " \
        "(1%); no peaking background is found in the signal region.")
  .note(:systematic_uncertainties,
        "tracking 2%, PID 2% (3% for anti-p- pi+), photon detection 3% (three photons), 1C " \
        "kinematic fit 2.9% (2.7% for the charge conjugate channel), alpha < 15 deg 1.8%, pi0 " \
        "reconstruction 0.7%, intermediate states, fitting procedure, number of psi(2S) events " \
        "4%, and the PDG branching fractions of psi(2S) -> gamma chi_cJ and pi0 -> gamma gamma.")

alg_pnbarpipi0.with_decay_card(decay_card_pnbarpipi0_chic0).apply(sel_pnbarpipi0)
alg_pnbarpipi0.execute_on([psip_data, psip_incMC, cont_data, cont_incMC,
                           exMC_pnbarpipi0_c0, exMC_pnbarpipi0_c1, exMC_pnbarpipi0_c2])

# =============================================================================
# ALGORITHM 4 (charge conjugate): psi(2S) -> gamma chi_cJ, chi_cJ -> pbar n pi+ pi0
# =============================================================================
alg_name_pbarnpipi0 = "PsipGammaChicJPbarNPiPi0"
alg_pbarnpipi0 = Algorithm.new(alg_name_pbarnpipi0)
alg_pbarnpipi0.set_header(["#{alg_name_pbarnpipi0}Alg/#{alg_name_pbarnpipi0}.h"])
             .set_constant({"ECMS" => [:double, 3.686]})

sel_pbarnpipi0 = Selection.new
sel_pbarnpipi0.select_track {
                cos_theta 0.93
                Vz        5.0
                Vr        0.5
                nChrp     "==1"   # exactly one pi+
                nChrn     "==1"   # exactly one anti-p-
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
              .pid(method: :probability) {
                prob_cut 0.001
                identify :proton, against: [:kaon, :pion]   # anti-p- hypothesis
                identify :pion,   against: [:kaon, :proton] # pi+ hypothesis
                nprm "==1"
                npip "==1"
              }
              .select_isolated_photon {
                angle_to_prp_track 30.0
                angle_to_prm_track 10.0
                nGam               ">=3"
              }
              .for_each(:gamma) {
                where { energy < 0.08 }
                remove
              }
              # Nominal 1C kinematic fit under the psi(2S) -> gamma gamma gamma
              # pbar n pi+ hypothesis, the neutron being the missing particle.
              .kinematic_fit([:gamma, :gamma, :gamma, :prm, :pip]) {
                nominal
                miss_track_of :n0
                constrain_four_momentum
                chi2_cut 200   # loose BOSS cut; the paper's chi2 < 10 is applied in ROOT
              }
              # Competing-hypothesis fits storing the chi2 of the
              # gamma gamma pbar n pi+ and gamma gamma gamma gamma pbar n pi+
              # hypotheses for the ROOT-level veto.
              .kinematic_fit([:gamma, :gamma, :prm, :pip]) {
                miss_track_of :n0
                constrain_four_momentum
              }
              .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :prm, :pip]) {
                miss_track_of :n0
                constrain_four_momentum
              }

alg_pbarnpipi0
  .note(:photon_energy,
        "the same photon isolation criteria as in the chi_cJ -> p nbar pi- channel are applied; " \
        "photon candidates with E < 80 MeV are removed before the kinematic fit.")
  .note(:pi0_selection,
        "pi0 candidates are formed from any pair of photon candidates kinematically fitted to the " \
        "pi0 mass with chi2 < 20 and at least one pi0 candidate is required; the 1C fit is " \
        "carried out over all possible three-photon combinations, the pair with invariant mass " \
        "closest to the pi0 mass being assigned to the pi0.")
  .note(:competing_hypotheses,
        "chi2(gamma gamma gamma pbar n pi+) < chi2(gamma gamma pbar n pi+) and " \
        "chi2(gamma gamma gamma pbar n pi+) < chi2(gamma gamma gamma gamma pbar n pi+) are " \
        "required; the chi2 values are stored by the two additional non-nominal kinematic fits " \
        "and the comparison is applied at the ROOT level.")
  .note(:background_veto,
        "the psi(2S) -> pi0 pi0 J/psi background is suppressed for events with at least four " \
        "photons by selecting the pi0 pi0 combination that minimises " \
        "Delta = sqrt((m(g1 g2) - m_pi0)^2 + (m(g3 g4) - m_pi0)^2) and requiring " \
        "|M_recoil(pi0 pi0) - M_J/psi| > 50 MeV/c^2. Applied at the ROOT level.")
  .note(:neutron_photon_angle,
        "alpha < 15 deg between the expected neutron direction and the nearest photon is " \
        "required; applied at the ROOT level.")
  .note(:transverse_momentum,
        "the transverse momentum of the antiproton is required to be > 0.3 GeV/c to reduce the " \
        "systematic uncertainty from the data-MC tracking efficiency difference at low " \
        "transverse momentum. Applied at the ROOT level.")
  .note(:charge_conjugate_channel,
        "this is the charge conjugate of chi_cJ -> p nbar pi- pi0; the selection is identical " \
        "with the roles of particle and antiparticle interchanged and is measured separately " \
        "(right column of Table II of the paper).")
  .note(:intermediate_states,
        "a significant rho+- signal is observed in the pi+- pi0 invariant mass spectrum while " \
        "the efficiency MC is generated according to phase space; a chi_cJ -> pbar n rho+ sample " \
        "is used to estimate the corresponding systematic uncertainty.")
  .note(:fit_details,
        "signal yields are extracted from unbinned maximum-likelihood fits to the pbar n pi+ pi0 " \
        "invariant mass in 3.30-3.64 GeV/c^2, as for the p nbar pi- pi0 channel. Applied at the " \
        "ROOT level.")
  .note(:systematic_uncertainties,
        "same sources as for the chi_cJ -> p nbar pi- pi0 channel; the PID uncertainty is 3% for " \
        "final states containing anti-p- pi+.")

alg_pbarnpipi0.with_decay_card(decay_card_pbarnpipi0_chic0).apply(sel_pbarnpipi0)
alg_pbarnpipi0.execute_on([psip_data, psip_incMC, cont_data, cont_incMC,
                           exMC_pbarnpipi0_c0, exMC_pbarnpipi0_c1, exMC_pbarnpipi0_c2])
