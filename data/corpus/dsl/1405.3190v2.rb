# Analysis: Search for the radiative transitions psi(3770) -> gamma eta_c and
# gamma eta_c(2S) through eta_c(eta_c(2S)) -> K_S0 K+ pi- (and psi(3770) -> gamma chi_c1,
# chi_c1 -> K_S0 K+ pi- as a cross-check), using 2.92 fb^-1 of data at sqrt(s) = 3.773 GeV.
# Three independent signal modes (eta_c, eta_c(2S), chi_c1) share the same final state and
# the same selection, but use different decay cards and algorithm objects.

### Datasets ###
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ---------------------------------------------------------------------------
# Decay cards (EvtGen format). The K_S0 K+- pi-+ sub-decays are generated
# according to the measured Dalitz-plot distributions in the analysis; PHSP is
# used here as the expressible fallback.
# ---------------------------------------------------------------------------

# Signal mode 1: psi(3770) -> gamma eta_c, eta_c -> K_S0 K+ pi-
decay_card_etac = <<~DECAYCARD
    Decay psi(3770)
    1.0000  gamma eta_c                         PHSP;
    Enddecay

    Decay eta_c
    1.0000  K_S0 K+ pi-                         PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+ pi-                             PHSP;
    Enddecay

    End
DECAYCARD

# Signal mode 2: psi(3770) -> gamma eta_c(2S), eta_c(2S) -> K_S0 K+ pi-
decay_card_etac2S = <<~DECAYCARD
    Decay psi(3770)
    1.0000  gamma eta_c(2S)                     PHSP;
    Enddecay

    Decay eta_c(2S)
    1.0000  K_S0 K+ pi-                         PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+ pi-                             PHSP;
    Enddecay

    End
DECAYCARD

# Cross-check mode: psi(3770) -> gamma chi_c1, chi_c1 -> K_S0 K+ pi-
decay_card_chic1 = <<~DECAYCARD
    Decay psi(3770)
    1.0000  gamma chi_c1                         PHSP;
    Enddecay

    Decay chi_c1
    1.0000  K_S0 K+ pi-                          PHSP;
    Enddecay

    Decay K_S0
    1.0000  pi+ pi-                              PHSP;
    Enddecay

    End
DECAYCARD

exMC_etac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3770_gamma_etac_KSKpi"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_etac
  config.cross_section   = :default
end

exMC_etac2S = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3770_gamma_etac2S_KSKpi"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_etac2S
  config.cross_section   = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3770_gamma_chic1_KSKpi"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_chic1
  config.cross_section   = :default
end

# Helper: the three modes share an identical final state gamma K_S0 K+- pi-+ and the
# same selection chain; only the decay card / algorithm name differ.
def build_selection
  sel = Selection.new
  sel.select_track {                 # at least four charged tracks (K_S0 -> pi+pi- + K + pi)
        cos_theta 0.93               # |cos(theta)| < 0.93 (MDC fiducial volume)
        Vz        10.0               # within +-10 cm of the IP along the beam direction
        Vr        1.0                # within 1 cm in the radial direction
        nChrp     ">=2"
        nChrn     ">=2"
        nNet      "==0"              # zero net charge
      }
     .select_photon {                # at least one good photon (the radiative photon)
        tdc_emc_start     0
        tdc_emc_end       14         # EMC timing within 700 ns of the collision
        angle_to_track    20.0       # more than 20 degrees from any charged track
        energyThreshold_b 0.025      # E > 25 MeV in the barrel (|cos(theta)| < 0.8)
        energyThreshold_e 0.050      # E > 50 MeV in the end-cap (0.86 < |cos(theta)| < 0.92)
        nGam              ">=1"
      }
     .assign({:chrgp => :pip, :chrgn => :pim})   # K_S0 daughters are assumed to be pions
     .secondary_vertex_fit([:pip, :pim]) {       # K_S0 -> pi+ pi- secondary vertex fit
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
      }
     .pid(method: :chi2_sum) {        # PID by minimizing chi2_PID over K/pi hypotheses
        chi_min_cut 4
        identify :kaon, :pion        # combinatorial K/pi assignment for the remaining tracks
      }
     # Nominal 4C kinematic fit to the gamma K_S0 K+- pi-+ hypothesis. The species of the
     # final-state particles and the best photon are chosen by minimizing
     # chi2_total = chi2_4C + chi2_PID(K) + chi2_PID(pi); the DSL selects the
     # smallest-chi2 combination automatically.
     .kinematic_fit([:gamma, :K_S0, :kp, :pim, :pip, :km]) {
        nominal
        constrain_four_momentum
        chi2_cut 200                 # loose BOSS cut; the paper's chi2_4C < 20 applied in ROOT
      }
     # Alternative 4C fit under the gamma gamma K_S0 K+- pi-+ hypothesis (no chi2_cut and no
     # nominal): stores chi2_4C of the extra-photon hypothesis, so that the pi0 background
     # suppression chi2_4C(gamma ...) < chi2_4C(gamma gamma ...) can be applied in ROOT.
     .kinematic_fit([:gamma, :gamma, :K_S0, :kp, :pim, :pip, :km]) {
        constrain_four_momentum
      }
  sel
end

# ---------------------------------------------------------------------------
# Signal mode 1: psi(3770) -> gamma eta_c -> gamma K_S0 K+ pi-
# ---------------------------------------------------------------------------
alg_etac = Algorithm.new("PsiPrimeToGammaEtac")
alg_etac.set_header(["PsiPrimeToGammaEtacAlg/PsiPrimeToGammaEtac.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
alg_etac.with_decay_card(decay_card_etac).apply(build_selection)

# ---------------------------------------------------------------------------
# Signal mode 2: psi(3770) -> gamma eta_c(2S) -> gamma K_S0 K+ pi-
# ---------------------------------------------------------------------------
alg_etac2S = Algorithm.new("PsiPrimeToGammaEtac2S")
alg_etac2S.set_header(["PsiPrimeToGammaEtac2SAlg/PsiPrimeToGammaEtac2S.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
alg_etac2S.with_decay_card(decay_card_etac2S).apply(build_selection)

# ---------------------------------------------------------------------------
# Cross-check mode: psi(3770) -> gamma chi_c1 -> gamma K_S0 K+ pi-
# ---------------------------------------------------------------------------
alg_chic1 = Algorithm.new("PsiPrimeToGammaChic1")
alg_chic1.set_header(["PsiPrimeToGammaChic1Alg/PsiPrimeToGammaChic1.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
alg_chic1.with_decay_card(decay_card_chic1).apply(build_selection)

# BOSS-side procedures that cannot be expressed with the current DSL capabilities
[alg_etac, alg_etac2S, alg_chic1].each do |alg|
  alg.note(:ks_mass_window, "K_S0 candidates: secondary vertex fits to all pairs of oppositely " \
        "charged tracks (assumed pions); the combination with the best fit quality is kept when " \
        "|M(pi+pi-) - M(K_S0)| < 10 MeV/c^2 and the decay length is more than twice the vertex " \
        "resolution. The fitted K_S0 information is used as input to the kinematic fit.")
     .note(:background_veto, "D0 veto: |M(K+-pi-+) - M(D0)| > 3 sigma to remove psi(3770) -> D0 D0bar, " \
        "D0 -> pi0 K_S0, D0 -> pi+ K- (and charge conjugate); pi0 veto: chi2_4C(gamma K_S0 K+-pi-+) < " \
        "chi2_4C(gamma gamma K_S0 K+-pi-+). Both are evaluated from kinematic-fit-corrected " \
        "quantities and applied in ROOT; the alternative-fit chi2_4C is stored for the pi0 veto.")
     .note(:helix_correction, "Track helix parameters (phi0, kappa, tan lambda) corrected for the " \
        "data/MC difference; correction factors extracted from pull distributions using the control " \
        "sample J/psi -> phi f0(980), phi -> K+K-, f0(980) -> pi+pi-. Efficiency differences with and " \
        "without the correction (3.9%, 5.5%, 5.3% for eta_c, eta_c(2S), chi_c1) are taken as " \
        "systematic uncertainties.")
     .note(:efficiency_curve, "Efficiencies determined from MC generated with the expected angular " \
        "distributions for psi(3770) -> gamma eta_c, gamma eta_c(2S), gamma chi_c1 and the measured " \
        "Dalitz-plot distributions for eta_c, eta_c(2S), chi_c1 -> K_S0 K+- pi-+ (from Belle). " \
        "eps = 27.87% (eta_c), 25.24% (eta_c(2S)), 28.46% (chi_c1).")
     .note(:pid_correction_method, "PID uses chi2_PID(i) = ((dE/dx_meas - dE/dx_exp)/sigma_dE/dx)^2 + " \
        "((TOF_meas - TOF_exp)/sigma_TOF)^2 evaluated for the pion, kaon and proton hypotheses; the " \
        "final state assignment minimizes chi2_total = chi2_4C + chi2_PID(K) + chi2_PID(pi).")
     .note(:background_veto, "Backgrounds from the continuum e+e- -> q qbar (smooth), from " \
        "e+e- -> pi0 K_S0 K+-pi-+ (measured from data), from e+e- -> (gamma_ISR/gamma_FSR) K_S0 K+-pi-+ " \
        "(MC normalized by luminosity) and from the radiative tail of the psi(3686) (peaking in the " \
        "eta_c, eta_c(2S) and chi_c1 mass regions) are subtracted/simulated in the fit to the " \
        "M(K_S0 K+-pi-+) spectrum.")
end

alg_etac.execute_on([psi3770_data, psi3770_incMC, exMC_etac])
alg_etac2S.execute_on([psi3770_data, psi3770_incMC, exMC_etac2S])
alg_chic1.execute_on([psi3770_data, psi3770_incMC, exMC_chic1])
