# =====================================================================
# BESIII: Search for the eta_c(2S)/h_c -> p pbar decays and measurements
#         of the chi_cJ -> p pbar branching fractions
# Data: 1.06 x 10^8 psi(3686) events at sqrt(s) = 3.686 GeV
#       plus a 44 pb^-1 continuum sample at sqrt(s) = 3.65 GeV
#
# Two final states are selected:
#   (A) psi(3686) -> gamma eta_c(2S) / gamma chi_cJ -> gamma p pbar
#   (B) psi(3686) -> pi0 h_c -> pi0 p pbar -> gamma gamma p pbar
# Final state (A) is shared by the eta_c(2S) and the three chi_cJ signals
# (identical topology and identical selection), so a single Algorithm with
# one common Selection serves all four; only the generated line shapes
# (E_gamma^3 x BW x damping factor for eta_c(2S)/chi_cJ) differ, which is
# handled at the exclusive-MC level.
# =====================================================================

### Dataset description ###
psip_data      = DatasetManager.real_data.find("709_3686")   # psi(3686) data at 3.686 GeV
psip_incMC     = DatasetManager.inclusive_mc.find("709_3686") # 1.06 x 10^8 inclusive psi(3686) MC events
cont_data_3650 = DatasetManager.real_data.find("709_3650")   # 44 pb^-1 continuum data at 3.65 GeV

# ---------------------------------------------------------------------
# Decay cards (EvtGen)
# ---------------------------------------------------------------------

# Channel A: psi(3686) -> gamma eta_c(2S), eta_c(2S) -> p pbar (phase space)
decay_card_etac2S = <<~DECAYCARD
    Decay psi(3686)
    1.0000 gamma eta_c(2S) PHSP;
    Enddecay

    Decay eta_c(2S)
    1.0000 p+ anti-p- PHSP;
    Enddecay

    End
DECAYCARD

# Channel A': psi(3686) -> gamma chi_cJ, chi_cJ -> p pbar
# The angular distribution of the protons follows 1 + alpha cos^2(theta) in
# the chi_cJ helicity frame, with alpha taken from measured data
# (alpha = 0.09 +- 0.11, 0.12 +- 0.20, -0.26 +- 0.17 for J = 0, 1, 2);
# the generation model is approximated by phase space here.
decay_card_chic0 = <<~DECAYCARD
    Decay psi(3686)
    1.0000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.0000 p+ anti-p- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic1 = <<~DECAYCARD
    Decay psi(3686)
    1.0000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 p+ anti-p- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic2 = <<~DECAYCARD
    Decay psi(3686)
    1.0000 gamma chi_c2 PHSP;
    Enddecay

    Decay chi_c2
    1.0000 p+ anti-p- PHSP;
    Enddecay

    End
DECAYCARD

# Channel B: psi(3686) -> pi0 h_c, h_c -> p pbar, pi0 -> gamma gamma
# (h_c -> p pbar generated according to phase space)
decay_card_hc = <<~DECAYCARD
    Decay psi(3686)
    1.0000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 p+ anti-p- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Background: psi(3686) -> pi0 p pbar (irreducible background in the
# pi0 h_c channel, and the source of the pi0 p pbar background in the
# gamma p pbar channel through epsilon(gamma p pbar)/epsilon(pi0 p pbar))
decay_card_bkg_pi0ppbar = <<~DECAYCARD
    Decay psi(3686)
    1.0000 pi0 p+ anti-p- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples for the signal channels
exMC_etac2S = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_etac2S_ppbar"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_etac2S
  config.cross_section   = :default
end

exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chic0_ppbar"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_chic0
  config.cross_section   = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chic1_ppbar"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_chic1
  config.cross_section   = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamma_chic2_ppbar"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_chic2
  config.cross_section   = :default
end

exMC_hc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0_hc_ppbar"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_hc
  config.cross_section   = :default
end

exMC_bkg_pi0ppbar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0_ppbar_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_pi0ppbar
  config.cross_section   = :default
end

# ---------------------------------------------------------------------
### Event selection (BOSS) ###
# ---------------------------------------------------------------------

# =====================================================================
# Algorithm A -- psi(3686) -> gamma eta_c(2S) / gamma chi_cJ -> gamma p pbar
#   two oppositely charged tracks (both identified as protons), at least
#   one good photon, 4C kinematic fit to gamma p pbar with chi2_4C < 40
# =====================================================================
algA_name = "gammappbar"
algA = Algorithm.new(algA_name)
algA.set_header(["#{algA_name}Alg/#{algA_name}.h"])
    .set_constant({"ECMS" => [:double, 3.686]})
    .set_alias({"std::vector<double>" => "Vdouble"})

selA = Selection.new
selA.select_track do
       cos_theta 0.93   # within the MDC polar-angle coverage, |cos(theta)| < 0.93
       Vz        10.0   # within 10 cm of the IP along the beam direction
       Vr        1.0    # within 1 cm of the beam line in the radial direction
       nChrp     "==1"  # one positive track
       nChrn     "==1"  # one negative track
       nNet      "==0"
     end
     # Photons: isolated EMC showers, >= 25 MeV in both barrel and endcap,
     # EMC timing in coincidence with the collision, and at least 10 degrees
     # from any charged track (the proton/anti-proton isolation of 15/25
     # degrees is applied in the select_isolated_photon step below)
     .select_photon do
       tdc_emc_start     0
       tdc_emc_end       14
       angle_to_track    10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.025
       nGam              ">=1"
     end
     # At least 15 (25) degrees away from the proton (anti-proton) candidate
     .select_isolated_photon do
       angle_to_prp_track 15.0
       angle_to_prm_track 25.0
       nGam              ">=1"
     end
     # Both tracks must be positively identified as protons (one p, one p-bar),
     # using the TOF (and dE/dx) likelihood
     .pid(method: :probability) do
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       nprp "==1"
       nprm "==1"
     end
     # 4C kinematic fit of the gamma p pbar candidates to the total initial
     # four-momentum of the colliding beams; when more than one photon
     # candidate exists, the combination minimizing chi2_4C is retained
     .kinematic_fit([:gamma, :prp, :prm]) do
       nominal
       constrain_four_momentum
       chi2_cut 40   # chi2_4C < 40
     end
     .note(:kinematic_fit,
           "a second, 3C kinematic fit (magnitude of the photon momentum allowed to float) is used to determine the signal yields; it keeps the psi(3686) -> p pbar background peak at its correct position and separates it from the eta_c(2S) signal. The psi(3686) -> p pbar / gamma_FSR p pbar and non-resonant backgrounds are described by MC shapes in the M(p pbar) fit")
     .note(:helix_correction,
           "corrections to the track-helix parameters are applied before the kinematic fit; the efficiency difference with and without the correction is taken as the systematic uncertainty")
     .note(:efficiency_curve,
           "the selection efficiency for chi_cJ depends on the proton helicity angular distribution 1 + alpha cos^2(theta); alpha is varied by +-1 sigma to estimate the systematic uncertainty")

algA.with_decay_card(decay_card_etac2S).apply(selA)

# =====================================================================
# Algorithm B -- psi(3686) -> pi0 h_c -> pi0 p pbar -> gamma gamma p pbar
#   two oppositely charged tracks, at least one identified as a proton,
#   at least two good photons, 4C kinematic fit to gamma gamma p pbar
#   with chi2_4C < 40 and 0.11 < M(gamma gamma) < 0.15 GeV/c^2
# =====================================================================
algB_name = "pi0ppbar"
algB = Algorithm.new(algB_name)
algB.set_header(["#{algB_name}Alg/#{algB_name}.h"])
    .set_constant({"ECMS" => [:double, 3.686]})
    .set_alias({"std::vector<double>" => "Vdouble"})

selB = Selection.new
selB.select_track do
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     "==1"
       nChrn     "==1"
       nNet      "==0"
     end
     # At least two photons to form the pi0
     .select_photon do
       tdc_emc_start     0
       tdc_emc_end       14
       angle_to_track    10.0
       energyThreshold_b 0.025
       energyThreshold_e 0.025
       nGam              ">=2"
     end
     # Proton/anti-proton isolation of the photons (15 / 25 degrees)
     .select_isolated_photon do
       angle_to_prp_track 15.0
       angle_to_prm_track 25.0
       nGam              ">=2"
     end
     # For the gamma gamma p pbar final state only one of the two tracks is
     # required to be positively identified as a proton or anti-proton
     .pid(method: :probability) do
       prob_cut 0.001
       identify :proton, against: [:kaon, :pion]
       nprp ">=1"
     end
     # 4C kinematic fit of the gamma gamma p pbar candidates to the total
     # initial four-momentum; the best photon pair is chosen by minimizing
     # chi2_4C
     .kinematic_fit([:gamma, :gamma, :prp, :prm]) do
       nominal
       constrain_four_momentum
       # pi0 mass window on the two selected photons
       invariant_mass_of(:gamma, :gamma).within(0.11, 0.15)
       chi2_cut 40   # chi2_4C < 40
     end
     .note(:background_veto,
           "the requirement M(p pbar gamma_high) < 3.66 GeV/c^2, where gamma_high is the higher-energy photon, is applied on the 3C-fit quantities to remove the psi(3686) -> gamma chi_cJ (J = 1, 2) background; the 3C fit floats the momentum of the lower-energy photon and is not expressible as a DSL constraint")
     .note(:background_veto,
           "the irreducible psi(3686) -> pi0 p pbar background is described by an ARGUS function in the M(p pbar) fit; the pi0 mass requirement contributes a 3% systematic uncertainty")
     .note(:kinematic_fit,
           "the h_c signal yield is extracted from an unbinned maximum-likelihood fit to M(p pbar) after the 4C selection, using a MC signal shape convolved with a smearing Gaussian whose parameters are determined from psi(3686) -> pi0 J/psi, J/psi -> p pbar")

algB.with_decay_card(decay_card_hc).apply(selB)

# ---------------------------------------------------------------------
# Execution on data + inclusive MC + exclusive MC
# ---------------------------------------------------------------------
datasets = [psip_data, psip_incMC, cont_data_3650,
            exMC_etac2S, exMC_chic0, exMC_chic1, exMC_chic2,
            exMC_hc, exMC_bkg_pi0ppbar]
algA.execute_on(datasets)
algB.execute_on(datasets)
