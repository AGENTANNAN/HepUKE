# BESIII paper arXiv:1208.2320v2
# Experimental study of psi' decays to K+ K- pi0 and K+ K- eta, based on
# (106 +/- 4) x 10^6 psi' events (156.4 pb^-1) collected with BESIII at
# sqrt(s) = 3.686 GeV. A 2.9 fb^-1 (43 pb^-1) sample at sqrt(s) = 3.773 GeV
# (3.65 GeV) is used for continuum / QED background studies.
# Final state: gamma gamma K+ K- (two charged tracks identified as kaons and
# two photons), no charged tracks other than the two kaons.
# Nominal fit: 4C kinematic fit to the psi' -> gamma gamma K+ K- hypothesis.
# All of psi' -> K+ K- pi0 (inclusive, incl. K*(892)+- K-+, K2*(1430)+- K-+,
# K*(1680)+- K-+, rho(1700) pi0 and the non-resonant K+ K- pi0 mode),
# psi' -> K+ K- eta, psi' -> eta phi and psi' -> pi0 phi share this identical
# final state and identical selection chain, so a single Algorithm is used.

### Dataset description ###
psip_data     = DatasetManager.real_data.find("709_3686")     # psi(2S) data, (106 +/- 4) M events (156.4 pb^-1)
psip_incMC    = DatasetManager.inclusive_mc.find("709_3686")  # psi(2S) inclusive MC (kkmc + besevtgen + lundcharm)
psi3770_data  = DatasetManager.real_data.find("712_3773")     # psi(3770) data at sqrt(s) = 3.773 GeV (2.9 fb^-1)
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")
off_data      = DatasetManager.real_data.find("709_3650")     # continuum data at sqrt(s) = 3.65 GeV (43 pb^-1)
off_incMC     = DatasetManager.inclusive_mc.find("709_3650")

### Decay cards (EvtGen format) ###
# Signal: psi' -> K+ K- pi0, pi0 -> gamma gamma (inclusive K+ K- pi0 final state)
decay_card_kk_pi0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000  K+  K-  pi0            PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma           PHSP;
    Enddecay

    End
DECAYCARD

# Signal: psi' -> K+ K- eta, eta -> gamma gamma
decay_card_kk_eta = <<~DECAYCARD
    Decay psi(2S)
    1.0000  K+  K-  eta            PHSP;
    Enddecay

    Decay eta
    1.0000  gamma  gamma           PHSP;
    Enddecay

    End
DECAYCARD

# Signal: psi' -> pi0 phi, phi -> K+ K-, pi0 -> gamma gamma.
# Generated with a 1 + cos^2(theta) angular distribution for psi' -> pi0 phi.
decay_card_pi0_phi = <<~DECAYCARD
    Decay psi(2S)
    1.0000  pi0  phi               PHSP;
    Enddecay

    Decay phi
    1.0000  K+  K-                 VSS;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma           PHSP;
    Enddecay

    End
DECAYCARD

# Signal: psi' -> eta phi, phi -> K+ K-, eta -> gamma gamma
decay_card_eta_phi = <<~DECAYCARD
    Decay psi(2S)
    1.0000  eta  phi               PHSP;
    Enddecay

    Decay phi
    1.0000  K+  K-                 VSS;
    Enddecay

    Decay eta
    1.0000  gamma  gamma           PHSP;
    Enddecay

    End
DECAYCARD

# QED background: e+e- -> gamma* -> K+ K- pi0 (estimated from the 3.773 and 3.65 GeV data)
decay_card_continuum_kk_pi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0000  K+  K-  pi0            PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma           PHSP;
    Enddecay

    End
DECAYCARD

# QED background: e+e- -> gamma* -> eta phi
decay_card_continuum_eta_phi = <<~DECAYCARD
    Decay psi(4260)
    1.0000  eta  phi               PHSP;
    Enddecay

    Decay phi
    1.0000  K+  K-                 VSS;
    Enddecay

    Decay eta
    1.0000  gamma  gamma           PHSP;
    Enddecay

    End
DECAYCARD

# Background (studied with the inclusive MC sample): the dominant background for
# psi' -> eta K+ K- comes from psi' -> gamma chi_c2, chi_c2 -> K+ K- pi0 / K+ K- eta,
# and from psi' -> gamma gamma_FSR K+ K- with an undetected / misreconstructed photon.
decay_card_bkg_chi_c2_kk_pi0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  chi_c2          PHSP;
    Enddecay

    Decay chi_c2
    1.0000  K+  K-  pi0            PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma           PHSP;
    Enddecay

    End
DECAYCARD

decay_card_bkg_chi_c2_kk_eta = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  chi_c2          PHSP;
    Enddecay

    Decay chi_c2
    1.0000  K+  K-  eta            PHSP;
    Enddecay

    Decay eta
    1.0000  gamma  gamma           PHSP;
    Enddecay

    End
DECAYCARD

# psi' -> gamma gamma_FSR K+ K-, where gamma_FSR is a final-state radiation photon.
decay_card_bkg_fsr_kk = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  K+  K-          PHOTOS PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples ###
exMC_kk_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_kk_pi0"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_kk_pi0
  config.cross_section   = :default
end
exMC_kk_pi0.save_to_config(format: :yaml, file_path: 'exMC_psip_kk_pi0')

exMC_kk_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_kk_eta"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_kk_eta
  config.cross_section   = :default
end
exMC_kk_eta.save_to_config(format: :yaml, file_path: 'exMC_psip_kk_eta')

exMC_pi0_phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_pi0_phi"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_pi0_phi
  config.cross_section   = :default
end

exMC_eta_phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_eta_phi"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_eta_phi
  config.cross_section   = :default
end

# QED continuum MC at the two off-resonance energy points (3.773 and 3.65 GeV).
exMC_continuum_kk_pi0 = DatasetManager.create_exclusive_mc_for([psi3770_data, off_data]) do |config|
  config.sample_name   = "continuum_gammastar_to_kk_pi0"
  config.events        = 100_000
  config.decay_card    = decay_card_continuum_kk_pi0
  config.cross_section = :default
end

exMC_continuum_eta_phi = DatasetManager.create_exclusive_mc_for([psi3770_data, off_data]) do |config|
  config.sample_name   = "continuum_gammastar_to_eta_phi"
  config.events        = 100_000
  config.decay_card    = decay_card_continuum_eta_phi
  config.cross_section = :default
end

exMC_bkg_chi_c2_kk_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_chi_c2_kk_pi0_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_chi_c2_kk_pi0
  config.cross_section   = :default
end

exMC_bkg_chi_c2_kk_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_chi_c2_kk_eta_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_chi_c2_kk_eta
  config.cross_section   = :default
end

exMC_bkg_fsr_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_fsr_kk_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_fsr_kk
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# psi' -> K+ K- pi0, psi' -> K+ K- eta, psi' -> eta phi and psi' -> pi0 phi all
# decay to the same gamma gamma K+ K- final state and share this selection chain.
alg_name = "PsipToKKGammaGamma"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.686]})
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93    # track polar angle within the MDC fiducial volume
                  Vz        10.0    # within +/- 10 cm of the IP along the beam direction
                  Vr        1.0     # within 1 cm of the beam line in the transverse plane
                  nChrp     "==1"   # exactly two charged tracks with net charge zero
                  nChrn     "==1"
                  nNet      "==0"
                }
               .select_photon {
                  tdc_emc_start     0     # EMC cluster timing suppresses electronic noise and
                  tdc_emc_end      14     # energy deposits from uncorrelated events
                  angle_to_track   10.0   # photon at least 10 degrees from any charged track
                                          # (rejects bremsstrahlung showers)
                  energyThreshold_b 0.025 # E_min = 25 MeV in the barrel (|cos theta| < 0.80)
                  energyThreshold_e 0.050 # E_min = 50 MeV in the endcaps (0.86 < |cos theta| < 0.92)
                  nGam             ">=2"  # two photons for the pi0 / eta -> gamma gamma decay
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  # TOF and dE/dx probabilities are combined and each track is assigned the
                  # hypothesis with the highest confidence level; both tracks must be kaons.
                  identify :kaon, against: [:pion, :proton]
                  nkp "==1"
                  nkm "==1"
                }
               # 4C kinematic fit under the psi' -> gamma gamma K+ K- hypothesis, constrained
               # to the sum of the initial e+e- beam four-momentum (energy-momentum conservation).
               # For events with more than two photon candidates the combination with the smallest
               # chi2 is retained automatically by the DSL.
               .kinematic_fit([:gamma, :gamma, :kp, :km]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200      # loose BOSS cut; the paper requires chi2 <= 20 (applied in ROOT)
               }

alg.note(:photon_multiplicity,
         "The number of photon candidates is required to satisfy 2 <= N_gamma <= 10. The DSL " \
         "photon multiplicity keyword expresses the lower bound only (nGam \">=2\"); the upper " \
         "bound of 10 cannot be expressed and is applied at the ROOT level.")
   .note(:photon_fiducial_region,
         "Photon candidates must satisfy the BESIII EMC fiducial and shower-quality criteria: " \
         "barrel photons (|cos theta| < 0.80) with E > 25 MeV, endcap photons " \
         "(0.86 < |cos theta| < 0.92) with E > 50 MeV; showers reconstructed in the angular range " \
         "between the barrel and the endcap (0.80 < |cos theta| < 0.86) are poorly reconstructed " \
         "and are excluded. The energy thresholds are expressed through energyThreshold_b/_e; the " \
         "excluded intermediate angular band and the shower-quality selection are not expressible " \
         "through the DSL keywords.")
   .note(:kinematic_fit_chi2,
         "The 4C kinematic fit is performed under the psi' -> gamma gamma K+ K- hypothesis; " \
         "candidates with chi^2 <= 20 are retained for further analysis. This cut uses " \
         "kinematic-fit-corrected four-momenta and is therefore applied at the ROOT level, where " \
         "only the loose chi^2 < 200 requirement is imposed in BOSS.")
   .note(:pi0_selection,
         "pi0 candidates are selected from the two photons with 0.117 GeV/c^2 <= M(gamma gamma) " \
         "<= 0.147 GeV/c^2 (six times the ~5 MeV/c^2 pi0 mass resolution). To suppress " \
         "psi' -> gamma chi_c0, chi_c0 -> K+ K-, the energy of the less energetic photon is " \
         "required to be larger than 70 MeV. Background from psi' -> pi0 J/psi, " \
         "J/psi -> K+ K- is removed by requiring |M(K+ K-) - m_J/psi| >= 7 MeV/c^2. All of these " \
         "quantities are obtained from kinematic-fit-corrected four-momenta and belong to the " \
         "ROOT-level analysis. A total of 1158 psi' -> K+ K- pi0 events are selected from the data.")
   .note(:eta_selection,
         "The eta candidates are reconstructed from the two selected photons, and the yields are " \
         "determined by a fit to the M(gamma gamma) distribution (~7 MeV/c^2 resolution). The " \
         "background from psi' -> eta J/psi, J/psi -> K+ K- is suppressed by requiring " \
         "M(K+ K-) < 3.05 GeV/c^2; the background from psi' -> gamma chi_{c0,2} with " \
         "chi_{c0,2} -> pi0/eta K+ K- is suppressed by requiring the lower-energy photon to lie " \
         "outside 115-185 MeV. These selections use post-fit quantities and are ROOT-level.")
   .note(:phi_selection,
         "For psi' -> pi0 phi, the phi candidates are selected with |M(K+ K-) - m_phi| < " \
         "10 MeV/c^2; background from the initial-state-radiation process e+e- -> gamma phi is " \
         "suppressed by requiring the energy of the energetic photon to be less than 1.6 GeV. " \
         "No significant pi0 signal is observed and an upper limit of N_up = 6 events at the 90% " \
         "confidence level is set (efficiency 35.63%, Br(psi' -> pi0 phi) < 4.0 x 10^-7). These " \
         "criteria are applied to kinematic-fit-corrected quantities at the ROOT level.")
   .note(:partial_wave_analysis,
         "A partial wave analysis of psi' -> K+ K- pi0 is performed with the relativistic " \
         "covariant tensor amplitude formalism, using an unbinned maximum likelihood fit " \
         "(minimization of S = -ln L with FUMILI). The decay amplitudes A(m) = psi_mu(m) sum_i " \
         "Lambda_i U_i^mu are built from the K+, K- and pi0 four-momenta, with Breit-Wigner line " \
         "shapes for the intermediate resonances. The best solution (significance > 5 sigma) " \
         "contains K*(892)+- K-+, K2*(1430)+- K-+, K*(1680)+- K-+ and rho(1700) pi0, plus the " \
         "non-resonant K+ K- pi0 mode (P-wave K+ K- system): 224 +/- 21, 251 +/- 22, 115 +/- 20, " \
         "59 +/- 10 and 721 +/- 60 events respectively. The global goodness of fit is " \
         "chi^2_all/ndf = 147.70/126 = 1.2. This fit is not expressible in the DSL.")
   .note(:eta_phi_fit,
         "For psi' -> eta K+ K-, a two-dimensional unbinned fit to the scatter plot of " \
         "M(K+ K-) versus M(gamma gamma) is performed assuming the two variables are independent. " \
         "The fit function includes the line shapes of eta phi(1020), eta phi_3(1850), " \
         "eta phi(2170) and the non-resonant eta K+ K- decay (resonances described by " \
         "non-relativistic Breit-Wigner functions with PDG masses and widths convolved with a " \
         "detector resolution function), plus polynomial backgrounds. The yield of " \
         "psi' -> eta phi is 232 +/- 16 events, and the net signals after QED subtraction are " \
         "216 +/- 16 (eta phi) and 284 +/- 27 (eta K+ K-) events. This fit is not expressible in " \
         "the DSL.")
   .note(:background_estimation,
         "Non-pi0 background in the selected K+ K- pi0 sample, estimated from a pi0 sideband " \
         "(M(gamma gamma) in [0.079, 0.109] and [0.165, 0.195] GeV/c^2), is 43 +/- 7 events; a " \
         "low level of non-K+ K- background (3 events) arises from psi' -> pi0 pi0 J/psi, " \
         "J/psi -> mu+ mu- with muons misidentified as kaons. QED background e+e- -> gamma* -> " \
         "K+ K- pi0 is estimated from data at sqrt(s) = 3.773 GeV (195 +/- 3 events) and " \
         "3.65 GeV (195 +/- 27 events) after luminosity normalisation, and is modelled by a PWA " \
         "fit to the 3.773 GeV data set. Backgrounds are subtracted from the likelihood through " \
         "ln L = ln L_dt - sum ln L_bg, assuming no interference between signal and irreducible " \
         "QED background.")
   .note(:signal_efficiency,
         "Detection efficiencies determined from MC (weighted by the fitted intensity): " \
         "21.52% for the inclusive pi0 K+ K- channel, 20.25% for K*(892)+- K-+, 20.28% for " \
         "K2*(1430)+- K-+, 22.10% for eta K+ K- (weighted average over eta K+ K-, " \
         "eta phi_3(1850) and eta phi(2170)), 33.53% for eta phi and 35.63% for pi0 phi.")
   .note(:branching_fractions,
         "Measured branching fractions: Br(psi' -> pi0 K+ K-) = (4.07 +/- 0.16 +/- 0.26) x 10^-5; " \
         "Br(psi' -> K*(892)+ K- + c.c.) = (3.18 +/- 0.30) x 10^-5; " \
         "Br(psi' -> K2*(1430)+ K- + c.c.) = (7.12 +/- 0.62) x 10^-5 (first observation, helicity " \
         "selection rule violation); Br(psi' -> eta K+ K-) = (3.08 +/- 0.29 +/- 0.25) x 10^-5 as " \
         "((2.97 +/- 0.28) x 10^-5 in Table III); Br(psi' -> eta phi) = (3.14 +/- 0.23 +/- 0.23) x " \
         "10^-5 as ((3.08 +/- 0.29) x 10^-5 in Table III); Br(psi' -> pi0 phi) < 4.0 x 10^-7 at " \
         "the 90% confidence level.")
   .note(:systematic_uncertainties,
         "Systematic uncertainties on the branching fractions (%): photon efficiency 2, pi0 mass " \
         "cut 1.1, kaon tracking 2, PID 2, kinematic fitting 1.9 / 3.2 / 4.3 / 2.1 / 1.7 / 2.1 " \
         "(pi0 K+ K- / K*+- K-+ / K2*+- K-+ / eta K+ K- / eta phi / pi0 phi), number of psi' " \
         "decays 4, background shape 1.6 / 0.4 (eta K+ K- / eta phi), fitting range 3.6 / 0.6, " \
         "Br[KJ* -> pi0 K] 2.4, Br[P -> gamma gamma] 0.5, Br[phi -> K+ K-] 1.2, QED background 4.5, " \
         "additional states 4; totals 6.3 / 6.9 / 6.2 / 8.0 / 7.3 / 5.8. Additional PWA-related " \
         "systematics for K*(892)+- K-+ and K2*(1430)+- K-+ (Breit-Wigner form, additional " \
         "resonances, non-K+ K- pi0 background, QED background, K*(1680)/rho(1700) widths, masses " \
         "and widths of K*(892)/K2*(1430)) give totals of +5.4/-7.3% and +14.3/-5.1% respectively.")

alg.with_decay_card(decay_card_kk_pi0).apply(event_selection)

datasets = [psip_data, psip_incMC, psi3770_data, psi3770_incMC, off_data, off_incMC,
            exMC_kk_pi0, exMC_kk_eta, exMC_pi0_phi, exMC_eta_phi,
            exMC_bkg_chi_c2_kk_pi0, exMC_bkg_chi_c2_kk_eta, exMC_bkg_fsr_kk] +
           exMC_continuum_kk_pi0 + exMC_continuum_eta_phi

root_files = alg.execute_on(datasets)
