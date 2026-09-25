# =============================================================================
# BESIII: observation of e+e- -> omega chi_c1,2 at sqrt(s) > 4.4 GeV
# arXiv:1511.08564v2
#
# Data: five energy points in the range 4.416 <= sqrt(s) <= 4.599 GeV
# (4.416, 4.467, 4.527, 4.574 and 4.599 GeV).
#
# Signal channels:
#   e+e- -> omega chi_c1,2, chi_c1,2 -> gamma J/psi, J/psi -> l+ l- (l = e, mu),
#           omega -> pi+ pi- pi0                 -> final state gamma pi+ pi- pi0 l+ l-
#   e+e- -> omega chi_c0, chi_c0 -> pi+ pi-  or  K+ K-,  omega -> pi+ pi- pi0
#
# Signal candidates must have exactly four charged tracks with zero net charge,
# a pi0 candidate and a photon. A 5C kinematic fit constrains the total
# four-momentum of the final state to the initial four-momentum of the colliding
# beams and constrains the invariant mass of the two photons from the pi0 to the
# nominal pi0 mass; the published requirement is chi2_5C < 60.
#
# BOSS part only: decay cards, exclusive MC and the event selection up to and
# including the final kinematic fit. The signal-region counting (omega, J/psi
# and chi_c1,2 mass windows), the background estimate from the non-omega
# sidebands, the fits to M(gamma J/psi) and the Born cross sections with their
# upper limits are ROOT-level and are captured in the notes below.
# =============================================================================

### Dataset description ###
# The five centre-of-mass energies of this analysis; the 4.416 GeV point is the
# high-luminosity 4420 sample (1043.9 pb^-1 in round-07, combined with the
# 46.8 pb^-1 of round-06 for the quoted 1074 pb^-1).
data_4416 = DatasetManager.real_data.find("703_4420")   # sqrt(s) = 4415.58 MeV, 1043.9 pb^-1
data_4467 = DatasetManager.real_data.find("703_4470")   # sqrt(s) = 4467.06 MeV,  111.09 pb^-1
data_4527 = DatasetManager.real_data.find("703_4530")   # sqrt(s) = 4527.14 MeV,  112.12 pb^-1
data_4574 = DatasetManager.real_data.find("703_4575")   # sqrt(s) = 4574.50 MeV,   48.93 pb^-1
data_4599 = DatasetManager.real_data.find("703_4600")   # sqrt(s) = 4599.53 MeV,  586.9 pb^-1

incMC_4416 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4467 = DatasetManager.inclusive_mc.find("703_4470")
incMC_4527 = DatasetManager.inclusive_mc.find("703_4530")
incMC_4574 = DatasetManager.inclusive_mc.find("703_4575")
incMC_4599 = DatasetManager.inclusive_mc.find("703_4600")

data_points = [data_4416, data_4467, data_4527, data_4574, data_4599]
inc_mc_points = [incMC_4416, incMC_4467, incMC_4527, incMC_4574, incMC_4599]

### Decay cards ###
# The Born cross section is measured as a function of sqrt(s) at five scan
# points; no ConExc mode covers the omega chi_cJ final states, so the BESIII
# KKMC convention (psi(4260) as top mother, PHSP generators) is used and the
# ISR / vacuum-polarisation corrections are taken from the dedicated QED
# calculation quoted in the paper (see the note :radiative_correction).

# Channel I: e+e- -> omega chi_c1,2 with chi_c1,2 -> gamma J/psi, J/psi -> e+ e-
decay_card_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 omega chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Same final state with chi_c2 and J/psi -> e+ e- (used for the chi_c2 signal)
decay_card_ee_c2 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 omega chi_c2 PHSP;
    Enddecay

    Decay chi_c2
    1.0000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Channel I with J/psi -> mu+ mu- (the two lepton channels are summed in the
# cross-section denominator through eps_e B_e + eps_mu B_mu)
decay_card_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.0000 omega chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu- PHOTOS VLL;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Channel II: e+e- -> omega chi_c0 with chi_c0 -> pi+ pi-
decay_card_c0_pipi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 omega chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Channel II: e+e- -> omega chi_c0 with chi_c0 -> K+ K-
decay_card_c0_kk = <<~DECAYCARD
    Decay psi(4260)
    1.0000 omega chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.0000 K+ K- PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive signal MC — the same card is run at every scan point ###
exMC_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "OmegaChiC1_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_ee
  config.cross_section = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "OmegaChiC1_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_mumu
  config.cross_section = :default
end

exMC_c2_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "OmegaChiC2_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_ee_c2
  config.cross_section = :default
end

exMC_c0_pipi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "OmegaChiC0_pipi"
  config.events        = 100_000
  config.decay_card    = decay_card_c0_pipi
  config.cross_section = :default
end

exMC_c0_kk = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "OmegaChiC0_KK"
  config.events        = 100_000
  config.decay_card    = decay_card_c0_kk
  config.cross_section = :default
end

### Event selection (BOSS) ###
# The two J/psi lepton channels share the identical final state and selection, so
# a single Algorithm is used for chi_c1 (both l = e and l = mu are handled by the
# lepton identification, and the two channels are summed through the published
# eps_e B_e + eps_mu B_mu denominator). The chi_c2 signal shares the same final
# state and selection and therefore uses the same algorithm; it is separated from
# chi_c1 only by the M(gamma J/psi) window applied in the ROOT analysis.
# The chi_c0 channels (chi_c0 -> pi+ pi- and chi_c0 -> K+ K-) have different
# final states and PID hypotheses and get their own algorithms (Rule T1).
#
# ECMS is intentionally not declared: the same algorithm runs at the five
# centre-of-mass energies and the per-job energy is injected through jobOptions.

# ---------------------------------------------------------------------------
# Channel I — e+e- -> omega chi_c1,2, chi_c1,2 -> gamma J/psi, J/psi -> l+ l-,
#             omega -> pi+ pi- pi0
# ---------------------------------------------------------------------------
alg_name_wc = "OmegaChiC1LL"
alg_wc = Algorithm.new(alg_name_wc)
alg_wc.set_header(["#{alg_name_wc}Alg/#{alg_name_wc}.h"])
      .set_alias({ "std::vector<double>" => "Vdouble" })

selection_wc = Selection.new
selection_wc.select_track {
              cos_theta   0.93    # |cos(theta)| < 0.93 for the MDC tracks
              Vz          10.0    # |Vz| < 10 cm along the beam direction
              Vr          1.0     # Vr < 1 cm in the plane perpendicular to the beam
              nChrp       "==2"   # exactly two positively charged tracks
              nChrn       "==2"   # exactly two negatively charged tracks
              nNet        "==0"   # zero net charge: exactly four tracks
            }
           .select_photon {
              tdc_emc_start     0      # EMC timing window, in units of 700 ns
              tdc_emc_end       14
              angle_to_track    10.0   # photon at least 10 deg from any charged track
              energyThreshold_b 0.025  # 25 MeV in the barrel
              energyThreshold_e 0.050  # 50 MeV in the end-cap
              nGam              ">=3"  # two photons from pi0 plus the chi_c1,2 photon
            }
           .pid(method: :probability) {
              prob_cut 0.001
              # The leptons from J/psi have momentum above 1 GeV/c; the EMC energy
              # deposit separates electrons from muons
              identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                             treat_as_electron_if_energy_above: 0.6
              identify :pion, against: [:kaon, :proton]  # the two pions from omega
              nlp "==1"    # one lepton from J/psi -> l+ l-
              nlm "==1"
              npip "==1"   # one pion from omega -> pi+ pi- pi0
              npim "==1"
            }
           # 5C kinematic fit: total four-momentum constrained to the initial
           # beams (4C) plus the invariant mass of the photon pair from the pi0
           # constrained to the nominal pi0 mass (1C). The chi_c1,2 photon is the
           # remaining photon; the fit selects the photon-pair combination with
           # the smallest chi2. Loose BOSS cut, the paper's chi2_5C < 60 is
           # applied in ROOT.
           .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim, :lp, :lm]) {
              nominal
              constrain_four_momentum
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
              chi2_cut 200
            }

alg_wc.note(:chi2_selection,
            "The published analysis performs a five-constraint (5C) kinematic fit " \
            "constraining the total four-momentum of the final state " \
            "gamma pi+ pi- pi0 l+ l- to the initial four-momentum of the colliding " \
            "beams and the invariant mass of the two photons from the pi0 to the " \
            "nominal pi0 mass, and requires chi2_5C < 60. The nominal BOSS pass " \
            "applies the loose default chi2_cut 200; the tight chi2_5C < 60 value " \
            "is applied at the ROOT stage.")
      .note(:signal_regions,
            "Signal regions are defined on invariants computed after the 5C fit " \
            "and are therefore applied in the ROOT analysis: " \
            "0.75 <= M(pi+ pi- pi0) <= 0.81 GeV/c^2 for the omega, " \
            "3.08 <= M(l+ l-) <= 3.12 GeV/c^2 for the J/psi, " \
            "[3.49, 3.53] GeV/c^2 for chi_c1 and [3.54, 3.58] GeV/c^2 for chi_c2 " \
            "(the J/psi mass resolution is 8 MeV/c^2 in the MC simulation); " \
            "[3.39, 3.47] GeV/c^2 is used as the chi_c1,2 sideband.")
      .note(:background_veto,
            "The dominant background is e+e- -> pi+ pi- pi0 chi_c1,2 with a " \
            "non-resonant pi+ pi- pi0, which peaks in the chi_c1,2 signal region. " \
            "It is estimated from the three non-omega regions (boxes A, B, C) of the " \
            "M(pi+ pi- pi0) versus M(gamma l+ l-) plane using " \
            "n_bkg = f (n_A + n_B - 0.5 n_C), with the normalisation factor f " \
            "obtained from phase-space pi+ pi- pi0 chi_c1,2 simulation. Other " \
            "backgrounds (eta' J/psi with eta' -> gamma omega, pi+ pi- psi(2S) with " \
            "psi(2S) -> pi0 pi0 J/psi or gamma chi_c1,2, and pi0 pi0 psi(2S) with " \
            "psi(2S) -> pi+ pi- J/psi) do not peak in the signal regions and are " \
            "negligible. These estimates are ROOT-level procedures.")
      .note(:radiative_correction,
            "The Born cross section is extracted as sigma_B = N_sig / " \
            "[ L (1+delta) (1/|1-Pi|^2) (eps_e B_e + eps_mu B_mu) B_1 ], where " \
            "(1+delta) is the radiative correction factor obtained from a QED " \
            "calculation that is iterated with the measured cross section until it " \
            "converges, and 1/|1-Pi|^2 is the vacuum-polarisation factor taken from " \
            "a QED calculation with 0.5% accuracy. These factors are provided " \
            "externally to the BOSS selection.")
      .note(:helix_correction,
            "The helix parameters of the simulated charged tracks are corrected so " \
            "that the MC momentum spectra match the data; the correction factors " \
            "for pi, e and mu are obtained from the control samples " \
            "e+e- -> pi+ pi- J/psi, J/psi -> e+ e- and J/psi -> mu+ mu-. The " \
            "difference in the MC efficiency with and without the correction is " \
            "taken as a systematic uncertainty.")
      .note(:angular_distribution,
            "The omega helicity angular distribution is generated as " \
            "1 +/- cos^2(theta_1) (theta_1 the polar angle of the omega in the " \
            "e+e- rest frame, z along the electron beam) and the chi_c1,2 photon " \
            "helicity angular distribution as 1 +/- cos^2(theta_2) (theta_2 the " \
            "polar angle of the photon in the chi_c1,2 rest frame with z along the " \
            "omega direction) instead of the PHSP model; the maximum change in the " \
            "MC efficiencies is taken as the systematic uncertainty.")
      .note(:upper_limit_method,
            "At the energy points where the signals are not significant the upper " \
            "limits on the Born cross section at 90% C.L. are computed with the " \
            "frequentist method of unbounded profile likelihood (the trolke package " \
            "in ROOT), assuming a Poisson-distributed background and a Gaussian " \
            "efficiency uncertainty, using the full denominator of the " \
            "cross-section formula as effective efficiency.")
      .note(:systematic_uncertainties,
            "Relative systematic uncertainties: luminosity 1.0%, tracking 4.0% " \
            "(1.0% per track), photon 1.0% per photon (from J/psi -> rho0 pi0), " \
            "J/psi mass window 1.6%, kinematic fit (1.3, 2.0, 2.1)% for " \
            "chi_c0, chi_c1, chi_c2, angular distribution and line shape as " \
            "tabulated, branching fractions B_e, B_mu and B_1 from the world " \
            "average, and the pi+ pi- / K+ K- cross feed for the chi_c0 channel.")

alg_wc.with_decay_card(decay_card_ee).apply(selection_wc)

# ---------------------------------------------------------------------------
# Channel II — e+e- -> omega chi_c0 with chi_c0 -> pi+ pi-
# ---------------------------------------------------------------------------
alg_name_c0pipi = "OmegaChiC0PiPi"
alg_c0pipi = Algorithm.new(alg_name_c0pipi)
alg_c0pipi.set_header(["#{alg_name_c0pipi}Alg/#{alg_name_c0pipi}.h"])
         .set_alias({ "std::vector<double>" => "Vdouble" })

selection_c0pipi = Selection.new
selection_c0pipi.select_track {
                  cos_theta   0.93
                  Vz          10.0
                  Vr          1.0
                  nChrp       "==2"   # two pi+ (one from omega, one from chi_c0)
                  nChrn       "==2"   # two pi- (one from omega, one from chi_c0)
                  nNet        "==0"
                }
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam              ">=2"   # the two photons from the pi0 of omega
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :pion, against: [:kaon, :proton]
                  npip "==2"
                  npim "==2"
                }
               # 5C kinematic fit: total four-momentum to the initial beams (4C)
               # plus the pi0 mass constraint on the photon pair (1C)
               .kinematic_fit([:gamma, :gamma, :pip, :pip, :pim, :pim]) {
                  nominal
                  constrain_four_momentum
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 200
                }

alg_c0pipi.note(:selection_from_ref19,
                "The e+e- -> omega chi_c0 selection follows the same criteria as " \
                "Ref. [19] (the BESIII observation of e+e- -> omega chi_c0 around " \
                "4.23 GeV): exactly four charged tracks with zero net charge, at " \
                "least two photons forming a pi0 candidate, a 5C kinematic fit with " \
                "chi2_5C < 60 (loose default chi2_cut 200 applied in BOSS). The " \
                "omega and chi_c0 signal windows are applied in ROOT.")
      .note(:signal_regions,
                "Signal regions applied in the ROOT analysis: " \
                "0.75 <= M(pi+ pi- pi0) <= 0.81 GeV/c^2 for the omega and the " \
                "chi_c0 window around the nominal chi_c0 mass in the pi+ pi- " \
                "invariant mass.")
      .note(:cross_feed,
                "The uncertainty due to the cross feed between the chi_c0 -> pi+ pi- " \
                "and chi_c0 -> K+ K- modes is estimated from the signal MC samples.")
      .note(:radiative_correction,
                "The Born cross section uses the same denominator structure as the " \
                "omega chi_c1,2 channels, with " \
                "D = (1+delta) (1/|1-Pi|^2) [eps_pi B(chi_c0 -> pi+ pi-) + " \
                "eps_K B(chi_c0 -> K+ K-)] B(omega -> pi+ pi- pi0) B(pi0 -> gamma gamma).")

alg_c0pipi.with_decay_card(decay_card_c0_pipi).apply(selection_c0pipi)

# ---------------------------------------------------------------------------
# Channel III — e+e- -> omega chi_c0 with chi_c0 -> K+ K-
# ---------------------------------------------------------------------------
alg_name_c0kk = "OmegaChiC0KK"
alg_c0kk = Algorithm.new(alg_name_c0kk)
alg_c0kk.set_header(["#{alg_name_c0kk}Alg/#{alg_name_c0kk}.h"])
       .set_alias({ "std::vector<double>" => "Vdouble" })

selection_c0kk = Selection.new
selection_c0kk.select_track {
                cos_theta   0.93
                Vz          10.0
                Vr          1.0
                nChrp       "==2"   # pi+ from omega and K+ from chi_c0
                nChrn       "==2"   # pi- from omega and K- from chi_c0
                nNet        "==0"
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
                prob_cut 0.001
                identify :kaon, against: [:pion, :proton]
                identify :pion, against: [:kaon, :proton]
                nkp "==1"    # the K+ from chi_c0 -> K+ K-
                nkm "==1"
                npip "==1"   # the pi+ from omega -> pi+ pi- pi0
                npim "==1"
              }
             .kinematic_fit([:gamma, :gamma, :kp, :km, :pip, :pim]) {
                nominal
                constrain_four_momentum
                invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                chi2_cut 200
              }

alg_c0kk.note(:selection_from_ref19,
              "Same selection strategy as the chi_c0 -> pi+ pi- channel, with the " \
              "charged pions of the chi_c0 decay replaced by a K+ K- pair " \
              "identified by the combined dE/dx and TOF PID. The 5C fit uses the " \
              "kaon mass hypothesis for those two tracks.")
        .note(:signal_regions,
              "Signal regions applied in the ROOT analysis: " \
              "0.75 <= M(pi+ pi- pi0) <= 0.81 GeV/c^2 for the omega and the " \
              "chi_c0 window in the K+ K- invariant mass.")
        .note(:cross_feed,
              "The pi+ pi- / K+ K- cross feed is estimated from the signal MC " \
              "samples and is one of the systematic uncertainties of the " \
              "omega chi_c0 measurement.")
        .note(:upper_limit_method,
              "The omega chi_c0 signals are not significant at any energy point, so " \
              "only upper limits on the Born cross section at 90% C.L. are derived " \
              "using the frequentist unbounded-profile-likelihood method of the " \
              "trolke package.")

alg_c0kk.with_decay_card(decay_card_c0_kk).apply(selection_c0kk)

### Execution ###
# Every algorithm runs over the same five energy points; the signal channels
# share the lepton-flavour exclusive samples.
root_files_wc = alg_wc.execute_on(data_points + inc_mc_points +
                                  exMC_ee + exMC_mumu + exMC_c2_ee)
root_files_c0pipi = alg_c0pipi.execute_on(data_points + inc_mc_points + exMC_c0_pipi)
root_files_c0kk = alg_c0kk.execute_on(data_points + inc_mc_points + exMC_c0_kk)
