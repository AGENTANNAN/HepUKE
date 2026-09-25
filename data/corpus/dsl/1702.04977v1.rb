# =============================================================================
# 1702.04977v1
# Measurement of the integrated luminosities of the BESIII data taken between
# 2.2324 and 4.5900 GeV (131 c.m. energy points) from the two QED processes
# e+e- -> (gamma) e+e-  (large-angle Bhabha) and e+e- -> gamma gamma (diphoton).
#
# BOSS-side spec: dataset preparation + event selection up to the final 4C
# kinematic fit. The event counting, the luminosity extraction
# L = (N_obs - N_bkg) / (sigma_QED * eps * eps_trig) and the background
# sideband subtraction are ROOT-level and are not part of this spec.
#
# The analysis is an ordinary (non-tag) analysis and the two QED processes have
# different final states (two charged tracks vs. two photons) and different
# selection criteria, so they get two separate Algorithm objects (Rule T1).
# Multi-energy scan: ECMS is deliberately NOT declared as an algorithm constant;
# every job takes its c.m. energy from the related dataset.
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets — data points of the 2012, 2013 and 2014 runs
### ---------------------------------------------------------------------------
# The measurement uses 131 c.m. energy points: 4 points of the 2012 run
# (2.2324, 2.4000, 2.8000, 3.4000 GeV), 104 points from 3.8500 to 4.5900 GeV
# (2013-2014 runs), 15 points near the J/psi production threshold, 4 points of
# the tau mass measurement and 4 points for charmonium studies. The sample table
# does not resolve every one of these points individually, so the available
# samples inside the c.m. energy window of the paper are collected with a range
# query; the selection below is energy independent, so the identical chain is
# applied at every point.
data_points  = DatasetManager.real_data.where(cms_energy: {value: 2232..4590})
incMC_points = DatasetManager.inclusive_mc.where(cms_energy: {value: 2232..4590})

### ---------------------------------------------------------------------------
### Decay cards
### ---------------------------------------------------------------------------
# Signal: large-angle Bhabha scattering e+e- -> (gamma) e+e-; the radiative
# photon is emitted as final-state radiation and is not required to be
# reconstructed. The real generator is Babayaga v3.5 (see the generator note);
# the declared card is the standard BOSS approximation.
decay_card_bhabha = <<~DECAYCARD
  Decay psi(4260)
  1.0000 e+ e- gamma           PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Signal: e+e- -> gamma gamma (diphoton luminosity measurement).
decay_card_gammagamma = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma gamma           PHSP;
  Enddecay

  End
DECAYCARD

# Background: e+e- -> tau+ tau- (generated with KKMC in the paper); the tau
# leptons can fake the two-track topology of the Bhabha selection.
decay_card_tautau = <<~DECAYCARD
  Decay psi(4260)
  1.0000 tau+ tau-             PHSP;
  Enddecay

  End
DECAYCARD

# Background: continuum e+e- -> hadrons (generated with LUARLW in the paper);
# the hadronic contamination is the dominant background of the diphoton
# selection and is estimated from the Delta phi sideband.
decay_card_hadrons = <<~DECAYCARD
  Decay psi(4260)
  1.0000 generic generic generic GENERIC;
  Enddecay

  End
DECAYCARD

### ---------------------------------------------------------------------------
### Exclusive MC — one signal/background sample per energy point
### ---------------------------------------------------------------------------
# "The e+e- -> (gamma) e+e-, gamma gamma and (gamma) mu+ mu- events are
# simulated with the generator Babayaga v3.5."
exMC_bhabha = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_bhabha_babayaga"    # auto-suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_bhabha
  config.cross_section = :default
end

exMC_gammagamma = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_gammagamma_babayaga"
  config.events        = 200_000
  config.decay_card    = decay_card_gammagamma
  config.cross_section = :default
end

exMC_tautau = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "bkg_tautau"
  config.events        = 100_000
  config.decay_card    = decay_card_tautau
  config.cross_section = :default
end

exMC_hadrons = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "bkg_continuum_hadrons"
  config.events        = 100_000
  config.decay_card    = decay_card_hadrons
  config.cross_section = :default
end

### ---------------------------------------------------------------------------
### Algorithm 1 — e+e- -> (gamma) e+e- (large-angle Bhabha)
### ---------------------------------------------------------------------------
alg_name_bhabha = "BhabhaLumi"
alg_bhabha = Algorithm.new(alg_name_bhabha)
alg_bhabha.set_header(["#{alg_name_bhabha}Alg/#{alg_name_bhabha}.h"])
# Multi-energy scan: ECMS is not set here, it is injected per energy point.

event_selection_bhabha = Selection.new
# --- Two good charged tracks with opposite charge ---------------------------
# Each track within +/-10 cm of the IP along the beam direction and 1 cm in the
# plane perpendicular to the beam, with |cos(theta)| < 0.8 in the MDC.
event_selection_bhabha
  .select_track do
    cos_theta 0.8      # |cos(theta)| < 0.8 (MDC polar angle)
    Vz        10.0     # |Vz| < 10 cm along the beam axis
    Vr        1.0      # Vr < 1 cm transverse to the beam axis
    nChrp     "==1"    # exactly one positively charged track
    nChrn     "==1"    # exactly one negatively charged track
    nNet      "==0"    # net charge zero
  end
  # --- EMC deposited-energy requirement (e/mu separation, mu+mu- suppression) --
  # The paper requires E_e+-(EMC) > 0.65 * E_beam.  The large-angle Bhabha
  # tracks carry p = E_beam, so the requirement is encoded energy-independently
  # as the EMC-energy-to-momentum ratio E_EMC/p > 0.65, valid at every c.m.
  # energy of the scan.
  .for_each(:charged) do
    define(:eraw_over_p) { eraw / p }
    where { eraw_over_p < 0.65 }
    remove
  end
  # --- No particle identification ---------------------------------------------
  # "Without applying further particle identification, the tracks are assigned
  # as electron and positron depending on their charges."
  .assign({chrgp: :ep, chrgn: :em})
  # --- Nominal 4C fit of the two leptons to the c.m. four-momentum -------------
  # The published luminosity is obtained by event counting and does not use a
  # kinematic fit; the fit is kept as the DSL end-point only.
  .kinematic_fit([:ep, :em]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_bhabha
  .note(:generator_model,
        "The Bhabha signal and the (gamma) e+e-, gamma gamma and (gamma) mu+mu- " \
        "background samples are generated with Babayaga v3.5, which models the " \
        "running of the electromagnetic coupling constant and final-state " \
        "radiation; the generator uncertainty is 0.5% for e+e- -> (gamma) e+e- and " \
        "1.0% for e+e- -> gamma gamma. Babayaga is a dedicated QED generator and is " \
        "not expressible as an EvtGen decay card; the declared card " \
        "(psi(4260) -> e+ e- gamma with PHOTOS VLL, the standard BOSS " \
        "approximation) is the closest available stand-in.")
  .note(:back_to_back_requirement,
        "The two charged tracks must be back-to-back in the c.m. system: " \
        "|Delta theta_e+-| = |theta_1 + theta_2 - 180 deg| < 10.0 deg and " \
        "|Delta phi_e+-| = ||phi_1 - phi_2| - 180 deg| < 5.0 deg. These are " \
        "pair-level variables built from both tracks and have no DSL construct; " \
        "they are applied as offline selections on the NTuple.")
  .note(:emc_energy_requirement,
        "The paper requires the EMC deposited energy of both the electron and the " \
        "positron to exceed 0.65 * E_beam, which suppresses the e+e- -> mu+mu- and " \
        "hadronic backgrounds. Because the large-angle Bhabha tracks satisfy " \
        "p = E_beam, the BOSS-level filter uses the equivalent, energy-independent " \
        "ratio E_EMC/p > 0.65; the exact E_EMC/E_beam form is kept in the ROOT " \
        "analysis.")
  .note(:no_kinematic_fit_in_paper,
        "The published measurement does NOT apply a kinematic fit: the number of " \
        "observed QED events N_QED^obs is obtained by counting the selected events " \
        "and the luminosity follows from " \
        "L = (N_QED^obs - N_QED^bkg) / (sigma_QED * eps_QED * eps_QED^trig). The " \
        "nominal 4C fit declared here is only the DSL end-point.")
  .note(:efficiency_curve,
        "The detection efficiency eps_QED is obtained from the signal MC by " \
        "applying the same selection as in data, at each c.m. energy individually. " \
        "The MC-statistics uncertainty is calculated as " \
        "1/sqrt(N) * sqrt((1 - eps)/eps), where N is the number of signal MC events " \
        "and eps the detection efficiency, giving 0.17% for the (gamma) e+e- process. " \
        "The trigger efficiency for barrel e+e- -> (gamma) e+e- events is 100% with " \
        "an uncertainty below 0.1%.")
  .note(:tracking_efficiency,
        "The tracking-efficiency systematic uncertainty is studied with a Bhabha " \
        "sample selected with EMC information only (no MDC tracking): at least two " \
        "EMC clusters with deposited energy above 0.65 * E_beam and " \
        "|cos(theta_EMC)| < 0.8, with the azimuthal difference " \
        "|Delta phi_e+-| in [-40 deg, -5 deg] or [5 deg, 40 deg] to remove " \
        "e+e- -> gamma gamma events (the two clusters are not back-to-back in the " \
        "xy-plane because the showers are bent by the magnetic field). The MDC " \
        "information is then applied and the ratio of surviving events is the " \
        "tracking efficiency; the data/MC difference of 0.41% is the systematic " \
        "uncertainty. The cluster-pair conditions have no DSL construct.")
  .note(:background_veto,
        "The background of the (gamma) e+e- selection (e+e- -> tau+tau-, " \
        "e+e- -> hadrons, e+e- -> e+e- + X with X hadrons or leptons) is estimated " \
        "by applying the identical requirements to the background MC samples " \
        "(tau+tau- with KKMC, hadrons with LUARLW, e+e- + X with BesTwogam) and " \
        "normalising to the data; the resulting background level after the " \
        "selection is 1e-5, and its uncertainty is neglected.")
  .note(:energy_points,
        "The measurement is performed for 131 c.m. energy points between 2.2324 " \
        "and 4.5900 GeV, individually. The declared dataset range query collects " \
        "the available samples inside that window; the per-point ECMS, the " \
        "integrated luminosities and the measured luminosities of Table III are " \
        "resolved in the analysis stage.")

### ---------------------------------------------------------------------------
### Algorithm 2 — e+e- -> gamma gamma (diphoton)
### ---------------------------------------------------------------------------
alg_name_gg = "GammaGammaLumi"
alg_gammagamma = Algorithm.new(alg_name_gg)
alg_gammagamma.set_header(["#{alg_name_gg}Alg/#{alg_name_gg}.h"])
# Multi-energy scan: ECMS is not set here, it is injected per energy point.

event_selection_gg = Selection.new
# --- Charged track veto ------------------------------------------------------
# "To select e+e- -> gamma gamma events, the number of good charged tracks is
# required to be zero" (the diphoton measurement uses EMC information only).
event_selection_gg
  .select_track do
    cos_theta 0.8      # |cos(theta)| < 0.8 (MDC polar angle)
    Vz        10.0     # |Vz| < 10 cm along the beam axis
    Vr        1.0      # Vr < 1 cm transverse to the beam axis
    nChrp     "==0"    # no positively charged track
    nChrn     "==0"    # no negatively charged track
    nNet      "==0"    # net charge zero
  end
  # --- Two neutral clusters ----------------------------------------------------
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14       # EMC timing window 0 < t < 700 ns
    energyThreshold_b 0.025    # E > 25 MeV in the barrel EMC (|cos(theta)| < 0.8)
    energyThreshold_e 0.050    # E > 50 MeV in the endcap EMC (0.86 < |cos(theta)| < 0.92)
    angle_to_track    10.0
    nGam              ">=2"    # the two most energetic clusters are the gamma gamma pair
  end
  # --- Nominal 4C fit of the two photons to the c.m. four-momentum -------------
  # The published luminosity is obtained by event counting and does not use a
  # kinematic fit; the fit is kept as the DSL end-point only.
  .kinematic_fit([:gamma, :gamma]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg_gammagamma
  .note(:generator_model,
        "The e+e- -> gamma gamma signal and the accompanying QED samples are " \
        "generated with Babayaga v3.5 (generator uncertainty 1.0% for " \
        "e+e- -> gamma gamma); the declared card (psi(4260) -> gamma gamma PHSP) " \
        "is the standard BOSS approximation of that generator.")
  .note(:photon_energy_windows,
        "The two selected neutral clusters must satisfy |cos(theta)| < 0.8, which " \
        "restricts them to the barrel EMC, and the deposited energy must satisfy " \
        "0.7 < E_gamma / E_beam < 1.16. Both are energy-dependent, cluster-level " \
        "conditions relative to the beam energy of each of the 131 points and are " \
        "applied on the NTuple, not in BOSS.")
  .note(:back_to_back_requirement,
        "The two photon candidates must be back to back in the azimuthal plane: " \
        "|Delta phi_gamma| = |phi_gamma1 - phi_gamma2| < 2.5 deg, where " \
        "phi_gamma1/2 are the azimuthal angles of the two photons. This is a " \
        "pair-level variable built from both clusters and has no DSL construct.")
  .note(:sideband_background,
        "The dominant background of the diphoton selection is the hadronic " \
        "contamination. The number of background events is extracted from the " \
        "Delta phi_gamma sideband region 2.5 deg < |Delta phi_gamma| < 5.0 deg, " \
        "whose distribution is verified to be flat with the background MC samples. " \
        "The background rates in data and MC are (1.53 +/- 0.03)% and " \
        "(1.31 +/- 0.04)%, respectively, and their difference 0.23% is taken as the " \
        "associated systematic uncertainty.")
  .note(:no_kinematic_fit_in_paper,
        "The published cross check does not apply a kinematic fit either; the " \
        "diphoton yield is obtained by event counting after the selection (all " \
        "energy points outside the J/psi region; for points from 3.0930 to " \
        "3.1200 GeV only this diphoton measurement is used, because " \
        "e+e- -> (gamma) J/psi events cannot be distinguished from " \
        "e+e- -> (gamma) e+e- experimentally). The nominal 4C fit declared here is " \
        "only the DSL end-point.")
  .note(:efficiency_curve,
        "The detection efficiency eps_QED is obtained from the diphoton signal MC " \
        "at each c.m. energy separately; the MC-statistics uncertainty is " \
        "1/sqrt(N) * sqrt((1 - eps)/eps) and amounts to 0.15% for this process. " \
        "The trigger efficiency for e+e- -> gamma gamma events is 100% with an " \
        "uncertainty below 0.1%.")
  .note(:background_veto,
        "The diphoton background is dominated by continuum hadronic events, " \
        "generated with LUARLW in the paper and approximated here by the standard " \
        "inclusive generic card; e+e- -> e+e- + X (BesTwogam) and " \
        "e+e- -> tau+tau- (KKMC) events are also generated for the background " \
        "study. The selection conditions used to define the sideband are applied " \
        "in the ROOT analysis.")
  .note(:systematic_uncertainty,
        "The luminosity uncertainty is obtained by re-measuring with the selection " \
        "conditions varied one at a time: |cos(theta)| < 0.8, |Delta theta_e+-| " \
        "< 10 deg, |Delta phi_e+-| < 5 deg, |Delta phi_gamma| < 2.5 deg, " \
        "E_e+-/E_beam > 0.65 and 0.7 < E_gamma/E_beam < 1.16. The further sources " \
        "are the tracking efficiency (0.41%), the cluster reconstruction " \
        "efficiency (0.05% per cluster), the c.m. energy (2 MeV, re-estimated with " \
        "an alternative MC sample 2 MeV above nominal), the MC statistics, the " \
        "background estimation, the trigger efficiency (< 0.1%) and the generator. " \
        "The total uncertainty is 0.7% for e+e- -> (gamma) e+e- and 1.1% for " \
        "e+e- -> gamma gamma at sqrt(s) = 2.2324 GeV; the tracking, cluster, " \
        "trigger and generator uncertainties are common to all energy points, the " \
        "others are determined per point.")

alg_bhabha.with_decay_card(decay_card_bhabha).apply(event_selection_bhabha)
alg_gammagamma.with_decay_card(decay_card_gammagamma).apply(event_selection_gg)

### ---------------------------------------------------------------------------
### Execution
### ---------------------------------------------------------------------------
alg_bhabha.execute_on(data_points + incMC_points +
                      exMC_bhabha + exMC_tautau + exMC_hadrons)

alg_gammagamma.execute_on(data_points + incMC_points +
                          exMC_gammagamma + exMC_tautau + exMC_hadrons)
