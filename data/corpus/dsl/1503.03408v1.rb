# =============================================================================
# BESIII precision measurement of the integrated luminosity of the data taken
# between 3.810 GeV and 4.600 GeV (arXiv:1503.03408v1)
#   --  BOSS part (dataset preparation + event selection)
#
# The time-integrated luminosity of ~5 fb^-1 of e+e- data collected at 21 c.m.
# energies (3.810-4.600 GeV) is measured with the large-angle Bhabha process
# e+e- -> (gamma) e+e-, and cross checked with the di-gamma process
# e+e- -> gamma gamma.  Two independent selections are therefore defined, one
# per process; each runs over every energy point.
#
# This is a multi-energy (scan) measurement, so ECMS is deliberately NOT
# declared as an algorithm constant here: every job takes its c.m. energy from
# the related dataset (the C++ property default acts as a placeholder).
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets (21 data samples between 3.810 and 4.600 GeV)
### ---------------------------------------------------------------------------
# The dataset table lists 21 samples but does not resolve the two separately
# acquired samples at 4.230, 4.260 and 4.420 GeV (Table I superscripts 1 and 2),
# so 18 distinct sample names are used below; the pairs are merged.
sample_names = %w[
  703_3810 703_3900 703_4009 703_4090 703_4190 703_4210 703_4220 703_4230
  703_4245 703_4260 703_4310 703_4360 703_4390 703_4420 703_4470 703_4530
  703_4575 703_4600
]

data_points = sample_names.map { |name| DatasetManager.real_data.find(name) }

# Inclusive MC corresponding to 500 pb^-1 at 4.260 GeV (QED processes, continuum
# hadrons and ISR to J/psi and psi(3686)); used for background studies and for
# the optimisation of the selection criteria (Sec. III).
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

### ---------------------------------------------------------------------------
### Decay cards
### ---------------------------------------------------------------------------
# Signal: large-angle Bhabha scattering e+e- -> (gamma) e+e-.  The visible final
# state is an electron-positron pair; the radiative photon is generated as FSR
# (PHOTOS) as required by the BABAYAGA3.5 settings (RunningAlpha 1, FSR switch 1).
# 'Particle vpho' is not used: the card relies on the per-dataset c.m. energy
# injected by execute_on, so the same card is valid at every energy point.
decay_card_bhabha = <<~DECAYCARD
  Decay psi(4260)
  1.0000 e+ e- gamma           PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Cross-check channel: e+e- -> gamma gamma (di-gamma luminosity measurement).
decay_card_gammagamma = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma gamma           PHSP;
  Enddecay

  End
DECAYCARD

### ---------------------------------------------------------------------------
### Exclusive MC: one Bhabha and one di-gamma sample at every energy point
### ---------------------------------------------------------------------------
# "At each energy point, one million Bhabha events were generated using the
# BABAYAGA3.5 generator" (Sec. III).  The same production is used for the
# di-gamma cross check.
exMC_bhabha = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name     = "sig_bhabha_babayaga"
  config.events          = 1_000_000
  config.decay_card      = decay_card_bhabha
  config.cross_section   = :default
end

exMC_gammagamma = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name     = "ck_gammagamma"
  config.events          = 1_000_000
  config.decay_card      = decay_card_gammagamma
  config.cross_section   = :default
end

### ---------------------------------------------------------------------------
### Algorithm 1 -- nominal measurement: e+e- -> (gamma) e+e-
### ---------------------------------------------------------------------------
alg_bhabha = Algorithm.new("BhabhaLumi")
alg_bhabha.set_header(["BhabhaLumiAlg/BhabhaLumi.h"])

sel_bhabha = Selection.new
# --- Charged track selection ------------------------------------------------
# Exactly two oppositely charged tracks, |cos(theta)| < 0.8 measured by the MDC,
# originating from a cylinder of radius 1 cm about the IP perpendicular to the
# beam axis and +/-10 cm along it.  No particle identification is applied: the
# tracks are assigned electron / positron by their charge alone.
sel_bhabha.select_track do
            cos_theta 0.8      # |cos(theta)| < 0.8 (MDC polar angle)
            Vz        10.0     # |Vz| < 10 cm along the beam axis
            Vr        1.0      # Vr < 1 cm transverse to the beam axis
            nChrp     "==1"
            nChrn     "==1"
            nNet      "==0"
          end
          # --- EMC deposited-energy requirement (di-muon suppression) ---------
          # E_dep(e+-) > sqrt(s)/4.26 x 1.55 GeV.  The threshold is relative to
          # the c.m. energy, so the value below is the reference point at
          # sqrt(s) = 4.26 GeV; the per-energy rescaling is handled at the
          # analysis level (see the :efficiency_curve note).
          .for_each(:charged) do
            where { eraw < 1.55 }
            remove
          end
          # --- Momentum requirement (ISR resonance suppression) ---------------
          # p(e+-) > sqrt(s)/4.26 x 2 GeV/c, to remove ISR production of light
          # vector resonances (J/psi, psi(3686), ...) decaying to e+e-; again the
          # reference value at sqrt(s) = 4.26 GeV is used in BOSS.
          .for_each(:charged) do
            where { p < 2.0 }
            remove
          end
          # Nominal 4C fit of the two leptons to the c.m. four-momentum.  The
          # published luminosity does not use a kinematic fit (the signal yield
          # is obtained by counting selected events); the fit is kept as the DSL
          # end-point only.
          .kinematic_fit([:ep, :em]) do
            nominal
            constrain_four_momentum
            chi2_cut 200
          end

alg_bhabha.with_decay_card(decay_card_bhabha).apply(sel_bhabha)

### ---------------------------------------------------------------------------
### Algorithm 2 -- cross check: e+e- -> gamma gamma
### ---------------------------------------------------------------------------
alg_gammagamma = Algorithm.new("GGammaLumi")
alg_gammagamma.set_header(["GGammaLumiAlg/GGammaLumi.h"])

sel_gammagamma = Selection.new
# --- Charged track veto -------------------------------------------------------
# The di-gamma selection uses the electromagnetic calorimeter only (no MDC
# tracking information): no good charged track is allowed in the event.
sel_gammagamma.select_track do
                 Vz        10.0
                 Vr        1.0
                 nChrp     "==0"
                 nChrn     "==0"
                 nNet      "==0"
               end
               # --- Photon selection -----------------------------------------
               # At least two EMC clusters; the two most energetic clusters are
               # taken as the e+e- -> gamma gamma pair.  The deposited-energy
               # requirement of the EMC-only selection,
               # E > sqrt(s)/4.26 x 1.8 GeV (x 1.55 GeV above 4.420 GeV), is
               # applied with the sqrt(s) = 4.260 GeV reference value.
               .select_photon do
                 tdc_emc_start     0
                 tdc_emc_end       14
                 energyThreshold_b 0.025
                 energyThreshold_e 0.050
                 angle_to_track    10.0
                 nGam              ">=2"
               end
               .for_each(:gamma) do
                 where { energy < 1.8 }
                 remove
               end
               # Nominal 4C fit of the two photons to the c.m. four-momentum
               # (DSL end-point; the published cross check counts events).
               .kinematic_fit([:gamma, :gamma]) do
                 nominal
                 constrain_four_momentum
                 chi2_cut 200
               end

alg_gammagamma.with_decay_card(decay_card_gammagamma).apply(sel_gammagamma)

### ---------------------------------------------------------------------------
### Inexpressible BOSS-side procedures (captured for the systematics step)
### ---------------------------------------------------------------------------
alg_bhabha
  .note(:generator_model,
        "The signal and background QED samples are generated with BABAYAGA3.5 " \
        "(options: Ebeam = sqrt(s)/2, MinThetaAngle 20 deg, MaxThetaAngle 160 deg, " \
        "MinimumEnergy 0.04 GeV, MaximumAcollinearity 180 deg, RunningAlpha 1, FSR " \
        "switch 1), which models the running of the electromagnetic coupling " \
        "constant and final-state radiation. BABAYAGA3.5 is a dedicated QED " \
        "generator and is not expressible as an EvtGen decay card; the declared " \
        "card (psi(4260) -> e+ e- gamma with PHOTOS VLL) is the closest standard " \
        "BOSS approximation. The theoretical uncertainty of the generator cross " \
        "section is 0.5% and is taken as a systematic uncertainty.")
  .note(:efficiency_curve,
        "Both selection requirements are relative to the c.m. energy -- " \
        "E_dep(e+-) > sqrt(s)/4.26 x 1.55 GeV and p(e+-) > sqrt(s)/4.26 x 2 GeV/c -- " \
        "which is what makes the detection efficiency nearly independent of the " \
        "c.m. energy. The BOSS-level thresholds are the sqrt(s) = 4.26 GeV " \
        "reference values; the per-energy-point rescaling (and the change of the " \
        "EMC-only energy requirement from sqrt(s)/4.26 x 1.8 GeV to x 1.55 GeV " \
        "above 4.420 GeV) is applied at the analysis level.")
  .note(:cm_frame_variables,
        "The deposited energies, momenta and polar angles are determined in the " \
        "initial e+e- c.m. frame; because of the crossing angle of the beams the " \
        "c.m. system is slightly boosted relative to the laboratory frame, and " \
        "the boost is not applied to the BOSS-level track/shower variables.")
  .note(:background_veto,
        "A Bhabha event sample selected with EMC information only (no MDC " \
        "tracking) is used for the tracking-efficiency systematic uncertainty: at " \
        "least two EMC clusters, the two most energetic clusters assumed to be the " \
        "e+e- pair, deposited energies above sqrt(s)/4.26 x 1.8 GeV (x 1.55 GeV " \
        "above 4.420 GeV), |cos(theta_EMC)| < 0.8 and the azimuthal correlation " \
        "Delta phi = |phi_1 - phi_2| - 180 deg required to be in [-40,-5] deg or " \
        "[5,40] deg to remove the di-photon background. The efficiency with which " \
        "these EMC-only events pass the nominal track requirements is compared " \
        "between data and MC and the difference (0.39%) is the associated " \
        "systematic uncertainty. The cluster-angle and Delta phi conditions are " \
        "pair-level quantities that have no DSL construct.")
  .note(:no_kinematic_fit_in_paper,
        "The published analysis does NOT apply a kinematic fit: the number of " \
        "observed Bhabha events is obtained by counting the selected events, and " \
        "the luminosity follows from L = N_Bhabha^obs / (sigma_Bhabha x epsilon). " \
        "The nominal 4C fit declared here is only the DSL end-point and is not " \
        "part of the published event selection.")
  .note(:fit_model_variation,
        "Systematic uncertainties are evaluated by varying the selection: the " \
        "polar-angle requirement from |cos(theta)| < 0.8 to < 0.7 (0.38%), the " \
        "EMC energy requirement from sqrt(s)/4.26 x 1.55 GeV to x 1.71 GeV (0.09%) " \
        "and the momentum requirement from sqrt(s)/4.26 x 2 GeV/c to x 2.06 GeV/c " \
        "(0.43%). The largest deviation among the sub-samples is assigned to all " \
        "of them. Further sources: MC statistics 0.25%, beam energy 0.42% " \
        "(c.m. energy changed by 2 MeV, the difference between the tabulated " \
        "energy and the value measured with e+e- -> (gamma) mu+mu- events), " \
        "trigger efficiency 0.10% and generator 0.50%; total 0.97%.")
  .note(:trigger_efficiency,
        "The trigger efficiency for the Bhabha process is 100% with an " \
        "uncertainty of less than 0.1%.")

alg_gammagamma
  .note(:generator_model,
        "The di-gamma cross check uses the same BABAYAGA3.5 QED generator as the " \
        "nominal Bhabha measurement; the declared card " \
        "(psi(4260) -> gamma gamma PHSP) is the standard BOSS approximation.")
  .note(:background_veto,
        "The di-gamma cross check uses the same EMC-only selection as the " \
        "tracking-efficiency study, except that, to reduce the Bhabha background, " \
        "the azimuthal correlation between the two clusters must satisfy " \
        "Delta phi = |phi_1 - phi_2| - 180 deg in [-0.8 deg, 0.8 deg] (photons are " \
        "not deflected by the magnetic field). The two most energetic clusters are " \
        "taken as the photon pair and each must satisfy |cos(theta_EMC)| < 0.8. " \
        "These cluster-pair conditions have no DSL construct and are applied in " \
        "the ROOT analysis.")
  .note(:no_kinematic_fit_in_paper,
        "The cross-check luminosity L_ck is obtained by event counting, exactly as " \
        "the nominal Bhabha measurement; the declared 4C fit is only the DSL " \
        "end-point.")

### ---------------------------------------------------------------------------
### Execution
### ---------------------------------------------------------------------------
alg_bhabha.execute_on(data_points +
                      [incMC_4260] + exMC_bhabha)

alg_gammagamma.execute_on(data_points +
                          [incMC_4260] + exMC_gammagamma)
