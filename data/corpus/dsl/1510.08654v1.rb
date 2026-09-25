# =============================================================================
# BESIII arXiv:1510.08654v1
# Measurement of the centre-of-mass energies of the BESIII data samples
# collected around 4 GeV (2011-2014) with the di-muon process
# e+e- -> gamma_ISR/FSR mu+ mu-.
#
#   E_cms = M(mu+ mu-) + DeltaM_ISR/FSR
#
# The invariant mass of the muon pair is fitted run-by-run (a run is about one
# hour of data taking) with a Gaussian in the (-1 sigma, +2 sigma) window about
# the peak; the peak positions versus run number then give the time stability
# and, for samples 4260_1 and 4230_2, a linear-in-run-number drift.  The
# ISR/FSR mass shift is taken from dedicated BABAYAGA3.5 samples generated with
# the radiative corrections turned on and off.
#
# The muon momentum scale is validated with e+e- -> gamma_ISR J/psi,
# J/psi -> mu+ mu- (gamma_FSR) in the same data, the low-momentum charged-track
# measurement with D0 -> K- pi+ / D0bar -> K+ pi-, and the whole E_cms
# determination with e+e- -> pi+ pi- K+ K- and e+e- -> pi+ pi- p pbar.
#
# BOSS part only: dataset preparation and event selection up to and including
# the final kinematic fit.  The per-run Gaussian fits, the mass-shift fits and
# the extraction of E_cms are ROOT-level and are captured in the notes below.
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets
### ---------------------------------------------------------------------------
# 25 data sets (18 distinct c.m. energies) taken between 3.810 and 4.600 GeV.
# The dataset table lists 18 sample names but does not resolve the samples that
# were acquired separately at the same energy (the first/second/third sub-samples
# at 4.009, 4.230, 4.260 and 4.420 GeV); the pairs/triples are merged into a
# single sample name below and the sub-sample splitting is done by run number in
# the ROOT analysis, exactly as in Table I of the paper.
sample_names = %w[
  703_3810 703_3900 703_4009 703_4090 703_4190 703_4210 703_4220 703_4230
  703_4245 703_4260 703_4310 703_4360 703_4390 703_4420 703_4470 703_4530
  703_4575 703_4600
]

data_points = sample_names.map { |name| DatasetManager.real_data.find(name) }

# Inclusive MC is available only for a subset of the energy points (the
# "4190scan"/"4210scan"/"4220scan" entries are the 2012 sub-samples whose run
# ranges match this paper, the plain 4190/4210/4220 entries are the later
# round-10 data).
inc_mc_names = %w[
  703_4009 703_4190scan 703_4210scan 703_4220scan 703_4230
  703_4260 703_4360 703_4420 703_4600
]

inc_mc_points = inc_mc_names.map { |name| DatasetManager.inclusive_mc.find(name) }

### ---------------------------------------------------------------------------
### Decay cards
### ---------------------------------------------------------------------------
# Signal: e+e- -> gamma_ISR/FSR mu+ mu-.  Two reference samples per energy point
# are needed, one with the radiative corrections switched on and one with them
# switched off, because DeltaM_ISR/FSR is the difference of the fitted mu+mu-
# peak positions between the two.  BABAYAGA3.5 models the ISR (and the running
# of alpha) and PHOTOS models the FSR; the generator itself is not expressible
# as an EvtGen card, so the closest standard BOSS decomposition is declared and
# the generator settings are recorded in the :generator_model note.
# 'Particle vpho' is deliberately omitted: the same card is valid at every
# energy point and the per-point sqrt(s) is injected by execute_on.
decay_card_dimuon_rad = <<~DECAYCARD
    Decay psi(4260)
    1.0000 mu+ mu- gamma                            PHOTOS VLL;
    Enddecay

    End
DECAYCARD

decay_card_dimuon_plain = <<~DECAYCARD
    Decay psi(4260)
    1.0000 mu+ mu-                                  VLL;
    Enddecay

    End
DECAYCARD

# Momentum-scale validation: e+e- -> gamma_ISR J/psi, J/psi -> mu+ mu- with FSR.
decay_card_jpsi_rad = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma J/psi                              PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu-                                  PHOTOS VLL;
    Enddecay

    End
DECAYCARD

decay_card_jpsi_plain = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma J/psi                              PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu-                                  VLL;
    Enddecay

    End
DECAYCARD

# Low-momentum validation: D0 -> K- pi+ and its charge conjugate.
decay_card_d0 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D0 anti-D0                               PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+                                   PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi-                                   PHSP;
    Enddecay

    End
DECAYCARD

# E_cms cross-check channels.
decay_card_pipikk = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- K+ K-                            PHSP;
    Enddecay

    End
DECAYCARD

decay_card_ppbarpipi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- p+ anti-p-                       PHSP;
    Enddecay

    End
DECAYCARD

### ---------------------------------------------------------------------------
### Exclusive MC -- 50,000 events per energy point per sample (as in the paper)
### ---------------------------------------------------------------------------
exMC_dimuon_rad = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_dimuon_isrfsr_on"    # auto-suffixed per energy point
  config.events        = 50_000
  config.decay_card    = decay_card_dimuon_rad
  config.cross_section = :default
end

exMC_dimuon_plain = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_dimuon_isrfsr_off"
  config.events        = 50_000
  config.decay_card    = decay_card_dimuon_plain
  config.cross_section = :default
end

exMC_jpsi_rad = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ck_jpsi_isr_fsr_on"
  config.events        = 50_000
  config.decay_card    = decay_card_jpsi_rad
  config.cross_section = :default
end

exMC_jpsi_plain = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ck_jpsi_isr_fsr_off"
  config.events        = 50_000
  config.decay_card    = decay_card_jpsi_plain
  config.cross_section = :default
end

exMC_d0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ck_d0_kpi"
  config.events        = 50_000
  config.decay_card    = decay_card_d0
  config.cross_section = :default
end

exMC_pipikk = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ck_pipikk"
  config.events        = 50_000
  config.decay_card    = decay_card_pipikk
  config.cross_section = :default
end

exMC_ppbarpipi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ck_ppbarpipi"
  config.events        = 50_000
  config.decay_card    = decay_card_ppbarpipi
  config.cross_section = :default
end

### ---------------------------------------------------------------------------
### Algorithm 1 -- nominal E_cms measurement: e+e- -> gamma_ISR/FSR mu+ mu-
### ---------------------------------------------------------------------------
# Multi-energy scan: ECMS is NOT declared as an algorithm constant here -- every
# job takes its c.m. energy from the related dataset / jobOptions.
alg_dimuon = Algorithm.new("DimuonEcms")
alg_dimuon.set_header(["DimuonEcmsAlg/DimuonEcms.h"])

sel_dimuon = Selection.new
# --- Charged track selection -------------------------------------------------
# Exactly two good oppositely charged tracks, each consistent with originating
# from the run-dependent IP within 1 cm in the transverse plane and 10 cm along
# the beam, and both inside the barrel region |cos(theta)| < 0.80.
sel_dimuon.select_track do
           cos_theta 0.80    # |cos(theta)| < 0.80 (barrel only)
           Vz        10.0    # |Vz| < 10 cm along the beam
           Vr        1.0     # Vxy < 1 cm transverse to the beam
           nChrp     "==1"   # exactly one positive track
           nChrn     "==1"   # exactly one negative track
           nNet      "==0"   # opposite charges
         end
         # --- EMC energy deposition (radiative-Bhabha suppression) -------------
         # The energy deposited in the EMC by each charged track must be below
         # 0.4 GeV; a muon leaves only a minimum-ionising deposit, an electron of
         # the same momentum a full shower.
         .for_each(:charged) do
           where { eraw > 0.4 }
           remove
         end
         # --- Lepton identification -------------------------------------------
         # The two tracks are the muons of the di-muon final state (high-momentum
         # leptons: EMC energy / MUC depth path).
         .pid(method: :probability) do
           identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                           treat_as_electron_if_energy_above: 0.6
           nlp "==1"   # one l+
           nlm "==1"   # one l-
         end
         # Nominal 4C fit of the two leptons to the c.m. four-momentum.  The
         # published analysis does NOT apply a kinematic fit (the mu+mu- mass is
         # measured from the tracks themselves and then corrected for the
         # radiative mass shift); the fit is kept only as the DSL end-point.
         .kinematic_fit([:lp, :lm]) do
           nominal
           constrain_four_momentum
           chi2_cut 200
         end

alg_dimuon
  .note(:generator_model,
        "The di-muon signal and the DeltaM_ISR/FSR reference samples are generated " \
        "with BABAYAGA3.5 (50,000 events per sample, at every energy point) with the " \
        "radiative corrections switched on and off. BABAYAGA3.5 is a dedicated QED " \
        "generator modelling ISR and the running of the electromagnetic coupling " \
        "constant; it is not expressible as an EvtGen decay card, so the declared " \
        "cards (psi(4260) -> mu+ mu- gamma PHOTOS VLL for 'on' and " \
        "psi(4260) -> mu+ mu- VLL for 'off') are the closest standard BOSS " \
        "approximation. The same samples, processed with the selection below, give " \
        "the mass shift DeltaM_ISR/FSR used in E_cms = M(mu+mu-) + DeltaM_ISR/FSR.")
  .note(:dimuon_selection,
        "The di-muon selection requires exactly two good oppositely charged tracks, " \
        "each with |cos(theta)| < 0.80 and the point of closest approach to the " \
        "run-dependent IP within 1 cm in the transverse plane and 10 cm along the " \
        "beam. The EMC energy deposition of each track must be below 0.4 GeV " \
        "(radiative-Bhabha suppression) and the TOF timing difference between the " \
        "two tracks must satisfy |Delta t| < 4 ns. The invariant mass M(mu+mu-) is " \
        "fitted run-by-run with a Gaussian in the (-1 sigma, +2 sigma) window about " \
        "the peak; the peak positions versus run number give the time stability of " \
        "E_cms and, for the samples 4260_1 and 4230_2, a linear drift.")
  .note(:opening_angle_cut,
        "The opening angle between the two tracks is required to satisfy " \
        "178.60 deg < theta_mumu < 179.64 deg, which suppresses cosmic rays and " \
        "di-muon events with high-energy radiative photons. This is a two-track " \
        "pair-level quantity and has no DSL construct; it is applied at the ROOT " \
        "level.")
  .note(:tof_time_difference_veto,
        "Cosmic-ray background is further suppressed by requiring the TOF timing " \
        "difference between the two tracks to satisfy |Delta t| < 4 ns. The " \
        "difference of two tracks' flight times is a per-event quantity, not a " \
        "per-candidate DSL property.")
  .note(:emc_energy_deposition,
        "The requirement E_EMC < 0.4 GeV per charged track suppresses radiative " \
        "Bhabha events (e+e- -> (gamma) e+ e-), where the electron deposits its full " \
        "energy in the calorimeter, while a muon of the same momentum deposits only " \
        "a minimum-ionising amount. It is expressed above as a for_each filter on " \
        "the charged-track list using the eraw property.")
  .note(:background_level,
        "After the full selection the residual background is less than 0.001% of " \
        "the signal and is neglected in the analysis.")
  .note(:isr_fsr_mass_shift,
        "DeltaM_ISR/FSR is obtained as the difference in the fitted M(mu+mu-) peak " \
        "positions between the MC samples with ISR/FSR on and off; the same event " \
        "selection as for data is imposed on those samples. DeltaM_ISR/FSR versus " \
        "E_cms is fitted with a linear function, " \
        "DeltaM_ISR/FSR = (-3.53 +- 1.11) + (1.67 +- 0.28)e-3 x E_cms/MeV " \
        "(chi2/n.d.f. = 6.3/13); the fit is a ROOT-level procedure.")
  .note(:fsr_mass_shift,
        "The FSR-only shift DeltaM_FSR, obtained by comparing di-muon MC samples " \
        "with FSR on and off, is parameterised as " \
        "DeltaM_FSR = (-1.34 +- 0.84) + (0.56 +- 0.21)e-3 x E_cms, i.e. " \
        "0.79 +- 0.09 MeV at 3.81 GeV and 1.24 +- 0.14 MeV at 4.60 GeV.")
  .note(:run_dependent_fit,
        "The peak position of M(mu+mu-) is extracted run-by-run (one run is about " \
        "one hour of data taking) and then fitted as a function of the run number. " \
        "For the samples 4260_1 and 4230_2 the measured M(mu+mu-) drifts slowly and " \
        "is described by a linear function, " \
        "(4367.37 +- 53.53) + (-3.75 +- 1.80)e-3 x N_run MeV/c^2 and " \
        "(4316.81 +- 7.76) + (-2.87 +- 0.25)e-3 x N_run MeV/c^2 respectively; for " \
        "all other samples the average value is used. The samples 4009_1 and 4009_2 " \
        "(4420_2 and 4420_3) are separated because of a sudden drop of the average " \
        "energy between them. Run-by-run extraction and the time-stability fit are " \
        "ROOT-level.")
  .note(:systematics,
        "Systematic uncertainties on E_cms: momentum measurement 0.011% (from a " \
        "first-order fit of M^cor(J/psi) versus E_cms, whose largest deviation from " \
        "the nominal J/psi mass is 0.34 MeV/c^2), ISR/FSR correction 0.37 MeV/c^2 " \
        "(standard deviation of DeltaM_ISR/FSR versus E_cms), generator " \
        "0.036 +- 0.067 MeV/c^2 (BABAYAGA3.5 versus BABAYAGA@NLO), and time " \
        "stability below 0.25 MeV on average (largest difference between the " \
        "linear-in-run-number fit and the average value, per sample). The total " \
        "systematic uncertainty, obtained by adding the items in quadrature, is " \
        "below 0.8 MeV for every sample.")
  .note(:results,
        "The 18 luminosity-weighted average E_cms values (statistical and " \
        "systematic) are listed in Table II of the paper; e.g. " \
        "3810 -> 3807.65 +- 0.10 +- 0.58 MeV, 4009 -> 4007.62 +- 0.05 +- 0.66 MeV, " \
        "4230 -> 4226.26 +- 0.04 +- 0.65 MeV, 4260 -> 4257.97 +- 0.04 +- 0.66 MeV " \
        "and 4600 -> 4599.53 +- 0.07 +- 0.74 MeV. The overall precision of the " \
        "measurement is 0.8 MeV, and E_cms is found to be stable during data taking " \
        "for most samples.")
  .note(:no_kinematic_fit_in_paper,
        "The published analysis does NOT apply a kinematic fit: E_cms is obtained " \
        "from the measured M(mu+mu-) peak plus the ISR/FSR mass shift. The nominal " \
        "4C fit declared in the selection is only the DSL end-point and is not part " \
        "of the published event selection.")

alg_dimuon.with_decay_card(decay_card_dimuon_rad).apply(sel_dimuon)

### ---------------------------------------------------------------------------
### Algorithm 2 -- muon momentum validation: e+e- -> gamma_ISR J/psi
### ---------------------------------------------------------------------------
alg_jpsi = Algorithm.new("JpsiMomentumCheck")
alg_jpsi.set_header(["JpsiMomentumCheckAlg/JpsiMomentumCheck.h"])

sel_jpsi = Selection.new
# Same charged-track and EMC requirements as the di-muon measurement; the
# opening-angle requirement is replaced by the looser cosmic-ray cut
# cos(theta_mumu) > -0.98 (the J/psi events are not back-to-back).
sel_jpsi.select_track do
         cos_theta 0.80
         Vz        10.0
         Vr        1.0
         nChrp     "==1"
         nChrn     "==1"
         nNet      "==0"
       end
       .for_each(:charged) do
         where { eraw > 0.4 }
         remove
       end
       .pid(method: :probability) do
         identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                         treat_as_electron_if_energy_above: 0.6
         nlp "==1"
         nlm "==1"
       end
       # Nominal 4C fit of the two muons to the c.m. four-momentum (DSL
       # end-point only: the ISR photon carries away energy and the paper fits
       # the M(mu+mu-) peak without any kinematic constraint).
       .kinematic_fit([:lp, :lm]) do
         nominal
         constrain_four_momentum
         chi2_cut 200
       end

alg_jpsi
  .note(:jpsi_momentum_validation,
        "The high-momentum muon measurement is validated with J/psi -> mu+ mu- " \
        "candidates selected via e+e- -> gamma_ISR J/psi. The event selection is " \
        "that of the di-muon measurement (two good oppositely charged tracks, " \
        "|cos(theta)| < 0.80, Vxy < 1 cm, |Vz| < 10 cm, E_EMC < 0.4 GeV per track), " \
        "with the opening-angle requirement replaced by the cosmic-ray cut " \
        "cos(theta_mumu) > -0.98. The remaining background, from the radiative " \
        "di-muon process with the same final state, shows a smooth M(mu+mu-) " \
        "distribution and is described by a linear function; the J/psi peak is " \
        "fitted with a crystal-ball function. Adjacent samples with small " \
        "statistics are combined to reduce the fluctuation.")
  .note(:fsr_correction,
        "Because of final-state radiation (J/psi -> mu+ mu- gamma_FSR) the measured " \
        "M^obs(mu+mu-) lies slightly below the nominal J/psi mass. The shift " \
        "DeltaM_FSR is estimated from simulated e+e- -> gamma_ISR J/psi samples " \
        "(50,000 events per energy point, generated with PHOTOS with FSR on and off) " \
        "and is found to be independent of E_cms: a weighted average " \
        "DeltaM_FSR = (0.59 +- 0.04) MeV/c^2 is obtained by fitting the " \
        "FSR-on/FSR-off difference versus E_cms.")
  .note(:jpsi_mass_result,
        "The FSR-corrected J/psi mass M^cor(mu+mu-) is consistent for all data " \
        "samples, with an average of 3096.79 +- 0.08 MeV/c^2, in agreement with the " \
        "nominal J/psi mass within errors; the small difference is taken as the " \
        "momentum-measurement systematic uncertainty.")
  .note(:no_kinematic_fit_in_paper,
        "No kinematic fit is applied in the published analysis (this channel is a " \
        "pure momentum-scale validation); the declared 4C fit is only the DSL " \
        "end-point.")

alg_jpsi.with_decay_card(decay_card_jpsi_rad).apply(sel_jpsi)

### ---------------------------------------------------------------------------
### Algorithm 3 -- low-momentum validation: D0 -> K- pi+
### ---------------------------------------------------------------------------
alg_d0 = Algorithm.new("D0KpiCheck")
alg_d0.set_header(["D0KpiCheckAlg/D0KpiCheck.h"])

sel_d0 = Selection.new
sel_d0.select_track do
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     "==1"
         nChrn     "==1"
         nNet      "==0"
       end
       .pid(method: :probability) do
         prob_cut 0.001
         identify :kaon, against: [:pion, :proton]
         identify :pion, against: [:kaon, :proton]
         nkm  "==1"   # K- of the D0
         npip "==1"   # pi+ of the D0
       end
       .kinematic_fit([:km, :pip]) do
         nominal
         constrain_four_momentum
         chi2_cut 200
       end

alg_d0
  .note(:d0_mass_validation,
        "The low-momentum charged-track measurement is validated with the decay " \
        "channels D0 -> K- pi+ and D0bar -> K+ pi- (the charge-conjugate mode is " \
        "given its own algorithm below, since the kaon/pion charge assignment " \
        "differs). The measured mass, M^obs(K- pi+ / K+ pi-) = 1864.00 +- 0.7 MeV " \
        "(statistical uncertainty only), is consistent with the nominal D0/D0bar " \
        "mass with a deviation of 0.84 +- 0.71 MeV.")
  .note(:no_kinematic_fit_in_paper,
        "The published validation simply counts the K pi invariant mass from the " \
        "measured track momenta; the declared 4C fit is only the DSL end-point.")

alg_d0.with_decay_card(decay_card_d0).apply(sel_d0)

### ---------------------------------------------------------------------------
### Algorithm 4 -- low-momentum validation: D0bar -> K+ pi-
### ---------------------------------------------------------------------------
alg_d0bar = Algorithm.new("D0barKpiCheck")
alg_d0bar.set_header(["D0barKpiCheckAlg/D0barKpiCheck.h"])

sel_d0bar = Selection.new
sel_d0bar.select_track do
            cos_theta 0.93
            Vz        10.0
            Vr        1.0
            nChrp     "==1"
            nChrn     "==1"
            nNet      "==0"
          end
          .pid(method: :probability) do
            prob_cut 0.001
            identify :kaon, against: [:pion, :proton]
            identify :pion, against: [:kaon, :proton]
            nkp  "==1"   # K+ of the D0bar
            npim "==1"   # pi- of the D0bar
          end
          .kinematic_fit([:kp, :pim]) do
            nominal
            constrain_four_momentum
            chi2_cut 200
          end

alg_d0bar
  .note(:d0_mass_validation,
        "Charge-conjugate counterpart of the D0 -> K- pi+ validation; the two " \
        "channels share the measured mass M^obs(K- pi+ / K+ pi-) = " \
        "1864.00 +- 0.7 MeV used in the paper.")
  .note(:no_kinematic_fit_in_paper,
        "No kinematic fit is applied in the published validation; the declared 4C " \
        "fit is only the DSL end-point.")

alg_d0bar.with_decay_card(decay_card_d0).apply(sel_d0bar)

### ---------------------------------------------------------------------------
### Algorithm 5 -- cross check: e+e- -> pi+ pi- K+ K-
### ---------------------------------------------------------------------------
alg_pipikk = Algorithm.new("EcmsCheckPiPiKK")
alg_pipikk.set_header(["EcmsCheckPiPiKKAlg/EcmsCheckPiPiKK.h"])

sel_pipikk = Selection.new
sel_pipikk.select_track do
            cos_theta 0.93
            Vz        10.0
            Vr        1.0
            nChrp     "==2"
            nChrn     "==2"
            nNet      "==0"
          end
          .pid(method: :probability) do
            prob_cut 0.001
            identify :kaon, against: [:pion, :proton]
            identify :pion, against: [:kaon, :proton]
            nkp  "==1"
            nkm  "==1"
            npip "==1"
            npim "==1"
          end
          .kinematic_fit([:pip, :pim, :kp, :km]) do
            nominal
            constrain_four_momentum
            chi2_cut 200
          end

alg_pipikk
  .note(:ecms_cross_check,
        "The E_cms values obtained from the di-muon process are cross checked with " \
        "the processes e+e- -> pi+ pi- K+ K- and e+e- -> pi+ pi- p pbar, in which " \
        "E_cms is estimated from the corrected invariant masses " \
        "M^cor(pi+pi-K+K-) and M^cor(pi+pi-p pbar) in the same way as for the " \
        "di-muon channel. Both are found to be consistent with the di-muon result, " \
        "the largest deviation being 0.53 +- 0.75 MeV in the 4420 sample.")
  .note(:no_kinematic_fit_in_paper,
        "The published cross check uses the measured four-momenta of the four " \
        "charged tracks directly (with the same ISR/FSR-type correction as for the " \
        "di-muon channel); the declared 4C fit is only the DSL end-point.")

alg_pipikk.with_decay_card(decay_card_pipikk).apply(sel_pipikk)

### ---------------------------------------------------------------------------
### Algorithm 6 -- cross check: e+e- -> pi+ pi- p pbar
### ---------------------------------------------------------------------------
alg_ppbarpipi = Algorithm.new("EcmsCheckPPbarPiPi")
alg_ppbarpipi.set_header(["EcmsCheckPPbarPiPiAlg/EcmsCheckPPbarPiPi.h"])

sel_ppbarpipi = Selection.new
sel_ppbarpipi.select_track do
               cos_theta 0.93
               Vz        10.0
               Vr        1.0
               nChrp     "==2"   # pi+ and p
               nChrn     "==2"   # pi- and pbar
               nNet      "==0"
             end
             .pid(method: :probability) do
               prob_cut 0.001
               identify :proton, against: [:kaon, :pion]
               identify :pion, against: [:kaon, :proton]
               nprp "==1"
               nprm "==1"
               npip "==1"
               npim "==1"
             end
             .kinematic_fit([:pip, :pim, :prp, :prm]) do
               nominal
               constrain_four_momentum
               chi2_cut 200
             end

alg_ppbarpipi
  .note(:ecms_cross_check,
        "Second E_cms cross-check channel, e+e- -> pi+ pi- p pbar, analysed through " \
        "the corrected invariant mass M^cor(pi+pi-p pbar) exactly as the di-muon " \
        "channel; the result is consistent with the di-muon E_cms, with the largest " \
        "deviation (0.53 +- 0.75 MeV) observed in the 4420 sample.")
  .note(:no_kinematic_fit_in_paper,
        "No kinematic fit is applied in the published cross check; the declared 4C " \
        "fit is only the DSL end-point.")

alg_ppbarpipi.with_decay_card(decay_card_ppbarpipi).apply(sel_ppbarpipi)

### ---------------------------------------------------------------------------
### Execution
### ---------------------------------------------------------------------------
common_datasets = data_points + inc_mc_points

alg_dimuon.execute_on(common_datasets + exMC_dimuon_rad + exMC_dimuon_plain)
alg_jpsi.execute_on(common_datasets + exMC_jpsi_rad + exMC_jpsi_plain)
alg_d0.execute_on(common_datasets + exMC_d0)
alg_d0bar.execute_on(common_datasets + exMC_d0)
alg_pipikk.execute_on(common_datasets + exMC_pipikk)
alg_ppbarpipi.execute_on(common_datasets + exMC_ppbarpipi)
