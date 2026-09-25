# =============================================================================
# BESIII measurement of the Born cross section of e+e- -> p pbar
# (arXiv:1504.02680v3)  --  BOSS part (dataset preparation + event selection)
#
# The process is a continuum production: the Born cross section is measured at
# 12 c.m. energies between 2232.4 and 3671.0 MeV, therefore the signal is
# modelled with the ConExc generator (ISR up to second order + the measured
# sigma_0(m)); mode 0 of the ConExc mode table is p pbar, valid 1.877-4.500 GeV.
# ConExc also provides the ISR / vacuum-polarisation factors f_ISR, f_vacuum
# used in sigma_Born = N_sig / (L * eps * f_ISR * f_vacuum).
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets (multi-energy scan -> cross-section measurement)
### ---------------------------------------------------------------------------
energy_points = [
  DatasetManager.real_data.find("713_Rscan_2232"),  # sqrt(s) = 2232.4 MeV
  DatasetManager.real_data.find("713_Rscan_2396"),  # sqrt(s) = 2400.0 MeV
  DatasetManager.real_data.find("713_Rscan_2800"),  # sqrt(s) = 2800.0 MeV
  DatasetManager.real_data.find("713_Rscan_3020"),  # sqrt(s) = 3050.0 / 3060.0 MeV
  DatasetManager.real_data.find("713_Rscan_3080"),  # sqrt(s) = 3080.0 MeV
  DatasetManager.real_data.find("709_3650"),        # sqrt(s) = 3500.0-3671.0 MeV (5 points)
]

incMC_points = [
  DatasetManager.inclusive_mc.find("713_Rscan_2232"),
  DatasetManager.inclusive_mc.find("713_Rscan_2396"),
  DatasetManager.inclusive_mc.find("713_Rscan_2800"),
  DatasetManager.inclusive_mc.find("713_Rscan_3020"),
  DatasetManager.inclusive_mc.find("713_Rscan_3080"),
  DatasetManager.inclusive_mc.find("709_3650"),
]

### ---------------------------------------------------------------------------
### Decay cards
### ---------------------------------------------------------------------------
# Signal: e+e- -> p pbar at the scan point. 'Particle vpho' is deliberately
# OMITTED so that execute_on injects the correct per-point sqrt(s) from each
# related dataset -- the same card is then valid at every energy point.
decay_card_signal = <<~DECAYCARD
  Decay vpho
  1 ConExc 0;
  Enddecay
  Decay vhdr
  1 p+ anti-p- PHSP;
  Enddecay
  End
DECAYCARD

# Hadronic background: e+e- -> pi+ pi- (uniform phase space, as in the paper)
decay_card_bkg_pipi = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Hadronic background: e+e- -> K+ K-
decay_card_bkg_kk = <<~DECAYCARD
  Decay psi(4260)
  1.0000 K+ K- PHSP;
  Enddecay
  End
DECAYCARD

# Hadronic background: e+e- -> p pbar pi0
decay_card_bkg_pppi0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 p+ anti-p- pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Hadronic background: e+e- -> p pbar pi0 pi0
decay_card_bkg_pppi0pi0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 p+ anti-p- pi0 pi0 PHSP;
  Enddecay
  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Hadronic background: e+e- -> Lambda anti-Lambda
decay_card_bkg_llbar = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda0 anti-Lambda0 PHSP;
  Enddecay
  Decay Lambda0
  1.000 p+ pi- PHSP;
  Enddecay
  Decay anti-Lambda0
  1.000 anti-p- pi+ PHSP;
  Enddecay
  End
DECAYCARD

### ---------------------------------------------------------------------------
### Exclusive MC (signal: >10x the data statistics at every energy point)
### ---------------------------------------------------------------------------
exMC_signal = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "sig_ppbar_conexc"   # auto-suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default             # required field; inert for ConExc
end

# Background MC: equivalent luminosities at least as large as the data.
exMC_bkg_pipi = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "bkg_pipi_phsp"
  config.events        = 50_000
  config.decay_card    = decay_card_bkg_pipi
  config.cross_section = :default
end

exMC_bkg_kk = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "bkg_kk_phsp"
  config.events        = 50_000
  config.decay_card    = decay_card_bkg_kk
  config.cross_section = :default
end

exMC_bkg_pppi0 = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "bkg_pppi0_phsp"
  config.events        = 50_000
  config.decay_card    = decay_card_bkg_pppi0
  config.cross_section = :default
end

exMC_bkg_pppi0pi0 = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "bkg_pppi0pi0_phsp"
  config.events        = 50_000
  config.decay_card    = decay_card_bkg_pppi0pi0
  config.cross_section = :default
end

exMC_bkg_llbar = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "bkg_llbar_phsp"
  config.events        = 50_000
  config.decay_card    = decay_card_bkg_llbar
  config.cross_section = :default
end

### ---------------------------------------------------------------------------
### Algorithm + event selection
### ---------------------------------------------------------------------------
# NOTE: this is a multi-energy cross-section measurement, so ECMS is NOT set as
# a constant here -- each job takes its c.m. energy from the dataset/jobOptions.
alg_ppbar = Algorithm.new("PpbarXS")
alg_ppbar.set_header(["PpbarXSAlg/PpbarXS.h"])

event_selection = Selection.new

# --- Charged track selection -------------------------------------------------
# Good track: within MDC coverage |cos(theta)| < 0.93, point of closest
# approach to the IP within 1 cm in the plane perpendicular to the beam and
# within +/-10 cm along the beam direction. Exactly two good charged tracks,
# one proton and one anti-proton.
event_selection.select_track do
  cos_theta 0.93     # |cos(theta)| < 0.93
  Vz        10.0     # |Vz| < 10 cm along the beam
  Vr        1.0      # |Vr| < 1 cm in the transverse plane
  nChrp     "==1"    # one positive track  (proton candidate)
  nChrn     "==1"    # one negative track  (anti-proton candidate)
  nNet      "==0"    # net charge zero
end

# --- Particle identification -------------------------------------------------
# Combined dE/dx (MDC) and TOF information gives PID probabilities for the
# pion, kaon and proton hypotheses; the type with the highest probability is
# assigned. One proton and one anti-proton are required.
event_selection.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprp "==1"
  nprm "==1"
end

# --- E/p requirement (Bhabha suppression) ------------------------------------
# E/p of each proton candidate must be smaller than 0.5, where E is the energy
# deposited in the EMC and p the momentum measured in the MDC.
event_selection.for_each(:prp) do
  define(:e_over_p) { eraw / p }
  where { e_over_p >= 0.5 }
  remove
end

event_selection.for_each(:prm) do
  define(:e_over_p) { eraw / p }
  where { e_over_p >= 0.5 }
  remove
end

# --- Final state reconstruction / endpoint -----------------------------------
# The published analysis extracts the signal yield by counting events inside a
# momentum window and does not perform a kinematic fit; the nominal 4C fit
# below is kept only as the DSL end-point (loose chi2, tight cut in ROOT).
event_selection.kinematic_fit([:prp, :prm]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

### ---------------------------------------------------------------------------
### Inexpressible BOSS-side procedures (captured for the systematics step)
### ---------------------------------------------------------------------------
alg_ppbar
  .note(:tof_time_difference_veto,
        "cosmic-ray background rejected by requiring |T_trk1 - T_trk2| < 4 ns, " \
        "the difference of the TOF-measured flight times of the two tracks; a " \
        "per-event quantity, not expressible as a per-candidate DSL property")
  .note(:bhabha_suppression_cut,
        "for data samples with sqrt(s) > 2400.0 MeV the proton is additionally " \
        "required to satisfy cos(theta_p) < 0.8 to suppress Bhabha background; " \
        "the cut is c.m.-energy dependent and is applied per energy point")
  .note(:opening_angle_cut,
        "c.m.-energy dependent requirement on the opening angle between proton and " \
        "anti-proton: theta_ppbar > 178 deg for sqrt(s) <= 2400.0 MeV and " \
        "theta_ppbar > 179 deg for sqrt(s) > 2400.0 MeV")
  .note(:momentum_window,
        "|p_meas - p_exp| < 5 sigma_p for the proton and the anti-proton, where " \
        "p_exp is the expected momentum in the c.m. system at the given sqrt(s) and " \
        "sigma_p the corresponding momentum resolution; used to count the signal yield")
  .note(:no_kinematic_fit_in_paper,
        "the published analysis does NOT apply a kinematic fit: signal yields are " \
        "obtained by event counting after the momentum window. The nominal 4C fit " \
        "in this selection is only the DSL end-point and is not part of the " \
        "published event selection")
  .note(:beam_associated_background,
        "beam-associated background (beam-gas, beam-pipe, Touschek) studied with " \
        "separated-beam data at sqrt(s) = 2400.0 and 3400.0 MeV, normalised by " \
        "data-taking time and efficiencies; found to be negligible at all energies")
  .note(:qed_background_generator,
        "QED backgrounds e+e- -> l+l- (l = e, mu) and e+e- -> gamma gamma are " \
        "generated with the Babayaga generator, which is not expressible as an " \
        "EvtGen exclusive decay card")
  .note(:efficiency_curve,
        "detection efficiency depends on |G_E/G_M| through the cos(theta_p) " \
        "distribution: for the energy points where |G_E/G_M| is measured the " \
        "measured ratio is used as MC input, otherwise the efficiency is the " \
        "average of the values obtained with |G_E| = 0 and |G_M| = 0")
  .note(:non_resonant_qed_subtraction,
        "the non-resonant QED contribution is estimated from data taken far from " \
        "any resonance and from psi(3770) data, and a simultaneous fit gives " \
        "N_phipi0(3.097) < 5.8 at 90% C.L.; the contribution is neglected")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

### ---------------------------------------------------------------------------
### Execution
### ---------------------------------------------------------------------------
alg_ppbar.execute_on(energy_points + incMC_points +
                     exMC_signal +
                     exMC_bkg_pipi + exMC_bkg_kk + exMC_bkg_pppi0 +
                     exMC_bkg_pppi0pi0 + exMC_bkg_llbar)
