# BESIII paper arXiv:1205.4284v1
# Study of psi' -> gamma chi_{c0,2}, chi_{c0,2} -> gamma gamma with 106 M psi' events,
# and measurement of the two-photon widths Gamma_{gamma gamma}(chi_{c0,2}) and
# f_{0/2} = Gamma^{lambda=0}/Gamma^{lambda=2} of chi_{c2} -> gamma gamma.
# Final state: three photons, no charged tracks. Nominal fit: 4C kinematic fit.

### Dataset description ###
psip_data   = DatasetManager.real_data.find("709_3686")     # psi(2S) data, 106 M psi' events (156.4 pb^-1)
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")  # psi(2S) inclusive MC
off_data    = DatasetManager.real_data.find("709_3650")     # off-psi' data at sqrt(s) = 3.65 GeV (44.1 pb^-1)
off_incMC   = DatasetManager.inclusive_mc.find("709_3650")  # corresponding continuum inclusive MC
psi3770_data  = DatasetManager.real_data.find("712_3773")   # psi(3770) data at sqrt(s) = 3.773 GeV (921.8 pb^-1)
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

### Decay cards (EvtGen format) ###
# Signal: psi(2S) -> gamma chi_c0, chi_c0 -> gamma gamma.
# The E1 radiative transition psi(2S) -> gamma chi_c0 is generated with a
# (1 + cos^2 theta) angular distribution; chi_c0 -> gamma gamma is generated
# with a uniform (PHSP) angular distribution.
decay_card_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  chi_c0            PHSP;
    Enddecay

    Decay chi_c0
    1.0000  gamma  gamma             PHSP;
    Enddecay

    End
DECAYCARD

# Signal: psi(2S) -> gamma chi_c2, chi_c2 -> gamma gamma.
# The full angular amplitudes for psi' -> gamma chi_c2 and chi_c2 -> gamma gamma
# are generated with the measured parameters x = A1/A0 = 1.55, y = A2/A0 = 2.10
# and f_{0/2} = 0 (pure helicity-two two-photon state), i.e. the angular
# distribution is not pure PHSP. The generator-level parameterisation is captured
# in the note below; the card keeps the decay topology.
decay_card_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  chi_c2            PHSP;
    Enddecay

    Decay chi_c2
    1.0000  gamma  gamma             PHSP;
    Enddecay

    End
DECAYCARD

# Control sample used to extract the signal line shapes in the E_gamma1 spectrum:
# psi(2S) -> gamma chi_{c0,2}, chi_{c0,2} -> K+ K- (purity > 99.2%).
decay_card_chic0_kk = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  chi_c0            PHSP;
    Enddecay

    Decay chi_c0
    1.0000  K+  K-                   PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic2_kk = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  chi_c2            PHSP;
    Enddecay

    Decay chi_c2
    1.0000  K+  K-                   PHSP;
    Enddecay

    End
DECAYCARD

# Peaking-background channels entering the chi_{c0,2} signal windows:
# chi_{c0,c2} -> pi0 pi0 and chi_{c0,2} -> eta eta, with pi0/eta -> gamma gamma
# where two of the photons are soft and escape detection.
decay_card_bkg_pi0pi0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  chi_c0            PHSP;
    Enddecay

    Decay chi_c0
    1.0000  pi0  pi0                 PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma             PHSP;
    Enddecay

    End
DECAYCARD

decay_card_bkg_etaeta = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma  chi_c2            PHSP;
    Enddecay

    Decay chi_c2
    1.0000  eta  eta                 PHSP;
    Enddecay

    Decay eta
    1.0000  gamma  gamma             PHSP;
    Enddecay

    End
DECAYCARD

# Continuum background e+e- -> gamma gamma (gamma), described by the Babayaga QED
# generator; simulated at the off-psi' energy point (3.65 GeV) and at 3.773 GeV.
decay_card_continuum = <<~DECAYCARD
    Decay psi(4260)
    1.0000  gamma  gamma             PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples ###
exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_chic0_chic0_to_gammagamma"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_chic0
  config.cross_section   = :default
end
exMC_chic0.save_to_config(format: :yaml, file_path: 'exMC_chic0_gammagamma')

exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_chic2_chic2_to_gammagamma"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_chic2
  config.cross_section   = :default
end
exMC_chic2.save_to_config(format: :yaml, file_path: 'exMC_chic2_gammagamma')

exMC_chic0_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_chic0_chic0_to_KK"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_chic0_kk
  config.cross_section   = :default
end
exMC_chic0_kk.save_to_config(format: :yaml, file_path: 'exMC_chic0_KK')

exMC_chic2_kk = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_chic2_chic2_to_KK"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_chic2_kk
  config.cross_section   = :default
end
exMC_chic2_kk.save_to_config(format: :yaml, file_path: 'exMC_chic2_KK')

exMC_bkg_pi0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_chic_pi0pi0_peaking_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_pi0pi0
  config.cross_section   = :default
end
exMC_bkg_pi0pi0.save_to_config(format: :yaml, file_path: 'exMC_peaking_pi0pi0')

exMC_bkg_etaeta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_gamma_chic_etaeta_peaking_bkg"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_etaeta
  config.cross_section   = :default
end
exMC_bkg_etaeta.save_to_config(format: :yaml, file_path: 'exMC_peaking_etaeta')

# Continuum e+e- -> gamma gamma (gamma) MC at the two off-resonance energy points.
exMC_continuum = DatasetManager.create_exclusive_mc_for([off_data, psi3770_data]) do |config|
  config.sample_name   = "continuum_gammagamma_babayaga"
  config.events        = 100_000
  config.decay_card    = decay_card_continuum
  config.cross_section = :default
end

### Event selection (BOSS) ###
# Single algorithm: chi_c0 and chi_c2 share the identical final state
# (gamma gamma gamma, no charged tracks) and identical selection chain.
alg_name = "PsipToGammaChiGammaGamma"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.686]})
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  nTot     "==0"     # no detected charged particles in the event
                }
               .select_photon {
                  tdc_emc_start     0     # EMC timing window
                  tdc_emc_end      14
                  angle_to_track   10.0   # min angle to the nearest charged track
                  energyThreshold_b 0.025 # E > 25 MeV for photons in the barrel
                  energyThreshold_e 0.050 # E > 25 MeV extrapolated to the end-cap region
                  nGam             "==3"  # exactly three photon candidates
                }
               .kinematic_fit([:gamma, :gamma, :gamma]) {
                  nominal                  # nominal fit: corrected four-momenta used downstream
                  constrain_four_momentum  # 4C energy-momentum conservation constraint
                  chi2_cut 200             # loose chi^2 cut in BOSS; chi^2_4C <= 80 applied in ROOT
               }

# Procedures that cannot be expressed in the DSL are recorded as notes.
alg.note(:leading_photon_assignment,
         "Of the three photon candidates, the smallest-energy photon is assigned to the " \
         "radiated photon gamma_1 from psi' -> gamma_1 chi_{c0,2}, while the second-largest and " \
         "the largest energy photons are assigned to gamma_2 gamma_3 from chi_{c0,2} -> gamma gamma. " \
         "The E_gamma1 (radiated photon energy) spectrum is the fit variable; the E_gamma1 " \
         "resolution is sigma = 6.74 +/- 0.29 MeV for chi_c0 and 3.91 +/- 0.09 MeV for chi_c2.")
   .note(:photon_fiducial_region,
         "Photon candidates are required to have |cos theta| < 0.75 with respect to the e+ beam " \
         "direction, in addition to the EMC fiducial/shower-quality criteria. This requirement " \
         "suppresses the continuum e+e- -> gamma gamma (gamma) background, whose two energetic " \
         "photons are mostly emitted in the forward and backward regions. The angular acceptance " \
         "cut is not expressible via the select_photon energy/angle-to-track keywords.")
   .note(:event_vertex,
         "The average event vertex of each run is used as the origin for the selected candidates.")
   .note(:signal_mc_angular_distributions,
         "Signal MC angular distributions: psi(2S) -> gamma_1 chi_c0 is generated with a " \
         "(1 + cos^2 theta) distribution (pure E1 transition); chi_c0 -> gamma gamma is generated " \
         "with a uniform angular distribution; psi(2S) -> gamma_1 chi_c2, chi_c2 -> gamma gamma is " \
         "generated with the full angular amplitudes using x = A1/A0 = 1.55, y = A2/A0 = 2.10 and " \
         "f_{0/2} = 0 (pure helicity-two two-photon state). These generator-level shapes are not " \
         "expressible in the EvtGen decay card used here.")
   .note(:signal_efficiency,
         "Efficiencies from MC: epsilon(chi_c0) = (35.4 +/- 0.06)% and epsilon(chi_c2) = (38.0 +/- 0.07)%; " \
         "the difference is due primarily to the different angular distributions.")
   .note(:background_lineshape,
         "The shape of the dominant non-peaking continuum background in the E_gamma1 spectrum is " \
         "described by the data-driven function f_bg(E_gamma1) = p0 + p1*E_gamma1 + p2*(E_gamma1)^a, " \
         "whose parameters are obtained from a fit to the psi(3770) data at sqrt(s) = 3.773 GeV " \
         "(921.8 pb^-1) and cross-checked against the 44.1 pb^-1 off-psi' sample at 3.65 GeV and " \
         "against Babayaga e+e- -> gamma gamma (gamma) MC.")
   .note(:signal_shape_from_control_sample,
         "The chi_c0 and chi_c2 signal line shapes in the E_gamma1 spectrum are taken from the " \
         "smoothed histograms of the nearly background-free psi' -> gamma_1 chi_{c0,2}, " \
         "chi_{c0,2} -> K+ K- control sample (purity > 99.2%).")
   .note(:peaking_background_subtraction,
         "Peaking backgrounds from chi_{c0,c2} -> pi0 pi0 and eta eta (with two soft photons " \
         "undetected or outside the fiducial volume) contribute 25.8 and 7.8 events to the chi_c0 " \
         "and chi_c2 signals respectively; they are subtracted from the fitted yields, giving " \
         "N(chi_c0) = 813 +/- 63 and N(chi_c2) = 1131 +/- 66. The yields are obtained from an " \
         "unbinned maximum-likelihood fit to the E_gamma1 spectrum.")
   .note(:helicity_component_analysis,
         "The helicity amplitude analysis of psi' -> gamma chi_c2, chi_c2 -> gamma gamma is performed " \
         "with an unbinned maximum-likelihood fit to the three-dimensional angular distribution " \
         "(cos theta_1, cos theta_2, phi_2) over 0.09 < E_gamma1 < 0.15 GeV, using the twelve angular " \
         "factors a_1..a_12 averaged over phase-space MC to account for detector acceptance, and with " \
         "the x, y parameters fixed to the BESIII values x = 1.55 +/- 0.05 +/- 0.07 and " \
         "y = 2.10 +/- 0.07 +/- 0.05. The background is subtracted via ln L_s = ln L - ln L_b, with " \
         "sidebands (0.07, 0.08) GeV and (0.16, 0.20) GeV, yielding f_{0/2} = 0.00 +/- 0.02 (stat.) " \
         "+/- 0.02 (syst.). This angular fit is not expressible in the DSL.")
   .note(:systematic_uncertainties,
         "Relative systematic uncertainties on the branching fractions (chi_c0 / chi_c2): number of " \
         "psi' events 4.0% (common), neutral trigger efficiency 0.1% (common), photon detection 1.5%, " \
         "kinematic fit 1.0% (common), resonance fitting 3.2% / 2.9%, peaking background 0.3% / 0.1%, " \
         "helicity-two assumption 0.4% (chi_c2 only).")

alg.with_decay_card(decay_card_chic0).apply(event_selection)

root_files = alg.execute_on([psip_data, psip_incMC, off_data, off_incMC,
                             psi3770_data, psi3770_incMC,
                             exMC_chic0, exMC_chic2, exMC_chic0_kk, exMC_chic2_kk,
                             exMC_bkg_pi0pi0, exMC_bkg_etaeta] + exMC_continuum)
