# BESIII paper 1607.00738v1
# Determination of the number of J/psi events with inclusive J/psi decays
# (2012 sample: (1086.9 +- 6.0) x 10^6; 2009 sample: (223.7 +- 1.4) x 10^6).
#
# This is an inclusive counting analysis: N_J/psi = (N_sel - N_bg) /
# (eps_trig * eps_data^{psi(3686)} * f_cor).  There is no exclusive final state
# and no kinematic fit; the BOSS-side scope is the inclusive J/psi event
# selection (charged tracks, photons, event-level QED/beam-gas suppression)
# applied both to the J/psi data and to the psi(3686) -> pi+ pi- J/psi sample
# used to measure the detection efficiency experimentally.

### Dataset preparation ###
# J/psi data and inclusive MC (2009 and 2012 samples) at sqrt(s) = 3.097 GeV
jpsi_data   = DatasetManager.real_data.find("708_3097")
jpsi_incMC  = DatasetManager.inclusive_mc.find("708_3097")

# psi(3686) data and inclusive MC at sqrt(s) = 3.686 GeV; the dedicated
# psi(3686) sample taken on May 26, 2012 is used to determine eps_data^{psi(3686)}
# via psi(3686) -> pi+ pi- J/psi, J/psi -> inclusive
psip_data   = DatasetManager.real_data.find("709_3686")
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")

# Exclusive MC: psi(3686) -> pi+ pi- J/psi with the J/psi decaying inclusively.
# Provides eps_MC^{psi(3686)} (denominator of the correction factor f_cor).
decay_card_psip_to_pipi_jpsi = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- J/psi                    VVPIPI;
    Enddecay

    Decay J/psi
    1.0000 anything                         PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: inclusive J/psi produced directly in the e+ e- collision
# (J/psi at rest), providing eps_MC^{J/psi} (numerator of f_cor).
decay_card_jpsi_inclusive = <<~DECAYCARD
    Decay J/psi
    1.0000 anything                         PHSP;
    Enddecay

    End
DECAYCARD

exMC_psip_to_pipi_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_pipi_jpsi_inclusive"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_psip_to_pipi_jpsi
  config.cross_section   = :default
end

exMC_jpsi_inclusive = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_inclusive_direct"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_jpsi_inclusive
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name       = "JpsiInclusiveEvents"
jpsi_algorithm = Algorithm.new(alg_name)
jpsi_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
              .set_constant({"ECMS" => [:double, 3.097]})
              .set_alias({"std::vector<double>" => "Vdouble"})

# Inclusive J/psi event selection (Sec. 2 of the paper).
#   charged tracks: |cos theta| < 0.93, Vr < 1 cm, |Vz| < 15 cm, p < 2.0 GeV/c,
#                   at least two charged tracks in the event;
#   photons: E_dep > 25 MeV (barrel, |cos theta| < 0.83) / > 50 MeV
#            (endcap, 0.86 < |cos theta| < 0.93), EMC cluster timing 0 < T <= 700 ns.
event_selection = Selection.new
event_selection.select_track {
                  cos_theta    0.93     # |cos(theta)| < 0.93 for MDC tracks
                  Vz          15.0      # |Vz| < 15 cm (beam direction)
                  Vr           1.0      # Vr < 1 cm (radial direction to run-dependent IP)
                  nTot        ">=2"     # at least two charged tracks in the event
                }
               .select_photon {
                  energyThreshold_b 0.025  # barrel EMC threshold 25 MeV (|cos(theta)| < 0.83)
                  energyThreshold_e 0.050  # endcap EMC threshold 50 MeV (0.86 < |cos(theta)| < 0.93)
                  tdc_emc_start     0      # EMC cluster timing T > 0
                  tdc_emc_end      14      # EMC cluster timing T <= 700 ns (14 x 50 ns)
               }

# Remaining event-level criteria of the inclusive J/psi selection, the
# efficiency/background bookkeeping and the psi(3686) tag side.  These have no
# dedicated DSL construct in the BOSS selection chain and are recorded as notes.
jpsi_algorithm
  .note(:track_momentum_cut,
        "Each good charged track is required to have momentum p < 2.0 GeV/c in the MDC. " \
        "This upper limit is not expressible through the select_track cos_theta/Vz/Vr keyword set.")
  .note(:visible_energy_cut,
        "Event-level cut E_vis > 1.0 GeV, where E_vis is the sum of charged-particle energies " \
        "(computed from the track momenta under the pion-mass hypothesis) and of the neutral " \
        "shower energies deposited in the EMC.  Removes about one third of the background while " \
        "retaining 99.4% of the signal.")
  .note(:two_prong_lepton_veto,
        "For events with exactly two charged tracks, both track momenta are required to be " \
        "< 1.5 GeV/c; this removes the Bhabha (e+e- -> e+e-) and dimuon (e+e- -> mu+mu-) " \
        "lepton-pair cluster at p ~ 1.55 GeV/c in the two-prong momentum scatter plot.")
  .note(:emc_deposit_per_track_veto,
        "For each good charged track the energy deposited in the EMC is required to be < 1 GeV, " \
        "further rejecting Bhabha events (the characteristic ~1.5 GeV peak in the per-track " \
        "EMC deposition spectrum).")
  .note(:selection_yield,
        "After the full inclusive selection: N_sel = (854.60 +- 0.03) x 10^6 candidates from the " \
        "2012 J/psi data and N_sel = (179.63 +- 0.01) x 10^6 from the 2009 J/psi data.")
  .note(:background_normalisation,
        "The background from QED processes, cosmic rays and beam-gas interactions is estimated " \
        "from the continuum data taken at sqrt(s) = 3.08 GeV: N_bg = N_3.08 * (L_J/psi / L_3.08) * " \
        "(s_3.08 / s_J/psi).  With N_3.08 = 1440376 +- 1200, this gives N_bg = (14.55 +- 0.02) x 10^6 " \
        "for 2012 (1.7% of the selected events) and N_bg = (6.58 +- 0.04) x 10^6 for 2009 " \
        "(total background 3.7%, including the ~1.5% QED contribution).")
  .note(:luminosity_from_gammagamma,
        "Integrated luminosities of the J/psi and of the sqrt(s) = 3.08 GeV samples are determined " \
        "from e+e- -> gamma gamma: at least two EMC showers, second most energetic shower energy in " \
        "[1.2, 1.6] GeV, both shower polar angles |cos theta| < 0.8; signal region " \
        "|Delta phi| < 2.5 deg with sideband 2.5 < |Delta phi| < 5 deg, where " \
        "Delta phi = |phi_gamma1 - phi_gamma2| - 180 deg.  For 2012 this gives " \
        "L_J/psi = 315.02 +- 0.14 pb^-1 and L_3.08 = 30.84 +- 0.04 pb^-1.")
  .note(:efficiency_from_psi2S,
        "The inclusive J/psi detection efficiency eps_data^{psi(3686)} is determined experimentally " \
        "from psi(3686) -> pi+ pi- J/psi, J/psi -> inclusive.  Tag-side soft-pion selection: at " \
        "least two oppositely charged pions in the MDC with |cos theta| < 0.93, Vr < 1 cm, " \
        "|Vz| < 15 cm and momentum < 0.4 GeV/c, with no further requirement on the remaining " \
        "tracks or showers.  The pi+ pi- recoil mass spectrum is fitted with a double Gaussian " \
        "for the J/psi signal plus a second-order Chebychev polynomial for the background, giving " \
        "N_inc = (1147.8 +- 1.9) x 10^3; after applying the same inclusive J/psi selection to the " \
        "remaining tracks and showers, N_inc^sel = (877.6 +- 1.7) x 10^3, hence " \
        "eps_data^{psi(3686)} = (76.46 +- 0.07)%.")
  .note(:correction_factor,
        "Correction factor for the kinematic difference between J/psi at rest and J/psi from " \
        "psi(3686) -> pi+ pi- J/psi: f_cor = eps_MC^{J/psi} / eps_MC^{psi(3686)}.  The two large " \
        "inclusive MC samples give eps_MC^{psi(3686)} = (75.76 +- 0.06)% and " \
        "eps_MC^{J/psi} = (76.58 +- 0.04)%, i.e. f_cor = 1.0109 +- 0.0009 (2012) and " \
        "f_cor = 1.0105 +- 0.0009 (2009).")
  .note(:trigger_efficiency,
        "Trigger efficiency eps_trig = 1.00 (100%), taken from the BESIII trigger study of " \
        "various reactions.")
  .note(:mc_model_uncertainty,
        "The MC model dependence is evaluated by generating alternative inclusive samples without " \
        "the LUNDCHARM model and comparing the resulting correction factor with its nominal value; " \
        "the changes 0.42% (2012) / 0.36% (2009) are assigned as systematic uncertainties.")
  .note(:tracking_uncertainty,
        "The track reconstruction efficiency difference between data and MC is < 1% per charged " \
        "track.  Varying the tracking efficiency by -1% in both J/psi and psi(3686) MC samples " \
        "changes the correction factor by 0.03%; an additional uncertainty of 0.30% is assigned to " \
        "the 2009 sample where the J/psi and psi(3686) data were taken at different times.")
  .note(:fit_to_jpsi_peak,
        "The fit to the pi+ pi- recoil mass spectrum contributes 0.19% (2012) / 0.17% (2009), " \
        "obtained by varying the fit range ([3.07, 3.13] -> [3.08, 3.12] GeV/c^2), the signal " \
        "shape (histogram from psi(3686) -> pi+ pi- J/psi, J/psi -> mu+ mu-) and the background " \
        "shape (second- -> first-order Chebychev polynomial).")
  .note(:noise_mixing_uncertainty,
        "Noise-realisation uncertainty estimated by reconstructing the psi(3686) MC sample with the " \
        "noise sample accompanying the J/psi data taking: 0.09% (2012) and 0.12% (2009, including " \
        "the 0.06% from splitting the data into three subsamples of different noise level).")
  .note(:soft_pion_efficiency_uncertainty,
        "The soft-pion selection efficiency eps_{pi+pi-} recoiling against the J/psi in " \
        "psi(3686) -> pi+ pi- J/psi depends on the multiplicity of the J/psi decay.  The MC " \
        "efficiency is reweighted with the multiplicity dependence extracted from the data " \
        "(psi(3686) sample versus J/psi at rest), giving 0.28% (2012) and 0.34% (2009).")

jpsi_algorithm.with_decay_card(decay_card_jpsi_inclusive).apply(event_selection)

root_files = jpsi_algorithm.execute_on([jpsi_data, jpsi_incMC,
                                        psip_data, psip_incMC,
                                        exMC_jpsi_inclusive, exMC_psip_to_pipi_jpsi])
