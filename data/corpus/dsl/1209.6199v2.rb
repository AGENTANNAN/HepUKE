# =============================================================================
# BESIII paper 1209.6199v2 — Determination of the number of psi(2S) events
#   collected with BESIII in 2009, using psi' -> inclusive hadrons.
# arXiv:1209.6199v2
#
# This is an inclusive-counting analysis: the signal channel is
# psi' -> hadrons, there is no exclusive final state, no PID and no kinematic
# fit. The BOSS-side scope covers the track-level and event-level selection.
# The 3.650 GeV off-resonance continuum data are used to subtract the
# non-resonant background.
# =============================================================================

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(2S) data, 2009 run
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # inclusive psi(2S) MC (efficiency)
off_data   = DatasetManager.real_data.find("709_3650")     # off-resonance data, 3.650 GeV (44 pb^-1)
off_incMC  = DatasetManager.inclusive_mc.find("709_3650")  # continuum inclusive MC

# Inclusive psi(2S) -> hadrons sample for the detection efficiency
decay_card_psip_incl = <<~DECAYCARD
    Decay psi(2S)
    1.0000 anything                PHSP;
    Enddecay

    End
DECAYCARD

exMC_psip_incl = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_inclusive_hadrons"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_psip_incl
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsipInclusiveCount"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.686]})
   .set_alias({"std::vector<double>" => "Vdouble"})

# Track-level selection:
#   charged tracks must pass within 1 cm of the beam line in the plane
#   perpendicular to the beam and within +-15 cm of the IP along the beam;
#   at least one good charged track is required at the event level (N_good >= 1).
# Photon selection: EMC barrel (|cos theta| < 0.8) E > 25 MeV, end-cap
#   (0.86 < |cos theta| < 0.92) E > 50 MeV; showers in the intermediate angular
#   region are poorly reconstructed and excluded.
sel = Selection.new
sel.select_track {
      cos_theta 0.93    # |cos theta| < 0.93 (MDC acceptance)
      Vz        15.0    # |Vz| < 15 cm from the IP along the beam direction
      Vr         1.0    # |Vr| < 1 cm from the beam line in the transverse plane
      nTot      ">=1"   # at least one good charged track (N_good >= 1)
    }
   .select_photon {
      energyThreshold_b 0.025   # barrel (|cos theta| < 0.8): E > 25 MeV
      energyThreshold_e 0.050   # end-cap (0.86 < |cos theta| < 0.92): E > 50 MeV
      angle_to_track    10.0
      tdc_emc_start     0       # EMC cluster timing window suppresses noise and
      tdc_emc_end       14      # energy deposits unrelated to the event
      nGam              ">=2"   # at least two photons (needed for the N_good = 1 pi0)
    }

alg.note(:event_level_mult_gt2,
         "for events with N_good > 2 charged tracks no additional selection is applied; " \
         "the event passes on the inclusive-hadrons topology.")
    .note(:two_track_bhabha_dimuon_veto,
         "for events with exactly two good charged tracks (N_good = 2), where Bhabha and " \
         "dimuon events dominate, the momentum of each track is required to be less than " \
         "1.7 GeV/c and the opening angle between the two tracks is required to be less than " \
         "176 degrees. These event-level requirements are applied in ROOT on the stored " \
         "track momenta and directions.")
    .note(:single_track_pi0_requirement,
         "for events with one good charged track (N_good = 1), at least two additional photons " \
         "are required; from all photon-pair combinations the one whose invariant mass M(gamma " \
         "gamma) is closest to the pi0 mass is selected and |M(gamma gamma) - M(pi0)| < 0.015 " \
         "GeV/c^2 is required. Applied in ROOT via the pi0 reconstruction from the photon list.")
    .note(:visible_energy_cut,
         "E_visible/E_cm > 0.4 is required to suppress the low-energy background (LEB, mostly " \
         "e+e- -> e+e- + X and double-ISR events). E_visible is the energy sum of all charged " \
         "tracks (computed with the track momentum under the pion-mass hypothesis) plus neutral " \
         "EMC showers. Applied in ROOT.")
    .note(:average_z_vertex_signal_sideband,
         "the average Z-direction vertex of the event is defined as Vbar_Z = " \
         "sum_i V_Z^i / N_good, where V_Z^i is the distance of the point of closest approach of " \
         "track i to the IP along the beam direction. Events with |Vbar_Z| < 4.0 cm are taken as " \
         "signal, while events with 6.0 cm < |Vbar_Z| < 10.0 cm form the non-collision background " \
         "sideband, giving N_obs = N_signal - N_sideband. Applied in ROOT.")
    .note(:background_subtraction,
         "the off-resonance data sample at E_cm = 3.650 GeV is used to subtract the continuum " \
         "QED and radiative-return J/psi backgrounds, with the scaling factor " \
         "f = (L_psi'/L_3.65) * (3.65^2/3.686^2) = 3.677, determined from the integrated " \
         "luminosities and the 1/s energy dependence of the cross-section. The luminosities are " \
         "measured from e+e- -> gamma gamma events. The final number of psi' events is " \
         "N_psi' = (N_peak^obs - f * N_off-resonance^obs) / epsilon, where epsilon is the " \
         "selection efficiency from the inclusive psi' MC sample (the B(psi' -> inclusive " \
         "hadrons) is folded into the efficiency). Applied in ROOT.")
    .note(:luminosity_from_gammagamma,
         "the integrated luminosities at the psi' peak and at 3.650 GeV are determined from " \
         "e+e- -> gamma gamma events with the same track- and event-level selection: no good " \
         "charged track and at least two showers, the most energetic shower E > 0.7 E_beam and " \
         "the second most energetic shower E > 0.4 E_beam, and the two most energetic showers " \
         "back-to-back in the psi' rest frame with 178 < |phi_1 - phi_2| < 182 degrees. The " \
         "luminosity systematic errors nearly cancel in the ratio; the alternative determination " \
         "with Bhabha events gives f = 3.685.")
    .note(:leb_check,
         "the validity of the continuum subtraction for the residual low-energy background is " \
         "checked with candidate LEB events (E_visible/E_cm < 0.35): the peak/off-resonance ratios " \
         "are 3.3752 (N_good = 1) and 3.652 (N_good = 2), to be compared with the scaling factor " \
         "f = 3.677; the ~10% difference for N_good = 1 is taken as a systematic error.")
    .note(:efficiency_determination,
         "the detection efficiency epsilon for inclusive psi' -> hadrons is determined from " \
         "106 x 10^6 psi' inclusive MC events after the same selection criteria. Typical values " \
         "are 92.912%, 89.860%, 74.624% and 58.188% for N_good >= 1, 2, 3 and 4, respectively.")
    .note(:systematics,
         "the total systematic error is 0.81%, obtained by summing in quadrature the following " \
         "sources: background contamination 0.10%, N_obs determination (counting vs. Vbar_Z fit) " \
         "0.28%, choice of sideband region 0.45%, vertex selection 0.35%, momentum and opening " \
         "angle 0.05%, scaling factor f 0.02%, 0-prong events 0.17%, tracking 0.03%, charged-track " \
         "multiplicity (unfolding) 0.40%, sigma(e+e- -> tau+tau-) 0.17%, B(psi' -> X + J/psi) " \
         "0.00%, pi0 mass requirement 0.11%, event start time determination 0.10%, trigger " \
         "efficiency (negligible) and B(psi' -> hadrons) 0.13%. The final result is " \
         "N_psi' = (106.41 +- 0.86) x 10^6, where the error is systematic only.")

alg.with_decay_card(decay_card_psip_incl).apply(sel)
root_files = alg.execute_on([psip_data, psip_incMC, off_data, off_incMC, exMC_psip_incl])
