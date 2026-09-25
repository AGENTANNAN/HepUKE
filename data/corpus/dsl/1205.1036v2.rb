# =============================================================================
# BESIII: Study of J/psi -> p pbar and J/psi -> n nbar
#         (branching fractions and angular distributions)
# arXiv:1205.1036v2
# Sample: 225.2 million J/psi events collected in 2009
# =============================================================================

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ---- Decay cards -------------------------------------------------------------
# Signal: J/psi -> p pbar (phase space)
decay_card_pp = <<~DECAYCARD
  Decay J/psi
  1.000 p+ anti-p-   PHSP;
  Enddecay

  End
DECAYCARD

# Signal: J/psi -> n nbar (all-neutral final state, phase space)
decay_card_nn = <<~DECAYCARD
  Decay J/psi
  1.000 n0 anti-n0   PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive background for the p pbar analysis: J/psi -> gamma eta_c, eta_c -> p pbar
decay_card_bkg_gammac_pp = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta_c   PHSP;
  Enddecay

  Decay eta_c
  1.000 p+ anti-p-    PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive background for the p pbar analysis: J/psi -> pi0 p pbar
decay_card_bkg_pi0pp = <<~DECAYCARD
  Decay J/psi
  1.000 pi0 p+ anti-p-   PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma      PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive background for the n nbar analysis: J/psi -> gamma eta_c, eta_c -> n nbar
decay_card_bkg_gammac_nn = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta_c   PHSP;
  Enddecay

  Decay eta_c
  1.000 n0 anti-n0    PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive background for the n nbar analysis: J/psi -> pi0 n nbar
decay_card_bkg_pi0nn = <<~DECAYCARD
  Decay J/psi
  1.000 pi0 n0 anti-n0   PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma      PHSP;
  Enddecay

  End
DECAYCARD

# ---- Exclusive MC samples ----------------------------------------------------
exMC_pp = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_ppbar"
  c.related_dataset = jpsi_data
  c.events          = 1000000
  c.decay_card      = decay_card_pp
  c.cross_section   = :default
end

exMC_nn = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_nnbar"
  c.related_dataset = jpsi_data
  c.events          = 1000000
  c.decay_card      = decay_card_nn
  c.cross_section   = :default
end

exMC_bkg_pp = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_gamma_etac_ppbar"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_bkg_gammac_pp
  c.cross_section   = :default
end

exMC_bkg_pi0pp = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_pi0_ppbar"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_bkg_pi0pp
  c.cross_section   = :default
end

exMC_bkg_nn = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_gamma_etac_nnbar"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_bkg_gammac_nn
  c.cross_section   = :default
end

exMC_bkg_pi0nn = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_pi0_nnbar"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_bkg_pi0nn
  c.cross_section   = :default
end

# =============================================================================
# ALGORITHM 1: J/psi -> p pbar
# =============================================================================
alg_pp = Algorithm.new("JpsiToPPbar")
alg_pp.set_header(["JpsiToPPbarAlg/JpsiToPPbar.h"])
      .set_constant({"ECMS" => [:double, 3.097]})   # J/psi center-of-mass energy (GeV)

sel_pp = Selection.new
sel_pp.select_track {
         # Exactly two good charged tracks with |cos(theta)| < 0.8 (the two endcap
         # regions are excluded to reduce systematic uncertainties in tracking and PID)
         cos_theta  0.80
         Vz         10.0   # within +-10 cm of the IP along the beam direction
         Vr         1.0    # within +-1 cm of the IP in the plane perpendicular to the beam
         nChrp      "==1"
         nChrn      "==1"
         nNet       "==0"
       }
       # Loose PID applied to the positive track only (probability of the p hypothesis
       # greater than the pion and kaon hypotheses); the negative track is NOT identified
       # in order to maximise the efficiency and minimise the systematic uncertainty.
       .pid(method: :probability) {
         prob_cut 0.001
         identify :prp, against: [:pion, :kaon]
       }
       .assign({chrgn: :prm})   # the negative track is treated as an antiproton without PID
       # Vertex fit to the two selected tracks, improving the momentum resolution.
       # The vertex-fit-corrected momenta are then used for the p pbar opening-angle and
       # momentum requirements, which are applied at the ROOT level.
       .kinematic_fit([:prp, :prm]) {
         nominal
         vertex_fit([0, 1])
         chi2_cut 200   # loose default; no chi2 requirement is quoted in the paper
       }

alg_pp
  .note(:track_quality,
        "the track selection is restricted to |cos(theta)| < 0.8 (both endcap regions excluded) " \
        "to reduce systematic uncertainties in tracking and particle identification.")
  .note(:pid_asymmetry,
        "a loose PID requirement is applied only to the positive track (probability of the p " \
        "hypothesis greater than those of the pi+ and K+ hypotheses); no particle identification " \
        "is required for the negative track, which is taken to be the antiproton. This maximises " \
        "the efficiency and minimises the systematic uncertainty.")
  .note(:opening_angle,
        "the angle between the proton and the antiproton is required to be greater than 178 degrees; " \
        "this uses the vertex-fit-corrected momenta and is applied at the ROOT level.")
  .note(:momentum_window,
        "the measured momentum magnitude of both tracks is required to be within 30 MeV/c (about " \
        "3 sigma) of the expected value of 1.232 GeV/c for J/psi -> p pbar; applied at the ROOT level.")
  .note(:background_estimation,
        "backgrounds (estimated at 0.02% - 0.2% of the signal) are evaluated with three independent " \
        "procedures: the inclusive J/psi MC sample, exclusive MC of potential background processes " \
        "(Bhabha events, J/psi -> e+ e-, mu+ mu-, K+ K-, gamma p pbar, pi0 p pbar and gamma eta_c " \
        "with eta_c -> p pbar), and a sideband technique. The largest estimate (sideband) is assigned " \
        "as a systematic uncertainty; no subtraction is applied.")
  .note(:efficiency_correction,
        "the MC-determined efficiency is corrected in 16 bins of cos(theta) using data/MC ratios of " \
        "the tracking efficiency (from J/psi -> p pbar events with 1 or 2 good tracks) and of the " \
        "particle-identification efficiency; the corrected efficiency is smoothed with a fifth-order " \
        "polynomial in cos(theta).")
  .note(:angular_distribution,
        "the cos(theta) distribution of the proton is fitted with A (1 + alpha cos^2 theta) eps(cos theta) " \
        "to extract the angular-distribution parameter alpha and the normalisation; the yield is then " \
        "extrapolated from |cos theta| < 0.8 to the full angular range (ROOT level).")
  .note(:continuum_interference,
        "the systematic effect of interference between the J/psi peak and the continuum p pbar " \
        "production is evaluated using the total cross-section expression with the continuum term " \
        "sigma_ppbar^cont. taken from the amplitude analysis; the difference with and without this " \
        "term is assigned as a systematic uncertainty.")
  .note(:trigger_efficiency,
        "trigger efficiency is close to 100% for hadronic events containing charged particles.")

alg_pp.with_decay_card(decay_card_pp).apply(sel_pp)
root_files_pp = alg_pp.execute_on([jpsi_data, jpsi_incMC, exMC_pp, exMC_bkg_pp, exMC_bkg_pi0pp])

# =============================================================================
# ALGORITHM 2: J/psi -> n nbar  (all-neutral final state, EMC-shower analysis)
# =============================================================================
alg_nn = Algorithm.new("JpsiToNNbar")
alg_nn.set_header(["JpsiToNNbarAlg/JpsiToNNbar.h"])
      .set_constant({"ECMS" => [:double, 3.097]})

sel_nn = Selection.new
sel_nn.select_track {
         # No good charged track originating in the interaction region is allowed
         nTot "==0"
       }
       .select_photon {
         # EMC showers are the only handle on this all-neutral final state: the
         # antineutron annihilation "star" is selected as the most energetic shower
         # (0.6 - 2.0 GeV) and the neutron as a shower on the opposite side (0.06 - 0.6 GeV)
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         angle_to_track    10.0
         tdc_emc_start     0
         tdc_emc_end       14
         nGam              ">=2"
       }
       # No kinematic fit is performed in this analysis: the antineutron and neutron
       # candidates are identified from the EMC shower topology and the event yields are
       # extracted from fits to the angle between the nbar and n in bins of cos(theta).
       .partial_miss([]) {
         # no additional recoil-mass constraint is applied
       }

alg_nn
  .note(:nbar_candidate,
        "the most energetic EMC shower in the event is assigned to be the antineutron candidate and " \
        "is required to have an energy in the range 0.6 - 2.0 GeV, with a fiducial cut of " \
        "|cos(theta)| < 0.8 ensuring that the energy is fully contained in the EMC.")
  .note(:shower_second_moment,
        "the nbar candidate shower must satisfy the second moment requirement " \
        "S = sum_i E_i r_i^2 / sum_i E_i > 20 cm^2, where E_i is the energy deposited in the i-th " \
        "crystal and r_i the distance from that crystal to the shower centre; this suppresses " \
        "photon backgrounds.")
  .note(:shower_topology,
        "the number of EMC hits in a 50 degree cone around the nbar candidate shower direction is " \
        "required to be greater than 40, exploiting the distinctive hadronic shower topology of " \
        "the antineutron annihilation star.")
  .note(:neutron_candidate,
        "EMC showers on the opposite side of the detector are searched for the neutron: the shower " \
        "energy must lie between 0.06 and 0.6 GeV; if several showers are present the one most " \
        "back-to-back with respect to the nbar candidate is selected.")
  .note(:extra_energy_veto,
        "E_extra = 0 is required, where E_extra is the total deposited energy in the EMC excluding " \
        "the nbar shower and any additional energy inside the 50 degree cone; this suppresses " \
        "all-neutral J/psi decays, continuum production and electromagnetic processes.")
  .note(:signal_extraction,
        "the J/psi -> n nbar signal appears as an enhancement near 180 degrees in the angle between " \
        "the nbar shower and the neutron direction; the number of events is obtained by fitting the " \
        "angle distribution in bins of cos(theta) with signal and background shapes (ROOT level).")
  .note(:efficiency_determination,
        "the efficiency is obtained from data using J/psi -> p nbar pi- (and c.c.) and " \
        "J/psi -> pbar n pi+ samples to determine the nbar and neutron selection efficiencies " \
        "separately, and the J/psi -> p pbar sample to determine the E_extra cut efficiency; the " \
        "product efficiency is computed in 16 bins of cos(theta) and smoothed with a fifth-order " \
        "polynomial.")
  .note(:background_estimation,
        "backgrounds are assessed with the generic inclusive J/psi MC sample and with exclusive " \
        "channels J/psi -> pi0 n nbar, J/psi -> gamma n nbar, e+e- -> gamma gamma, " \
        "J/psi -> Sigma+ Sigma-, J/psi -> Sigma+ anti-Sigma-, J/psi -> p pbar and J/psi -> gamma eta_c " \
        "(eta_c -> n nbar); none of these exhibits peaking in the nbar-n angle distribution.")
  .note(:signal_shape,
        "the signal shape for fitting the angle between nbar and n is derived from the J/psi -> p pbar " \
        "sample, with the p pbar angular distribution corrected for the 1.0 T magnetic field before " \
        "being applied to J/psi -> n nbar.")
  .note(:trigger_efficiency,
        "for the all-neutral n nbar final state the trigger efficiency is a potential source of " \
        "uncertainty; the efficiency curve is corrected with the MC-determined trigger efficiency.")
  .note(:continuum_interference,
        "the systematic effect from interference between the J/psi peak and the continuum n nbar " \
        "production is evaluated using the total cross-section expression with the continuum term " \
        "sigma_nnbar^cont.; the difference with and without this term is assigned as a systematic error.")

alg_nn.with_decay_card(decay_card_nn).apply(sel_nn)
root_files_nn = alg_nn.execute_on([jpsi_data, jpsi_incMC, exMC_nn, exMC_bkg_nn, exMC_bkg_pi0nn])
