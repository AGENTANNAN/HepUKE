# BESIII arXiv:1507.06734v1
# Observation of eta' -> omega e+ e- and measurement of B(eta' -> omega gamma) in
# J/psi -> gamma eta' with (1310.6 +/- 10.5) x 10^6 J/psi events at sqrt(s) = 3.097 GeV.
#
# Two signal modes (Rule T1), each with its own Algorithm object:
#   (A) eta' -> omega gamma,     omega -> pi+ pi- pi0, pi0 -> gamma gamma
#       final state: gamma gamma gamma pi+ pi-      (>=4 photon candidates, 2 tracks)
#   (B) eta' -> omega e+ e-,     omega -> pi+ pi- pi0, pi0 -> gamma gamma
#       final state: gamma gamma gamma pi+ pi- e+ e- (>=3 photons, 4 tracks)
# The J/psi is produced with kkmc (ISR), the decays are generated with evtgen
# (FSR for charged particles with PHOTOS); the inclusive J/psi sample models the
# known decays with PDG branching fractions and the remainder with lundcharm.

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # 1.31 x 10^9 J/psi events
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # inclusive J/psi MC (background study)

### Decay cards (EvtGen format) ###
# (A) J/psi -> gamma eta', eta' -> omega gamma, omega -> pi+ pi- pi0 (OMEGA_DALITZ)
decay_card_omega_gamma = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'                  PHSP;
    Enddecay

    Decay eta'
    1.0000 omega gamma                 PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0                 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                 PHSP;
    Enddecay

    End
DECAYCARD

# (B) eta' -> omega e+ e-; the Dalitz decay is generated with the VMD / chiral
# transition amplitude of the paper; no dedicated EvtGen model exists for it, so
# the phase-space generator is used (the form-factor dependence is a systematic).
decay_card_omega_ee = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta'                  PHSP;
    Enddecay

    Decay eta'
    1.0000 omega e+ e-                 PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0                 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                 PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples ###
exMC_omega_gamma = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "etap_omega_gamma"
  config.related_dataset = jpsi_data
  config.events          = 200000
  config.decay_card      = decay_card_omega_gamma
  config.cross_section   = :default
end

# 6.0 x 10^5 signal events for the omega e+ e- mode (paper).
exMC_omega_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "etap_omega_ee"
  config.related_dataset = jpsi_data
  config.events          = 600000
  config.decay_card      = decay_card_omega_ee
  config.cross_section   = :default
end

### ---------------------------------------------------------------------------
### Event selection (BOSS) — Mode A: eta' -> omega gamma
### ---------------------------------------------------------------------------
alg_name_A = "EtapOmegaGamma"
alg_A = Algorithm.new(alg_name_A)
alg_A.set_header(["#{alg_name_A}Alg/#{alg_name_A}.h"])
     .set_constant({"ECMS" => [:double, 3.097]})   # sqrt(s) = 3.097 GeV (J/psi)
     .set_alias({"std::vector<double>" => "Vdouble"})

sel_A = Selection.new
sel_A.select_track {
        cos_theta 0.93   # |cos(theta)| < 0.93, theta w.r.t. the beam direction
        Vz        20.0   # closest approach to the IP within +-20 cm along the beam
        Vr        2.0    # closest approach within 2 cm in the plane perpendicular to the beam
        nChrp     "==1"  # pi+ from the omega -> pi+ pi- pi0
        nChrn     "==1"  # pi- from the omega
        nNet      "==0"
      }
      .select_photon {
        tdc_emc_start     0     # EMC cluster time inside the 700 ns window around the
        tdc_emc_end       14    # event start time (suppresses electronic noise)
        angle_to_track    10.0  # shower at least 10 degrees from the nearest charged track,
                                # rejecting bremsstrahlung showers
        energyThreshold_b 0.025 # E > 25 MeV for barrel showers (|cos(theta)| < 0.80)
        energyThreshold_e 0.050 # E > 50 MeV for end-cap showers (0.86 < |cos(theta)| < 0.92)
        nGam              ">=4" # radiative photon + the two photons from pi0 -> gamma gamma
      }
      # No PID is applied in this channel: the two oppositely charged tracks are taken
      # directly as the pi+ and pi- from the omega.
      .assign({:chrgp => :pip, :chrgn => :pim})
      # pi0 -> gamma gamma: the photon pair is constrained to the nominal pi0 mass. The
      # pair with the two-photon mass closest to m(pi0) is chosen by the fit.
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      }
      # 4C kinematic fit under the J/psi -> gamma gamma gamma pi+ pi- hypothesis
      # (one radiative photon, the pi0 from two photons and the two pions); the
      # combination with the smallest chi2_4C is retained.
      .kinematic_fit([:pi0, :gamma, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }

alg_A
  .note(:radiative_photon_selection,
        "Since J/psi -> gamma eta' is a two-body decay, the radiative photon carries a fixed " \
        "energy of 1.4 GeV; the photon candidate with the maximum energy is taken as the " \
        "radiative photon and is required to have E > 1.0 GeV. Applied at ROOT level on the " \
        "selected photon candidates.")
  .note(:pi0_selection,
        "The photon pair whose two-photon invariant mass is closest to the nominal pi0 mass is " \
        "taken as the pi0 candidate, and its invariant mass must satisfy " \
        "|M(gamma gamma) - m(pi0)| < 0.015 GeV/c^2. The nominal-pi0 Kalman constraint is applied " \
        "in BOSS; the 0.015 GeV/c^2 window is applied at ROOT level.")
  .note(:chi2_4c_cut,
        "The paper retains only candidates with chi2_4C < 80. A loose chi2_cut of 200 is applied " \
        "in BOSS and the published value is applied on the stored chi2_4C in the ROOT analysis.")
  .note(:signal_extraction,
        "The eta' -> omega gamma yield is determined from the mass difference " \
        "M(pi0 pi+ pi- gamma) - M(pi0 pi+ pi-), fitted with an unbinned maximum-likelihood fit " \
        "in which the signal is described by the MC shape convoluted with a Gaussian and the " \
        "background by a 3rd-order Chebyshev polynomial. Performed at ROOT level.")
  .note(:efficiency_curve,
        "The detection efficiency is (21.87 +/- 0.02)%, obtained from the signal MC sample after " \
        "applying the full selection; the branching fraction is extracted after subtraction of " \
        "the non-peaking background in the vertical (omega) and horizontal (eta') bands.")
  .with_decay_card(decay_card_omega_gamma)
  .apply(sel_A)

### ---------------------------------------------------------------------------
### Event selection (BOSS) — Mode B: eta' -> omega e+ e-
### ---------------------------------------------------------------------------
alg_name_B = "EtapOmegaEE"
alg_B = Algorithm.new(alg_name_B)
alg_B.set_header(["#{alg_name_B}Alg/#{alg_name_B}.h"])
     .set_constant({"ECMS" => [:double, 3.097]})
     .set_alias({"std::vector<double>" => "Vdouble"})

sel_B = Selection.new
sel_B.select_track {
        cos_theta 0.93
        Vz        20.0
        Vr        2.0
        nChrp     "==2"  # pi+ from the omega and the e+
        nChrn     "==2"  # pi- from the omega and the e-
        nNet      "==0"
      }
      # Combinatorial PID: for every charged track the pion / electron / muon hypothesis is
      # tested with the combined TOF and dE/dx chi^2; the particle-type assignment that
      # minimises the total chi2_PID is adopted. Together with the 4C fit this forms
      # chi2_4C+PID = chi2_4C + sum_j chi2_PID(j), used to select the event and the best
      # photon combination.
      .pid(method: :chi2_sum) {
        chi_min_cut 4
        identify :pion, :electron
      }
      .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=3"  # radiative photon + the two photons from pi0 -> gamma gamma
      }
      # pi0 -> gamma gamma with a nominal-mass constraint; of the three photon candidates the
      # pair closest to m(pi0) is taken as the pi0 and the remaining one as the radiative photon.
      .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
      }
      # 4C kinematic fit under the J/psi -> gamma gamma gamma pi+ pi- e+ e- hypothesis. The
      # event is kept only if the smallest chi2_4C+PID assignment corresponds to two
      # oppositely charged pions plus an electron and a positron.
      .kinematic_fit([:pi0, :gamma, :ep, :em, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
      }

alg_B
  .note(:pid_correction_method,
        "The particle type of each charged track (pion, electron or muon) is assigned by " \
        "minimising chi2_4C+PID = chi2_4C + sum_j chi2_PID(j), where chi2_PID(j) combines the " \
        "TOF and dE/dx information of track j. The combinatorial assignment is expressed in BOSS " \
        "by pid(method: :chi2_sum); the joint minimisation together with the 4C chi^2 (i.e. the " \
        "use of the smallest chi2_4C+PID combination as the event selection) is applied at ROOT level.")
  .note(:radiative_photon_selection,
        "The selected photon with the maximum energy is taken as the radiative photon from " \
        "eta' -> omega e+ e- and is required to have E > 1.0 GeV. Applied at ROOT level.")
  .note(:pi0_selection,
        "The other two photons are required to be consistent with a pi0 candidate, " \
        "|M(gamma gamma) - m(pi0)| < 0.015 GeV/c^2. The nominal-pi0 Kalman constraint is applied " \
        "in BOSS; the 0.015 GeV/c^2 window is applied at ROOT level.")
  .note(:background_veto,
        "The peaking background from J/psi -> gamma eta', eta' -> omega gamma in which the photon " \
        "from the eta' converts into an electron-positron pair is removed by requiring the " \
        "distance of the reconstructed e+ e- vertex from the z axis, R_xy, to be less than 2 cm. " \
        "This gamma-conversion veto is applied at ROOT level.")
  .note(:background_veto,
        "The residual peaking background from eta' -> omega gamma with gamma -> e+ e- that " \
        "survives R_xy < 2 cm is estimated from MC to be 2.6 +/- 0.3 events, and is subtracted in " \
        "the branching-fraction calculation.")
  .note(:background_veto,
        "The dominant non-resonant background eta' -> pi+ pi- eta, eta -> pi0 pi+ pi- with the " \
        "pion pair misidentified as an electron-positron pair is described by the inclusive-MC " \
        "shape with its magnitude fixed to the PDG branching fraction.")
  .note(:chi2_4c_cut,
        "The paper retains only candidates with chi2_4C < 80. A loose chi2_cut of 200 is applied " \
        "in BOSS and the published value is applied on the stored chi2_4C in the ROOT analysis.")
  .note(:form_factor,
        "The nominal signal MC is generated with the VMD / chiral-perturbation-theory amplitude of " \
        "the paper; alternative monopole and dipole form factors give a maximum efficiency " \
        "variation of 1.3%, taken as a systematic uncertainty.")
  .note(:signal_extraction,
        "The eta' -> omega e+ e- yield is determined from an unbinned maximum-likelihood fit to " \
        "M(pi0 pi+ pi- e+ e-) - M(pi0 pi+ pi-), with the signal described by the MC shape " \
        "convoluted with a Gaussian and the remaining background by a 2nd-order Chebyshev " \
        "polynomial. Performed at ROOT level.")
  .note(:efficiency_curve,
        "The detection efficiency is (5.45 +/- 0.03)%, obtained from the VMD-based signal MC.")
  .with_decay_card(decay_card_omega_ee)
  .apply(sel_B)

### Execute on data, inclusive MC and the signal MC samples ###
alg_A.execute_on([jpsi_data, jpsi_incMC, exMC_omega_gamma])
alg_B.execute_on([jpsi_data, jpsi_incMC, exMC_omega_ee])
