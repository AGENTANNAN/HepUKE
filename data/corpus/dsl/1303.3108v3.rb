# ============================================================================
# arXiv:1303.3108v3
# "Study of J/psi -> omega p pbar with omega -> gamma pi0"
# BESIII, 225.3e6 J/psi events at sqrt(s) = 3.097 GeV
#
# Signal chain : J/psi -> omega p pbar, omega -> gamma pi0, pi0 -> gamma gamma
# Final state  : gamma gamma gamma p pbar  (3 photons, 2 charged tracks)
# ============================================================================

### ---------------------------- Dataset preparation ---------------------------- ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Inclusive J/psi MC

# Signal: J/psi -> omega p pbar, omega -> gamma pi0, pi0 -> gamma gamma.
# Phase-space generation is used for the efficiency, as in the paper.
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 omega p+ anti-p- PHSP;
  Enddecay

  Decay omega
  1.0000 gamma pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi2omegappbar_omega2gammapi0"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Irreducible background: the phase-space process J/psi -> omega p pbar x
# gamma pi0 p pbar (same final state, no X(ppbar) resonant sub-structure);
# its ppbar mass shape is used as the non-resonant template in the fit.
decay_card_bkg_phsp = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma pi0 p+ anti-p- PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_bkg_phsp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi2gammapi0ppbar_phsp"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_phsp
  config.cross_section   = :default
end

### ---------------------------- Event selection (BOSS) ---------------------------- ###
alg_name = "Jpsitoomegappbar"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.097]})

sel = Selection.new
sel.select_track {                    # two well reconstructed charged tracks, net charge zero
      # tracks are reconstructed over |cos(theta)| < 0.93 in the MDC but only the
      # barrel region (|cos(theta)| < 0.8) is used, to reduce the systematic
      # uncertainties in tracking and particle identification
      cos_theta 0.8
      Vz        10.0                  # within +-10 cm of the IP along the beam direction
      Vr        1.0                   # within 1 cm of the IP transverse to the beam
      nChrp     "==1"
      nChrn     "==1"
      nNet      "==0"
    }
   # TOF and dE/dx are combined into confidence levels for the pi, K and p(pbar)
   # hypotheses; the hypothesis with the highest confidence level is assigned.
   # One proton and one anti-proton are required.
   .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:pion, :kaon]
      nprp "==1"
      nprm "==1"
    }
   .select_photon {                   # at least three photons from omega -> gamma pi0
      tdc_emc_start     0
      tdc_emc_end       14            # EMC timing window against noise
      energyThreshold_b 0.025         # barrel  E > 25 MeV (|cos(theta)| < 0.8)
      energyThreshold_e 0.050         # endcap  E > 50 MeV (0.86 < |cos(theta)| < 0.92)
      angle_to_track    10.0          # photon at least 10 deg from the nearest track
      nGam              ">=3"
    }
   # additional suppression of showers produced by the proton / anti-proton:
   # the photon must be more than 30 deg away from each of them
   .select_isolated_photon {
      angle_to_prp_track 30.0
      angle_to_prm_track 30.0
      nGam               ">=3"
    }
   # to reduce the data/MC difference of the tracking efficiency at low momentum
   .remove(:prp) { condition "three_momentum_of(:prp) < 0.3" }
   .remove(:prm) { condition "three_momentum_of(:prm) < 0.3" }
   # Nominal 4C energy-momentum conserving fit under the gamma gamma gamma p pbar
   # hypothesis. When more than three photons are selected, all combinations are
   # iterated and the one with the minimum chi2_4C is kept (automatic in
   # kinematic_fit). The paper requires chi2_4C < 30; the tight value is optimised
   # in ROOT and applied there, so only the loose default is used here.
   .kinematic_fit([:gamma, :gamma, :gamma, :prp, :prm]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
    }

alg.note(:signal_combination, "the pi0 candidate is formed from the two of the
    three selected photons whose invariant mass is closest to the pi0 mass; the
    omega is then reconstructed from the remaining photon and the pi0. The
    two-photon choice has no dedicated DSL combinatorics primitive (nested
    for_each is experimental) and is applied in ROOT, together with the
    requirement |M(gamma gamma) - M_pi0| < 15 MeV/c^2")
  .note(:mass_window, "the omega signal region is 0.753 GeV/c^2 < M(gamma pi0) <
    0.813 GeV/c^2; the omega sidebands are 0.663 < M(gamma pi0) < 0.693 GeV/c^2
    and 0.873 < M(gamma pi0) < 0.903 GeV/c^2, used to constrain the non-omega
    background in the ppbar mass fits. Applied in ROOT")
  .note(:analysis, "the gamma pi0 invariant mass is fitted with a Breit-Wigner
    (omega) convolved with a Novosibirsk resolution function plus a second-order
    Chebychev background; the ppbar mass near threshold is fitted with an
    acceptance-weighted S-wave Breit-Wigner for X(ppbar), the non-omega
    background f(delta) = N(delta^(1/2) + a1 delta^(3/2) + a2 delta^(5/2)) with
    delta = M(ppbar) - 2 m_p, and the sideband-constrained non-omega ppbar
    contribution. The yield upper limit is obtained with a Bayesian approach")
  .note(:efficiency_curve, "the detection efficiency for the signal is
    (16.1 +- 1.7)%, obtained from a MC sample of J/psi -> omega p pbar events
    generated with a phase-space distribution; the X(ppbar) acceptance enters
    through the acceptance-weighted Breit-Wigner used in the fit")

alg.with_decay_card(decay_card_signal).apply(sel)

### --------------------------------- Execution --------------------------------- ###
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_bkg_phsp])
