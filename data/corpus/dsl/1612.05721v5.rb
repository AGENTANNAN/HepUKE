# ============================================================================
# BESIII: first study of the doubly radiative decay eta' -> gamma gamma pi0
# (eta' produced via J/psi -> gamma eta')
# arXiv:1612.05721v5,  1.31 x 10^9 J/psi events at sqrt(s) = 3.097 GeV
#
# Signal chain: J/psi -> gamma eta', eta' -> gamma gamma pi0, pi0 -> gamma gamma
#               -> exactly five photons in the final state, no charged tracks
# The inclusive eta' -> gamma gamma pi0 yield, the eta' -> gamma omega
# (omega -> gamma pi0) component and the nonresonant component are extracted
# from ROOT-level fits to M(gamma gamma pi0) and M(gamma pi0); those fits, the
# M(gamma gamma)^2 dependent partial widths and the 5C kinematic fit's final
# chi2_5C < 30 requirement are all outside the BOSS scope (Rule T3).
# ============================================================================

### ---------------------------- Dataset preparation ---------------------------- ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # inclusive J/psi MC (background study)

# ---------------------------------------------------------------------------
# Signal: J/psi -> gamma eta', eta' -> gamma gamma pi0 (pi0 -> gamma gamma).
# The doubly radiative decay is generated with a phase-space model here; the
# paper uses the VMD model (rho(1450) / omega(1650) exchange) for the
# nonresonant component and helicity amplitudes for eta' -> gamma omega(rho)
# with omega(rho) -> gamma pi0.
# ---------------------------------------------------------------------------
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'            PHSP;
  Enddecay

  Decay eta'
  1.0000 gamma gamma pi0       PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma           PHSP;
  Enddecay

  End
DECAYCARD

# Dominant intermediate process of the same final state:
# J/psi -> gamma eta', eta' -> gamma omega, omega -> gamma pi0
decay_card_signal_omega = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'            PHSP;
  Enddecay

  Decay eta'
  1.0000 gamma omega           PHSP;
  Enddecay

  Decay omega
  1.0000 gamma pi0             PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma           PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etap_etap2gammagammapi0"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_signal_omega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etap_etap2gammaomega"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal_omega
  config.cross_section   = :default
end

# ---------------------------------------------------------------------------
# Background class I: J/psi -> gamma eta' with eta' decaying to a different
# final state.  These accumulate near the lower side of the eta' signal region
# and are dominated by eta' -> pi0 pi0 eta (eta -> gamma gamma), eta' -> 3 pi0
# and eta' -> gamma gamma.
# ---------------------------------------------------------------------------
decay_card_bkg_2pi0eta = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'            PHSP;
  Enddecay

  Decay eta'
  1.0000 pi0 pi0 eta           PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma           PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma           PHSP;
  Enddecay

  End
DECAYCARD

decay_card_bkg_3pi0 = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'            PHSP;
  Enddecay

  Decay eta'
  1.0000 pi0 pi0 pi0           PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma           PHSP;
  Enddecay

  End
DECAYCARD

exMC_bkg_2pi0eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etap_etap2pi0pi0eta"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_2pi0eta
  config.cross_section   = :default
end

exMC_bkg_3pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etap_etap23pi0"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_3pi0
  config.cross_section   = :default
end

# ---------------------------------------------------------------------------
# Background class II: J/psi decays without an eta', e.g.
# J/psi -> gamma pi0 pi0 and J/psi -> omega eta (omega -> gamma pi0, eta -> gamma gamma),
# which give a smooth distribution under the eta' signal region.
# ---------------------------------------------------------------------------
decay_card_bkg_g2pi0 = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma pi0 pi0         PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma           PHSP;
  Enddecay

  End
DECAYCARD

decay_card_bkg_omegaeta = <<~DECAYCARD
  Decay J/psi
  1.0000 omega eta             PHSP;
  Enddecay

  Decay omega
  1.0000 gamma pi0             PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma           PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma           PHSP;
  Enddecay

  End
DECAYCARD

exMC_bkg_g2pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_pi0pi0"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_g2pi0
  config.cross_section   = :default
end

exMC_bkg_omegaeta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_omega_eta"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_omegaeta
  config.cross_section   = :default
end

### ---------------------------- Event selection (BOSS) ---------------------------- ###
# One signal channel only (J/psi -> gamma eta', eta' -> gamma gamma pi0), so a
# single Algorithm object is used.
alg_name = "JpsiGammaEtapTo2GammaPi0"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.097]})

sel = Selection.new
sel.select_track {                 # no charged particles are allowed in the event
      nTot  "==0"
    }
   .select_photon {                # exactly five photon candidates
      tdc_emc_start     0          # EMC timing suppresses electronic noise and
      tdc_emc_end       14         # showers unrelated to the event (0, 700) ns
      angle_to_track    10.0       # > 10 degrees from any charged track
      energyThreshold_b 0.025      # barrel (|cos(theta)| < 0.80):  E > 25 MeV
      energyThreshold_e 0.050      # end-cap (0.86 < |cos(theta)| < 0.92): E > 50 MeV
      nGam              "==5"      # exactly five photons: 1 radiative + 2 from eta'
                                   # + 2 from pi0 -> gamma gamma
    }
   # 5C kinematic fit to the gamma gamma gamma pi0 hypothesis: 4-momentum
   # conservation (the initial e+e- / J/psi four-momentum) plus the pi0 mass
   # constraint on one photon pair.  When more than one pi0 (photon-pair)
   # combination is possible, the combination with the smallest chi2_5C is
   # retained.  The paper requires chi2_5C < 30; the loose BOSS default is used
   # here and the tight value is applied in ROOT (Rule T3).
   .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 200
    }

alg.note(:vertex_origin,
         "Only events without charged particles are analysed; the average event vertex of each " \
         "run is assumed to be the origin for the selected candidates.")
   .note(:photon_background,
         "Photons in the barrel/end-cap gap region are poorly measured and excluded; the energy " \
         "deposited in the nearby TOF counters is included in the photon reconstruction to improve " \
         "the efficiency and energy resolution.")
   .note(:pi0_veto,
         "To suppress multi-pi0 backgrounds and remove miscombined pi0 candidates, an event is " \
         "vetoed if any two of the five selected photons other than the pair forming the pi0 " \
         "candidate satisfy |M(gamma gamma) - M_pi0| < 18 MeV/c^2.  This is an all-pairs veto " \
         "applied after the 5C fit and is not expressible with a dedicated DSL primitive.")
   .note(:pi0_candidate_selection,
         "The pi0 candidate is reconstructed from a pair of photons and, if more than one pi0 " \
         "candidate exists, the combination with the smallest chi2_5C is selected - this " \
         "corresponds to the default minimum-chi2 combination selection performed inside the " \
         "kinematic fit block.")
   .note(:etap_reconstruction,
         "After the selection requirements the most energetic photon is taken as the primary " \
         "radiative photon from J/psi -> gamma eta', and the remaining two photons together with " \
         "the pi0 are used to reconstruct the eta' candidate; the eta' signal region is " \
         "|M(gamma gamma pi0) - M_eta'| < 25 MeV/c^2.  The eta' -> gamma omega component and the " \
         "nonresonant component are separated by a fit to M(gamma pi0) in the range " \
         "0.20-0.92 GeV/c^2, with the rho-omega interference considered.  ROOT-level procedure.")
   .note(:background_veto,
         "Class I background (J/psi -> gamma eta' with eta' decaying to a non-signal final state, " \
         "dominated by eta' -> pi0 pi0 eta, eta' -> 3 pi0 and eta' -> gamma gamma) accumulates near " \
         "the lower side of the eta' signal region and is fixed in shape and yield from MC.  " \
         "Class II background (J/psi decays without an eta', e.g. J/psi -> gamma pi0 pi0 and " \
         "J/psi -> omega eta) gives a smooth distribution under the eta' signal region and is " \
         "described by a third-order Chebyshev polynomial with free parameters in the fit to " \
         "M(gamma gamma pi0); for the M(gamma pi0) fit its shape is taken from the eta' mass " \
         "sidebands (738-788 and 1008-1058 MeV/c^2) with a fixed normalization.  No peaking " \
         "background remains after all selection criteria according to MC studies.")
   .note(:signal_mc_model,
         "The J/psi -> gamma eta' MC sample is generated with an angular distribution of " \
         "1 + cos^2(theta_gamma) for the radiative photon; eta' -> gamma omega(rho) with " \
         "omega(rho) -> gamma pi0 is generated with the helicity amplitude formalism, and the " \
         "nonresonant eta' -> gamma gamma pi0 is generated with the VMD model with rho(1450) or " \
         "omega(1650) exchange.  The inclusive signal shape used in the fit is an incoherent " \
         "mixture of the rho, omega and nonresonant components according to the fitted fractions, " \
         "and the nonresonant shape uncertainty is assessed with a phase-space alternative.")
   .note(:partial_widths,
         "The M(gamma gamma)^2 dependent partial widths dGamma(eta' -> gamma gamma pi0)/dM(gamma " \
         "gamma)^2 are obtained by bin-by-bin fits to the gamma gamma pi0 invariant mass " \
         "distribution (bin width 1.0 x 10^4 (MeV/c^2)^2, well above the 5 x 10^2 (MeV/c^2)^2 " \
         "resolution, so no unfolding is applied) and the efficiency-corrected, " \
         "background-subtracted yields.  ROOT-level procedure.")
   .note(:branching_fractions,
         "B(eta' -> gamma gamma pi0)_Incl. = N / (N_J/psi x epsilon x B(J/psi -> gamma eta') x " \
         "B(pi0 -> gamma gamma)) = (3.20 +- 0.07 +- 0.23) x 10^-3 with N_J/psi = (1310.6 +- 10.5) " \
         "x 10^6, epsilon = 16.1%; B(eta' -> gamma omega) x B(omega -> gamma pi0) = " \
         "(23.7 +- 1.4 +- 1.8) x 10^-4 with epsilon = 14.8%; the nonresonant " \
         "B(eta' -> gamma gamma pi0)_NR = (6.16 +- 0.64 +- 0.67) x 10^-4 with epsilon = 15.9%.")
   .with_decay_card(decay_card_signal)
   .apply(sel)

### --------------------------------- Execution --------------------------------- ###
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_signal_omega,
                exMC_bkg_2pi0eta, exMC_bkg_3pi0, exMC_bkg_g2pi0, exMC_bkg_omegaeta])
