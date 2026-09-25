# ============================================================================
# arXiv:1301.0053v1
# "Partial wave analysis of J/psi -> gamma eta eta"
# BESIII, 2.25e8 J/psi events at sqrt(s) = 3.097 GeV
#
# Signal chain : J/psi -> gamma eta eta, eta -> gamma gamma (5 photons, no tracks)
# The PWA itself (GPUPWA, covariant tensor amplitudes, unbinned maximum
# likelihood) is a ROOT-level procedure and lies outside the BOSS scope.
# ============================================================================

### ---------------------------- Dataset preparation ---------------------------- ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Inclusive J/psi MC (background study)

# Signal process: J/psi -> gamma eta eta (the scalar/tensor resonances X are
# not fixed by the generator; a phase-space generation is used for efficiency)
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta eta PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi2gammaetaeta"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Background channel explicitly identified in the paper (studied with MC)
decay_card_bkg_etapi0pi0 = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta pi0 pi0 PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_bkg_etapi0pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi2gammaetapi0pi0"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg_etapi0pi0
  config.cross_section   = :default
end

### ---------------------------- Event selection (BOSS) ---------------------------- ###
alg_name = "Jpsigammaetaeta"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.097]})

sel = Selection.new
sel.select_track {                    # no charged tracks in the final state
      cos_theta 0.92
      Vz        10.0
      Vr        1.0
      nChrp     "==0"
      nChrn     "==0"
    }
   .select_photon {                   # five (or six) good photons
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025         # barrel  E > 25 MeV
      energyThreshold_e 0.050         # endcap  E > 50 MeV
      nGam              ">=5"
    }
   # Nominal 4C kinematic fit under the J/psi -> 5 gamma hypothesis. When more
   # than five photons are selected, all permutations are iterated and the
   # combination with the smallest chi2_4C is kept (automatic in kinematic_fit).
   # The paper requires chi2_4C < 50; the loose default is used here and the
   # optimal value applied in ROOT.
   .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
    }

alg.note(:background_veto, "pi0 veto: events with any photon pair satisfying
    |M(gamma gamma) - m_pi0| < 0.015 GeV/c^2 are removed to suppress
    J/psi -> gamma pi0 pi0; the veto is applied on all photon-pair combinations
    before the 4C kinematic fit (no dedicated DSL primitive for an all-pairs veto)")
  .note(:signal_combination, "the two eta candidates are formed from the photon
    pairing that minimises delta = sqrt((M(gamma1 gamma2) - m_eta)^2 +
    (M(gamma3 gamma4) - m_eta)^2); the pairing is fixed here / in ROOT since
    nested for_each combination loops are experimental")
  .note(:mass_window, "both photon-pair invariant masses must satisfy
    |M(gamma gamma) - m_eta| < 0.04 GeV/c^2 (eta mass region, resolution about
    10 MeV/c^2); applied in ROOT after the nominal kinematic fit")
  .note(:background_veto, "mis-combination suppression: among all photon
    pairings only one combination must satisfy delta < 0.05 GeV/c^2, which
    reduces the fraction of events with a mis-combined photon from 5.3% to 0.8%")
  .note(:background_veto, "phi veto: |M(gamma eta) - m_phi| > 30 MeV/c^2 removes
    background from J/psi -> phi eta with phi -> gamma eta")
  .note(:background_veto, "background level of 6% estimated from the two
    dimensional eta mass sidebands 0.07 GeV/c^2 < |M(gamma gamma) - m_eta| <
    0.15 GeV/c^2; sideband events are included in the PWA with negative weights")
  .note(:efficiency_curve, "the detection efficiency is determined from a PWA
    amplitude-weighted MC sample (GPUPWA / covariant tensor amplitudes), not from
    a flat phase-space sample")

alg.with_decay_card(decay_card_signal).apply(sel)

### --------------------------------- Execution --------------------------------- ###
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_bkg_etapi0pi0])
