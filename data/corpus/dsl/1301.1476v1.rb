# ============================================================================
# arXiv:1301.1476v1
# "Search for the M1 radiative transition psi(3686) -> gamma eta_c(2S) with
#  eta_c(2S) -> K_S0 K+- pi-+ pi+ pi-"
# BESIII, 1.06e8 psi(3686) events at sqrt(s) = 3.686 GeV
#
# Signal chain : psi(3686) -> gamma eta_c(2S),
#                eta_c(2S) -> K_S0 K+- pi-+ pi+ pi-, K_S0 -> pi+ pi-
# Final state  : gamma pi+ pi- K+- pi-+ pi+ pi-  (6 charged tracks, net charge 0)
# ============================================================================

### ---------------------------- Dataset preparation ---------------------------- ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(3686) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # Inclusive psi(3686) MC
cont_data  = DatasetManager.real_data.find("709_3650")      # 3.65 GeV continuum (non-resonant background)

# Signal: psi(3686) -> gamma eta_c(2S), eta_c(2S) -> K_S0 K+- pi-+ pi+ pi-
# eta_c(2S) is neutral (C = +), so the two charge-conjugate K/pi arrangements
# are listed explicitly.
decay_card_signal = <<~DECAYCARD
  Decay psi(3686)
  1.0000 gamma eta_c(2S) PHSP;
  Enddecay

  Decay eta_c(2S)
  0.5000 K_S0 K+ pi- pi+ pi- PHSP;
  0.5000 K_S0 K- pi+ pi- pi+ PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip2gammaetac2S_KsKpipipi"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Irreducible background: the phase-space process has the same final state as
# the signal (branching fraction calculated to be 1.73e-4 in the paper) and its
# K_S0 K 3pi mass shape is used as a template in the final fit.
decay_card_bkg_phsp = <<~DECAYCARD
  Decay psi(3686)
  0.5000 gamma K_S0 K+ pi- pi+ pi- PHSP;
  0.5000 gamma K_S0 K- pi+ pi- pi+ PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

exMC_bkg_phsp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip2gammaKsKpipipi_phsp"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_bkg_phsp
  config.cross_section   = :default
end

# Background psi(3686) -> pi0 K_S0 K+- pi-+ pi+ pi- (one pi0 photon missed)
decay_card_bkg_pi0 = <<~DECAYCARD
  Decay psi(3686)
  0.5000 pi0 K_S0 K+ pi- pi+ pi- PHSP;
  0.5000 pi0 K_S0 K- pi+ pi- pi+ PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

exMC_bkg_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip2pi0KsKpipipi"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_bkg_pi0
  config.cross_section   = :default
end

### ---------------------------- Event selection (BOSS) ---------------------------- ###
alg_name = "Psip2gammaetac2S"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.686]})

sel = Selection.new
sel.select_track {                    # exactly six charged tracks, net charge zero
      cos_theta 0.93                  # MDC angular coverage |cos(theta)| < 0.93
      Vz        10.0                  # |Vz| < 10 cm along the beam axis
      Vr        1.0                   # Vr < 1 cm transverse to the beam axis
      nChrp     "==3"
      nChrn     "==3"
      nNet      "==0"
    }
   .select_photon {                   # at least one good photon (the M1 transition photon)
      tdc_emc_start     0
      tdc_emc_end       14            # EMC timing window against noise
      energyThreshold_b 0.025         # barrel  E > 25 MeV (|cos(theta)| < 0.8)
      energyThreshold_e 0.025         # endcap  E > 25 MeV (0.86 < |cos(theta)| < 0.92)
      nGam              ">=1"
    }
   # K_S0 -> pi+ pi- from secondary vertex fits over all oppositely charged track
   # pairs (tracks assumed to be pions). The combination with the best vertex fit
   # quality is kept; the K_S0 invariant mass must be within 10 MeV/c^2 of the
   # nominal mass and the secondary vertex well separated from the IP.
   .assign({:chrgp => :pip, :chrgn => :pim})
   .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
   # The remaining four charged tracks are one kaon and three pions. Combinatorial
   # chi2-sum PID (dE/dx + TOF) assigns the particle types by minimising the total
   # chi-square; here only the hadron part is applied, the chi2_4C part and the
   # chi2_total < 60 requirement are applied in ROOT.
   .pid(method: :chi2_sum) {
      chi_min_cut 4
      identify :kaon, :pion
    }
   # Nominal 4C kinematic fit under psi(3686) -> gamma pi+ pi- K+- pi-+ pi+ pi-
   # (participant list shown for the K+ pi- pi+ pi- arrangement; the three
   # remaining particle-combination hypotheses are resolved by the chi2-sum PID
   # assignment together with the fit's automatic combination iteration). The
   # paper requires chi2_total (chi2_4C + chi2_K + chi2_pi) < 60, and uses a 3C
   # fit in which the photon energy is left free; the tight requirement is
   # optimised and applied in ROOT.
   .kinematic_fit([:gamma, :K_S0, :kp, :pip, :pim, :pim]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
    }

alg.note(:pid_method, "PID combines dE/dx and TOF into
    chi2_PID(i) = ((dE/dx_meas - dE/dx_exp)/sigma_dEdx)^2 +
    ((TOF_meas - TOF_exp)/sigma_TOF)^2 for the pion, kaon and proton hypotheses.
    Four particle-combination hypotheses (K+ pi- pi+ pi-, pi+ K- pi+ pi-,
    pi+ pi- K+ pi-, pi+ pi- pi+ K-) are each subjected to the 4C fit, and the
    species assignment is the one giving the smallest total chi-square
    chi2_total = chi2_4C + chi2_K + chi2_pi; the event is kept if
    chi2_total < 60. Applied in ROOT after the nominal fit.")
  .note(:signal_combination, "for each event the M1 photon is the good photon
    giving the minimum chi2_4C when looping over all good photons; this
    photon-selection step has no dedicated DSL primitive and is applied in ROOT")
  .note(:kinematic_fit_variant, "because a fake photon carries no useful
    information, the photon energy is set free in the kinematic fit (3C fit) to
    avoid the mass distortion caused by the 25 MeV photon energy threshold; the
    4C fit result is used for comparison. The 3C fit has no dedicated DSL
    primitive and is applied in ROOT")
  .note(:mass_window, "K_S0 candidates must satisfy |M(pi+ pi-) - M_K_S0| <
    10 MeV/c^2; the final K_S0 K 3pi mass spectrum is fitted over
    [3.30, 3.70] GeV/c^2 with the eta_c(2S), chi_cJ (J = 0, 1, 2) signals and the
    psi(3686) -> K_S0 K+- pi-+ pi+ pi-, psi(3686) -> pi0 K_S0 K+- pi-+ pi+ pi-,
    ISR and phase-space backgrounds. Applied in ROOT")
  .note(:background_veto, "three ROOT-level background vetoes: (1) reject events
    for which the recoil mass of any pi+ pi- pair is within 15 MeV/c^2 of the
    J/psi nominal mass (psi(3686) -> pi+ pi- J/psi, J/psi ->
    gamma K_S0 K+- pi-+); (2) reject events with M(K_S0 K+- pi-+) > 3.05 GeV/c^2
    (psi(3686) -> eta J/psi, eta -> gamma pi+ pi-); (3) reject events for which
    the mass of any gamma pi+ pi- combination is within 20 MeV/c^2 of the
    nominal eta' mass (psi(3686) -> eta' K_S0 K+- pi-+, eta' -> gamma pi+ pi-)")
  .note(:background_sources, "from inclusive MC the dominant background is
    psi(3686) -> K_S0 K+- pi-+ pi+ pi- in which a fake photon or an FSR photon
    enters the final state; further backgrounds are psi(3686) ->
    pi0 K_S0 K+- pi-+ pi+ pi- with a missing photon and ISR. The FSR component
    is scaled by the ratio of FSR fractions in data and MC (scale factor 1.46),
    determined from the control samples psi(3686) -> gamma pi+ pi- K+ K- and
    psi(3686) -> gamma pi+ pi- pi+ pi-")
  .note(:continuum, "the continuum (including ISR) background is estimated from
    the 3.65 GeV data scaled by f_continuum = 3.6 for luminosity and cross
    section differences, with particle momenta and energies rescaled to the
    beam-energy difference; its K_S0 K 3pi mass shape is used in the final fit")
  .note(:efficiency_curve, "the detection efficiency for the signal selection is
    11.1%, determined from MC in which the psi(3686) -> gamma eta_c(2S) M1
    transition and the eta_c(2S) -> K_S0 K+- pi-+ pi+ pi- decay are generated
    with the proper angular distributions")

alg.with_decay_card(decay_card_signal).apply(sel)

### --------------------------------- Execution --------------------------------- ###
alg.execute_on([psip_data, psip_incMC, cont_data,
                exMC_signal, exMC_bkg_phsp, exMC_bkg_pi0])
