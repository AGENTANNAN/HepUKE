# Observation of chi_cJ decays to phi phi
#   [arXiv:2301.12922]
#
# psi(3686) -> gamma chi_cJ, chi_cJ -> phi phi, phi -> K+ K-
# Single energy: sqrt(s) = 3.686 GeV, (447.9 +/- 2.3) x 10^6 psi(2S) events
# Ordinary analysis: single Algorithm covering all chi_cJ states (shared final state)

### Dataset preparation ###
psi2s_data  = DatasetManager.real_data.find("709_3686")
psi2s_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: psi(3686) -> gamma chi_cJ, chi_cJ -> phi phi, phi -> K+ K-
decay_card = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_cJ   VSP_PWAVE;
  Enddecay

  Decay chi_cJ
  1.0000 phi phi   PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K-   VSS;
  Enddecay

  End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psip_gamma_chicJ_phiphi"
  c.related_dataset = psi2s_data
  c.events          = 1_000_000
  c.decay_card      = decay_card
  c.cross_section   = :default
end

### Event selection (BOSS) — Ordinary analysis ###
alg = Algorithm.new("ChicJToPhiPhi")
alg.set_header(["ChicJToPhiPhiAlg/ChicJToPhiPhi.h"])
   .set_constant({"ECMS" => [:double, 3.686]})

sel = Selection.new

# Charged tracks: exactly four kaon tracks, net charge zero
sel.select_track {
       cos_theta 0.93
       Vz        10.0
       Vr        1.0
       nChrp     ">=2"
       nChrn     ">=2"
       nTot      4
       nNet      0
     }
     .pid(method: :probability) {
       prob_cut 0.001
       identify :kp, :km, against: [:pion]
     }
     # At least one photon, E > 25 MeV (barrel) / 50 MeV (endcap)
     .select_photon {
       energyThreshold_b 0.025
       energyThreshold_e 0.050
       angle_to_track    10.0
       tdc_emc_start     0
       tdc_emc_end       700
       nGam              ">=1"
     }

# 4C kinematic fit: psi(3686) -> gamma 2(K+ K-), chi2_4C < 60
sel.kinematic_fit([:gamma, :kp, :km, :kp, :km]) do
  constrain_four_momentum
  chi2_cut 60
  nominal
end

# The paper requires an ndof check (ndof = 4) which is implicit in 4C fit.
# After the 4C fit, if multiple photon candidates survive, the one with the
# smallest chi2_4C is kept. The four kaons are paired into two phi candidates
# using the combination that minimizes:
#   delta = sqrt((M(1)_{K+K-} - m_phi)^2 + (M(2)_{K+K-} - m_phi)^2)
# These steps are applied at the ROOT level.

alg.with_decay_card(decay_card).apply(sel)

alg.note(:phi_pair_selection,
        "The four charged kaon tracks are paired into two phi candidates. " \
        "The combination that minimizes delta = sqrt((M(1)_{K+K-} - m_phi)^2 + " \
        "(M(2)_{K+K-} - m_phi)^2) is selected. If multiple radiative photons survive " \
        "the 4C fit, the one giving the smallest chi2_4C is retained.")
  .note(:photon_selection,
        "Photon selection: E > 25 MeV in the barrel (|cos(theta)| < 0.80), " \
        "E > 50 MeV in the endcap (0.86 < |cos(theta)| < 0.92). " \
        "EMC timing: 0 <= T <= 700 ns. Opening angle between the photon and " \
        "any charged track > 10 degrees.")
  .note(:kaon_pid,
        "Kaon PID uses the probability method. Each track is identified as a kaon " \
        "if its kaon probability exceeds its pion probability. " \
        "The tracking efficiency is corrected using data control samples.")
  .note(:phi_signal_region,
        "Phi signal region: |M(K+K-) - m_phi| < 0.015 GeV/c^2. " \
        "Sideband regions (0.045, 0.060) and (0.030, 0.060) GeV/c^2 from the " \
        "nominal phi mass are used for background estimation (ROOT level).")
  .note(:background,
        "Dominant backgrounds: continuum e+e- -> phi K+ K-, " \
        "psi(3686) -> phi K+ K- (non-phi), " \
        "psi(3686) -> omega phi (omega -> pi+ pi- pi0), " \
        "psi(3686) -> eta phi (eta -> gamma gamma). " \
        "Estimated using inclusive MC and sideband methods (ROOT level).")
  .note(:chicJ_extraction,
        "Signal yields for chi_c0, chi_c1, chi_c2 are extracted from simultaneous " \
        "unbinned ML fits to M(K+K-K+K-) distributions with the photon energy " \
        "constrained to the radiative transition photon region " \
        "E_gamma in [0.14, 0.22] GeV for chi_c0 and [0.08, 0.22] GeV for chi_c1,2. " \
        "The chi_cJ mass region [3.40, 3.56] GeV/c^2 is used (ROOT level).")
  .note(:efficiency,
        "Detection efficiencies are determined from signal MC samples generated " \
        "with a phase-space model for chi_cJ -> phi phi. Systematic uncertainties " \
        "include tracking (1.0% per track), PID (1.0% per track), photon detection " \
        "(1.0%), kinematic fit (1.0%), and signal model dependence.")

alg.execute_on([psi2s_data, psi2s_incMC, exMC])