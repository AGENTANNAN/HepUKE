# ============================================================
# Paper: 2006.05588v4
# Inclusive chi_cJ -> anything / J/psi -> anything distributions
# Uses psi(3686) -> gamma chi_cJ and chi_cJ -> gamma J/psi to study
# charged track, EMC shower, and pi0 multiplicity distributions
# ============================================================

data_3686 = DatasetManager.real_data.find("709_3686")
incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")

# Inclusive decay: psi(3686) -> gamma chi_cJ; chi_cJ -> anything
# For MC generation of inclusive chi_cJ decay distribution study
decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0  HELAMP 1 0 0 1 0 0;
    Enddecay
    Decay chi_c0
    1.0000 anything        PHSP;
    Enddecay
    End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psip_gamma_chi_cJ"
  config.related_dataset = data_3686
  config.events          = 500000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

alg = Algorithm.new("ChiCJInclusive")
alg.set_header(["ChiCJInclusiveAlg/ChiCJInclusive.h"])
    .set_constant({ "ECMS" => [:double, 3.686] })

sel = Selection.new
sel.select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end

alg
  .note(:analysis_type, "This is NOT a conventional hadronic-decay reconstruction.
    It measures INCLUSIVE multiplicity distributions (N_ch, N_sh, N_pi0)
    for chi_cJ -> anything and J/psi -> anything via psi(3686) radiative
    transitions. The analysis extracts inclusive event yields by fitting
    shower-energy (E_sh) distributions rather than reconstructing specific
    decay chains. The DSL event selection above captures only the basic
    track/photon/PID quality criteria applied at BOSS level.")
  .note(:radiative_transition_fit, "chi_cJ and J/psi event yields extracted
    by fitting the inclusive E_sh (shower energy) distribution. For
    psi(3686)->gamma chi_cJ: five peaks at ~172, 268, 396 MeV (chi_c2,1,0)
    plus chi_c1,2 -> gamma J/psi peaks at ~575, 780 MeV. Each E_sh fit is
    performed per N_ch (or N_sh, N_pi0) bin. Signal shapes from MC
    (convolved with bifurcated Gaussians for data resolution differences).
    Background from inclusive psi(3686) MC with radiative photons removed
    plus a 2nd-order polynomial.")
  .note(:event_selection, "Event selection: >=1 charged track (except
    N_ch=0 study where relaxed), >=1 neutral shower, minimum event energy.
    Background filters: non-psi(3686) event removal, psi(3686)->pi pi J/psi
    veto. For N_ch=0: additional requirements remove high-background regions:
    |Px_neu| > 1.0 GeV/c and |Py_neu| > 1.0 GeV/c vetoed.")
  .note(:detector_eff_correction, "Efficiency-corrected distributions
    determined using MC: N = D / epsilon, where D = detected data events,
    epsilon = MC efficiency vs N_ch/N_sh/N_pi0. MC events weighted by
    wt_pi0 (correction for data-MC pi0 multiplicity difference) and
    wt_trans (E_gamma^3 energy dependence of E1 radiative transitions).")
  .note(:response_matrices, "MC truth information used to construct matrices
    connecting detected distributions to produced (predetection) distributions.
    Assuming matrices also apply to data, produced charged particle, photon,
    and pi0 multiplicity distributions are obtained. Produced means compared
    with MARK I results for e+e- -> hadrons vs sqrt(s).")
  .note(:shower_photon_matching, "MC truth photon tagging: good shower match
    defined as D_theta < 0.1 rad between EMCSH and MC photon. Matching
    efficiency is 91.2% (|cos(theta)|<0.8, 0.25<E_sh<2 GeV).
    Good-match fraction varies from 60% (lowest E) to 89% (highest E).")
  .note(:pi0_reconstruction, "pi0 candidates: gamma-gamma pair with
    0.120 < M_gg < 0.145 GeV/c^2. Fraction of valid pi0s (R) determined
    by fitting M_gg with signal shape (MC truth) convolved with bifurcated
    Gaussian + 1st-order Chebychev background. N_pi0 may include
    miscombinations (e.g. N=3 could be 3 real or 2+1 fake, etc.).")
  .note(:lundcharm_model, "Unmeasured charmonium decays (85% for chi_c0,
    57% chi_c1, 72% chi_c2) simulated with LUNDCHARM model. Parameters
    optimized using 20M J/psi events. Comparison of data vs MC
    multiplicity distributions used to validate/improve LUNDCHARM.
    Data-MC disagreement in photon distributions suggests G-parity
    conservation should be included in MC simulation.")
  .note(:systematic_uncertainties, "Sources: E_sh fit systematics (signal
    shapes, background modeling, fit range), MC generator modeling,
    efficiency determination, pi0 miscombination fraction, continuum
    subtraction. Total systematic uncertainty varies by N_ch bin from
    1-15% range depending on statistics.")
  .with_decay_card(decay_card)
  .apply(sel)

alg.execute_on([data_3686, incMC_3686, exMC])