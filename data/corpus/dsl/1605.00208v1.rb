# ============================================================================
# 1605.00208v1  Measurement of the absolute branching fraction of
#               D+ -> K0 e+ nu_e  via K0 -> K_S0 -> pi0 pi0
#               (double-tag technique at sqrt(s) = 3.773 GeV)
#
# Data: 2.93 fb^-1 collected at the psi(3770) peak.
# Tag side: the D- single tag (ST) reconstructed in six hadronic modes:
#   K+ pi- pi-, K_S0 pi-, K+ pi- pi- pi0, K_S0 pi- pi0,
#   K_S0 pi+ pi- pi-, K+ K- pi-.
# Signal side: D+ -> K0 e+ nu_e (DT) with K0 -> K_S0 -> pi0 pi0 (four photons).
#
# Tag-based analysis: TagAnalysis (DTagAlg / DTagTool), no Selection object.
# BOSS part only: decay cards, exclusive MC, tag declarations and the
# kinematic fit. U_miss / M(K0 pi0 pi0) / E_max^extra gamma windows are ROOT-level.
# ============================================================================

### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # psi(3770) data, 2.93 fb^-1
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # generic psi(3770) inclusive MC

# ----------------------------------------------------------------------------
# Signal decay card: psi(3770) -> D+ D-, the D- decaying to a hadronic tag
# mode and the D+ to the semileptonic signal with K0 -> K_S0 -> pi0 pi0.
# The signal is simulated with the modified pole model for the q^2 dependence.
# ----------------------------------------------------------------------------
decay_card_K0enu = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 K0 e+ nu_e PHSP;
    Enddecay

    Decay K0
    1.0000 K_S0 pi0 pi0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi0 pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# ----------------------------------------------------------------------------
# Exclusive signal MC sample
# ----------------------------------------------------------------------------
exMC_K0enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3770_Dplus_to_K0_e_nu"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_K0enu
  config.cross_section   = :default
end

datasets = [psi3770_data, psi3770_incMC]

# Six hadronic D- single-tag modes (the tag is the D-):
#   K+ pi- pi-, K_S0 pi-, K+ pi- pi- pi0, K_S0 pi- pi0,
#   K_S0 pi+ pi- pi-, K+ K- pi-.
# K_S0 -> pi+ pi- and pi0 -> gamma gamma are handled inside the tag.
dtag_modes = [
  :DptoKPiPi,
  :DptoKsPi,
  :DptoKPiPiPi0,
  :DptoKsPiPi0,
  :DptoKsPiPiPi,
  :DptoKKPi
]

########################################################################
# D+ -> K0(-> K_S0 -> pi0 pi0) e+ nu_e
########################################################################
alg_name = "DplusToK0eNuDTag"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(decay_card_K0enu)

# Tag side: the hadronically decaying D- (single tag), charge pinned.
alg.tag_side(:Dplus) do |t|
  t.modes(*dtag_modes)
  t.charm -1
end

# Signal side: K0 -> K_S0 -> pi0 pi0 (four photons) plus the electron; the
# neutrino is undetectable and enters as one massless missing particle.
alg.signal_side do |s|
  s.photons 4                   # two pi0 -> four photons from the showers used by neither tag
  s.charged(ep: 1)
  s.require_charge 1            # the D+ charge is carried by the e+ alone (K0 is neutral)
  s.min_photon_angle 10.0
  s.missing :nu_e               # MASSLESS (resolved mass 0.0)
end

alg.fit do |f|
  f.constrain_four_momentum     # 4C: tag + gamma gamma gamma gamma + e+ + nu = ecms_lab
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # first pi0
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # second pi0
  f.chi2_cut 200                # loose BOSS-pass cut; the tight selection is applied in ROOT
end

alg
  .note(:tag_side_selection,
        "Charged tracks: |cos(theta)| < 0.93; tracks not from K_S0 must satisfy " \
        "V_xy < 1.0 cm and |V_z| < 10.0 cm. Kaon/pion separation uses the combined " \
        "dE/dx and TOF confidence levels (CL_K > CL_pi for a kaon, CL_pi > CL_K for a pion). " \
        "K_S0 daughters only require |V_z| < 20 cm, are assigned pi+ pi- without PID, " \
        "must form a common vertex with |M(pi+pi-) - M_K_S0| < 12 MeV/c^2 and a decay " \
        "length L with L/sigma_L > 2. Photons: shower time within 700 ns of the event " \
        "start, E > 25 (50) MeV in the EMC barrel (endcap), opening angle to the closest " \
        "charged track > 10 degrees. pi0: (0.115, 0.150) GeV/c^2 in M(gamma gamma) plus a " \
        "mass-constrained fit on the gamma gamma pair.")
  .note(:tag_deltaE_windows,
        "For each ST mode only the combination with the minimum |DeltaE| is kept; " \
        "DeltaE is required within (-25, +25) MeV for the K+ pi- pi-, K_S0 pi-, " \
        "K_S0 pi+ pi- pi- and K+ K- pi- final states, and within (-55, +40) MeV for " \
        "K+ pi- pi- pi0 and K_S0 pi- pi0.")
  .note(:tag_mBC_signal_region,
        "The ST yields are obtained from a fit to the M_BC distributions (MC signal shape " \
        "convoluted with a double Gaussian plus an ARGUS combinatorial background). " \
        "Candidates in 1.863 < M_BC < 1.877 GeV/c^2 form the ST signal region used for the " \
        "DT analysis; the total ST yield is N_ST^tot = 1522474 +- 2215.")
  .note(:background_veto,
        "Backgrounds associated with fake photons are suppressed by requiring " \
        "E_max^extra gamma < 0.300 GeV, where E_max^extra gamma is the maximum energy of " \
        "any photon not used in the DT event selection; the window is applied at the ROOT " \
        "stage on the stored observable.")
  .note(:k0_pi0pi0_selection,
        "The K0 -> K_S0 -> pi0 pi0 candidate is required to have M(pi0 pi0) within " \
        "(0.45, 0.51) GeV/c^2; if more than one combination survives, the one with the " \
        "minimum chi2_1(pi0 -> gamma gamma) + chi2_2(pi0 -> gamma gamma) is kept, where " \
        "chi2_1 and chi2_2 are the chi-squares of the mass-constrained fits on the two " \
        "pi0 -> gamma gamma candidates.")
  .note(:signal_side_selection,
        "At least four good photons and only one good charged track, not used in the ST " \
        "selection. Electron identification combines dE/dx, TOF and EMC confidence levels: " \
        "CL_e > 0.001 and CL_e/(CL_e + CL_pi + CL_K) > 0.8, with charge opposite to the " \
        "tag D-. To partially recover FSR and bremsstrahlung, the four-momenta of " \
        "photon(s) within 5 degrees of the initial electron direction are added to the " \
        "electron four-momentum.")
  .note(:missing_mass_fit,
        "The neutrino is undetectable, so U_miss = E_miss - |p_miss| is constructed with " \
        "E_miss = E_beam - E_K0 - E_e+ and p_miss = |p_D+ - p_K0 - p_e+|, where " \
        "p_D+ = (-phat_D-ST) * sqrt(E_beam^2 - m_D+^2) uses the direction of the ST D- and " \
        "the nominal D+ mass. The DT yield N_DT = 5013 +- 78 comes from a fit to the U_miss " \
        "distribution (signal and combinatorial background shapes from MC); this fit is " \
        "part of the ROOT stage.")
  .note(:efficiency_curve,
        "The absolute branching fraction is extracted as " \
        "B = N_DT / (N_ST^tot * epsbar * B(K0 -> pi0 pi0) * B(pi0 -> gamma gamma)^2) with the " \
        "averaged reconstruction efficiency epsbar = (25.58 +- 0.11)%, which does not include " \
        "the branching fractions of K0 -> pi0 pi0 and pi0 -> gamma gamma; K_S0 is assumed to " \
        "constitute half of the neutral-kaon decays.")

alg.apply

# Execute on real data, inclusive MC and the signal exclusive MC sample
root_files = alg.execute_on(datasets + [exMC_K0enu])
