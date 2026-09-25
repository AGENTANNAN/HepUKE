# ============================================================================
# 1605.00068v2  Measurement of the absolute branching fraction of
#               D+ -> K0bar mu+ nu_mu  (double-tag technique at sqrt(s) = 3.773 GeV)
#
# Data: 2.93 fb^-1 collected at the psi(3770) peak.
# Tag side: the D- single tag (ST) reconstructed in six hadronic modes:
#   K+ pi- pi-, K_S0 pi-, K+ pi- pi- pi0, K_S0 pi- pi0,
#   K_S0 pi+ pi- pi-, K+ K- pi-.
# Signal side: D+ -> K0bar mu+ nu_mu (DT), with K0bar -> pi+ pi- or K0 -> pi0 pi0.
#
# Tag-based analysis: TagAnalysis (DTagAlg / DTagTool), no Selection object.
# BOSS part only: decay cards, exclusive MC, tag declarations and the
# kinematic fit. U_miss / M(K0bar mu+) / E_max^extra gamma windows are ROOT-level.
# ============================================================================

### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")   # psi(3770) data, 2.93 fb^-1
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # generic psi(3770) inclusive MC

# ----------------------------------------------------------------------------
# Signal decay cards: psi(3770) -> D+ D-, one D- decaying to a hadronic tag
# mode, the other (D+) to the semileptonic signal. The signal is simulated with
# the modified pole model for the q^2 dependence.
# ----------------------------------------------------------------------------
# D+ -> K_S0 mu+ nu_mu (K0bar -> K_S0 -> pi+ pi-)
decay_card_KSpimu = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 K_S0 mu+ nu_mu PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# D+ -> K0(pi0 pi0) mu+ nu_mu
decay_card_K00mupi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 K0 mu+ nu_mu PHSP;
    Enddecay

    Decay K0
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
# Exclusive signal MC samples
# ----------------------------------------------------------------------------
exMC_KSpimu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3770_Dplus_to_KS_mu_nu"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_KSpimu
  config.cross_section   = :default
end

exMC_K00mupi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi3770_Dplus_to_K0pi0pi0_mu_nu"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_K00mupi
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
# Algorithm I : D+ -> K0bar(-> pi+ pi-) mu+ nu_mu
########################################################################
alg_name_KS = "DplusToKSmuNuDTag"
alg_KS = TagAnalysis.new(alg_name_KS)
alg_KS.set_header(["#{alg_name_KS}Alg/#{alg_name_KS}.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_KSpimu)

# Tag side: the hadronically decaying D- (single tag), charge pinned.
alg_KS.tag_side(:Dplus) do |t|
  t.modes(*dtag_modes)
  t.charm -1
end

# Signal side: K0bar -> pi+ pi- (the two pions) plus the muon; the neutrino is
# undetectable and enters as one massless missing particle.
alg_KS.signal_side do |s|
  s.charged(pip: 1, pim: 1, mup: 1)
  s.require_charge 1            # K_S0 -> pi+ pi- is neutral; the D+ charge is carried
                                # by the mu+ alone
  s.min_photon_angle 10.0
  s.missing :nu_mu              # MASSLESS (resolved mass 0.0)
end

alg_KS.fit do |f|
  f.constrain_four_momentum     # 4C: tag + pi+ + pi- + mu+ + nu = ecms_lab
  f.chi2_cut 200                # loose BOSS-pass cut; the tight selection is applied in ROOT
end

alg_KS
  .note(:tag_side_selection,
        "Charged tracks: |cos(theta)| < 0.93; tracks not from K_S0 must satisfy " \
        "V_xy < 1.0 cm and |V_z| < 10.0 cm. Kaon/pion separation uses the combined " \
        "dE/dx and TOF confidence levels (CL_K > CL_pi for a kaon, CL_pi > CL_K for a pion). " \
        "K_S0 daughters only require |V_z| < 20 cm, are assigned pi+ pi- without PID, " \
        "must form a common vertex with |M(pi+pi-) - M_K_S0| < 12 MeV/c^2 and a decay " \
        "length larger than 2 standard deviations of the vertex resolution away from the IP. " \
        "Photons: shower time within 700 ns of the event start, E > 25 (50) MeV in the " \
        "EMC barrel (endcap), opening angle to any charged track > 10 degrees. " \
        "pi0: (0.115, 0.150) GeV/c^2 in M(gamma gamma) plus a mass-constrained fit on the " \
        "gamma gamma pair.")
  .note(:tag_deltaE_windows,
        "For each ST mode only the combination with the minimum |DeltaE| is kept; " \
        "DeltaE(Gamma) is required within (-25, +25) MeV for the K+ pi- pi-, K_S0 pi-, " \
        "K_S0 pi+ pi- pi- and K+ K- pi- final states, and within (-55, +40) MeV for " \
        "K+ pi- pi- pi0 and K_S0 pi- pi0.")
  .note(:tag_mBC_signal_region,
        "The ST yields are obtained from a fit to the M_BC distributions (MC signal shape " \
        "convoluted with a double Gaussian plus an ARGUS combinatorial background). " \
        "Candidates in 1.863 < M_BC < 1.877 GeV/c^2 form the ST signal region used for the " \
        "DT analysis; the total ST yield is N_ST^tot = 1522474 +- 2215.")
  .note(:tag_side_selection,
        "The same charged-track, photon, pi0 and K_S0 criteria are used on both the ST and " \
        "the DT sides.")
  .note(:signal_side_selection,
        "Exactly one additional good charged track, with charge opposite to the tag D-, is " \
        "required. Muon identification combines dE/dx, TOF and EMC confidence levels: " \
        "CL_mu > CL_K, CL_mu > CL_e and CL_mu > 0.001; the EMC energy deposited by the " \
        "muon must lie within (0.1, 0.3) GeV to reject mis-identified pions.")
  .note(:background_veto,
        "Backgrounds dominated by D+ -> K0bar pi+ (pi0) are suppressed with " \
        "M(K0bar mu+) < 1.6 GeV/c^2 and E_max^extra gamma < 0.15 GeV, where the latter is " \
        "the maximum energy of any photon not used in the DT event selection; both windows " \
        "are applied at the ROOT stage on the stored observables.")
  .note(:missing_mass_fit,
        "The neutrino is undetectable, so U_miss = E_miss - |p_miss| is constructed with " \
        "E_miss = E_beam - E_K0bar - E_mu+ and p_miss = |p_D+ - p_K0bar - p_mu+|, where " \
        "p_D+ = (-phat_D-ST) * sqrt(E_beam^2 - m_D+^2) uses the direction of the ST D- and " \
        "the nominal D+ mass. The signal yield comes from a simultaneous fit to the two " \
        "U_miss distributions (K0bar -> pi+ pi- and K0 -> pi0 pi0); this fit is part of the " \
        "ROOT stage.")

alg_KS.apply

########################################################################
# Algorithm II : D+ -> K0(-> pi0 pi0) mu+ nu_mu
########################################################################
alg_name_00 = "DplusToK0pi0pi0muNuDTag"
alg_00 = TagAnalysis.new(alg_name_00)
alg_00.set_header(["#{alg_name_00}Alg/#{alg_name_00}.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_K00mupi)

alg_00.tag_side(:Dplus) do |t|
  t.modes(*dtag_modes)
  t.charm -1
end

# Signal side: K0 -> pi0 pi0 (four photons) plus the muon and the missing neutrino
alg_00.signal_side do |s|
  s.photons 4                   # two pi0 -> four photons from the showers used by neither tag
  s.charged(mup: 1)
  s.require_charge 1
  s.min_photon_angle 10.0
  s.missing :nu_mu              # MASSLESS (resolved mass 0.0)
end

alg_00.fit do |f|
  f.constrain_four_momentum     # 4C: tag + gamma gamma gamma gamma + mu+ + nu = ecms_lab
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # first pi0
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # second pi0
  f.chi2_cut 200
end

alg_00
  .note(:tag_side_selection,
        "Identical six hadronic D- single-tag modes and ST criteria as Algorithm I " \
        "(charged tracks |cos(theta)| < 0.93, V_xy < 1.0 cm, |V_z| < 10.0 cm for non-K_S0 " \
        "tracks; CL_K > CL_pi / CL_pi > CL_K kaon-pion separation; K_S0: |V_z| < 20 cm, " \
        "pi+ pi- common vertex with |M(pi+pi-) - M_K_S0| < 12 MeV/c^2 and decay length " \
        "> 2 sigma; photons E > 25 (50) MeV barrel (endcap), shower time within 700 ns, " \
        "angle to charged tracks > 10 degrees).")
  .note(:tag_deltaE_windows,
        "As in Algorithm I: minimum |DeltaE| combination per ST mode, with " \
        "|DeltaE| windows (-25, +25) MeV (K+ pi- pi-, K_S0 pi-, K_S0 pi+ pi- pi-, K+ K- pi-) " \
        "and (-55, +40) MeV (K+ pi- pi- pi0, K_S0 pi- pi0).")
  .note(:tag_mBC_signal_region,
        "ST signal region 1.863 < M_BC < 1.877 GeV/c^2; N_ST^tot = 1522474 +- 2215.")
  .note(:k0_pi0pi0_selection,
        "The K0(pi0 pi0) candidate is required to have M(pi0 pi0) within " \
        "(0.45, 0.51) GeV/c^2; if more than one combination survives, the one with the " \
        "minimum chi2_1(pi0 -> gamma gamma) + chi2_2(pi0 -> gamma gamma) is kept, where " \
        "chi2_1 and chi2_2 are the chi-squares of the mass-constrained fits on the two " \
        "pi0 -> gamma gamma candidates.")
  .note(:signal_side_selection,
        "Exactly one additional good charged track with charge opposite to the tag D-, " \
        "identified as a muon: CL_mu > CL_K, CL_mu > CL_e, CL_mu > 0.001, and EMC energy " \
        "deposition within (0.1, 0.3) GeV.")
  .note(:background_veto,
        "M(K0bar mu+) < 1.6 GeV/c^2 and E_max^extra gamma < 0.15 GeV, applied in ROOT on " \
        "the stored observables.")
  .note(:missing_mass_fit,
        "Simultaneous fit of the two U_miss distributions (K0bar -> pi+ pi- and " \
        "K0 -> pi0 pi0) in the ROOT stage; the observed DT yields are " \
        "N_DT^{+-,obs} = 16516 +- 130 and N_DT^{00,obs} = 4198 +- 33, corresponding to " \
        "N_DT^prd = 132712 +- 1041 after efficiency and daughter-branching-fraction " \
        "correction.")

alg_00.apply

# Execute on real data, inclusive MC and the two signal exclusive MC samples
root_files_KS = alg_KS.execute_on(datasets + [exMC_KSpimu])
root_files_00 = alg_00.execute_on(datasets + [exMC_K00mupi])