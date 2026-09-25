# ================================================================
# Datasets at sqrt(s) = 3.773 GeV
# ================================================================
data_3773  = DatasetManager.real_data.find("712_3773")     # psi(3770) real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

# ----------------------------------------------------------------
# Decay cards (EvtGen format)
# ----------------------------------------------------------------
# Channel 1: psi(3770) -> D0 D0bar,
#   signal D0 -> K- pi+ pi- e+ nu_e,
#   tag    D0bar -> K+ pi- / K+ pi- pi0 / K+ pi- pi+ pi-
decay_card_d0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ pi- e+ nu_e PHSP;
    Enddecay

    Decay anti-D0
    0.3333 K+ pi-          PHSP;
    0.3333 K+ pi- pi0      PHSP;
    0.3334 K+ pi- pi+ pi-  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Channel 2: psi(3770) -> D+ D-,
#   signal D+ -> K- pi+ pi0 e+ nu_e,
#   tag    D- -> K+ pi- pi- / K_S0 pi- / K+ pi- pi- pi0 /
#                K_S0 pi- pi0 / K_S0 pi- pi+ pi- / K+ K- pi-   (K_S0 -> pi+ pi-)
decay_card_dplus = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 K- pi+ pi0 e+ nu_e PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay D-
    0.1667 K+ pi- pi-         PHSP;
    0.1667 K_S0 pi-           PHSP;
    0.1666 K+ pi- pi- pi0     PHSP;
    0.1666 K_S0 pi- pi0       PHSP;
    0.1667 K_S0 pi- pi+ pi-   PHSP;
    0.1667 K+ K- pi-          PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# ----------------------------------------------------------------
# Exclusive MC: 1,000,000 events per channel
# ----------------------------------------------------------------
exMC_d0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_d0_kpipi_e_nu"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_d0
  config.cross_section   = :default
end

exMC_dplus = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_dplus_kpipi0_e_nu"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_dplus
  config.cross_section   = :default
end

# ================================================================
# Channel 1 :  D0 -> K- pi+ pi- e+ nu_e   (tag side = D0bar)
# ================================================================
alg_d0 = TagAnalysis.new("D0SemilepTag")
alg_d0.set_header(["D0SemilepTagAlg/D0SemilepTag.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_d0)

# Hadronic tag side built from pre-stored DTag candidates (single tag)
alg_d0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi   # D0bar tags: K+pi-, K+pi-pi0, K+pi-pi+pi-
  t.charm -1                                    # pin the tagged side to D0bar
end

# Semileptonic signal side: one K-, one pi+, one pi-, one e+, missing nu_e,
# net charge 0, and NO additional charged tracks (exact multiset).
alg_d0.signal_side do |s|
  s.charged(km: 1, pip: 1, pim: 1, ep: 1)
  s.missing :nu_e
  s.require_charge 0
end

# Four-momentum constrained kinematic fit
alg_d0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# BOSS-side procedures with no DSL construct
alg_d0.note(:pid_correction_method,
  "signal-side electron identification: CL_e > 0.001 and CL_e/(CL_e+CL_pi+CL_K) > 0.8, " \
  "with hadron suppression E/p - 0.05*chi2(e-dE/dx) > 0.60 (D0 channel). These criteria are " \
  "not tunable through the tag DSL; the signal-side e+ uses the built-in electron PID.")
alg_d0.note(:efficiency_curve,
  "final-state-radiation photons within 5 degrees of the signal e+ are merged into the " \
  "electron four-momentum before the kinematic fit.")

alg_d0.apply
alg_d0.execute_on([data_3773, incMC_3773, exMC_d0])

# ================================================================
# Channel 2 :  D+ -> K- pi+ pi0 e+ nu_e   (tag side = D-)
# ================================================================
alg_dp = TagAnalysis.new("DplusSemilepTag")
alg_dp.set_header(["DplusSemilepTagAlg/DplusSemilepTag.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_dplus)

# Hadronic tag side from pre-stored DTag candidates (single tag)
alg_dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1   # pin the tagged side to D-
end

# Semileptonic signal side: one K-, one pi+, one e+, a pi0 from two photons,
# missing nu_e, net charge +1, and NO additional charged tracks.
alg_dp.signal_side do |s|
  s.charged(km: 1, pip: 1, ep: 1)
  s.photons 2                 # two photons feed the pi0 -> gamma gamma
  s.missing :nu_e
  s.require_charge 1
end

# Four-momentum constrained fit; the photon pair is mass-constrained to the pi0
alg_dp.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

# BOSS-side procedures with no DSL construct
alg_dp.note(:pid_correction_method,
  "signal-side electron identification: CL_e > 0.001 and CL_e/(CL_e+CL_pi+CL_K) > 0.8, " \
  "with hadron suppression E/p - 0.05*chi2(e-dE/dx) > 0.53 (D+ channel). Not tunable through " \
  "the tag DSL; the signal-side e+ uses the built-in electron PID.")
alg_dp.note(:efficiency_curve,
  "final-state-radiation photons within 5 degrees of the signal e+ are merged into the " \
  "electron four-momentum before the kinematic fit.")
alg_dp.note(:D_momentum_constraint,
  "during the kinematic fit the signal D+ momentum is constrained to " \
  "p_D = -p_hat_Dbar * sqrt(E_beam^2 - m_Dbar^2), i.e. back-to-back with the tag D- with " \
  "magnitude fixed by the beam energy; this constraint has no dedicated DSL construct.")

alg_dp.apply
alg_dp.execute_on([data_3773, incMC_3773, exMC_dplus])