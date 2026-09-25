# BESIII: Improved measurements of the coherence factors and strong-phase differences
# in D → K-π+π+π- and D → K-π+π0 with quantum-correlated DDbar decays.
# Data: 7.93 fb^-1 at √s = 3.773 GeV (ψ(3770) → D0 Dbar0).
# Analysis: exclusive double-tag reconstruction — one D0 in a signal mode
#   (K-π+π+π- or K-π+π0), the other in a tag mode (flavour like-sign /
#   opposite-sign, CP-even, CP-odd, self-conjugate KS0π+π- / KL0π+π-).

### Dataset ###
psipp_data  = DatasetManager.real_data.find("712_3773")
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")

######################################################################
# Decay cards — signal decay chains generated once per (signal, tag)
# combination. To keep the spec compact, one card per signal mode is
# written with the tag decay set to a phase-space model; per-mode signal
# MC samples for each double-tag configuration are produced by the tag
# skeleton at execute_on time.
######################################################################

decay_card_K3pi = <<~DECAYCARD
  Decay psi(3770)
  1.0000  D0  anti-D0                          PHSP;
  Enddecay

  # Signal side: D0 -> K- pi+ pi+ pi- (LHCb CF model recommended; PHSP fallback)
  Decay D0
  1.0000  K-  pi+  pi+  pi-                    PHSP;
  Enddecay

  Decay anti-D0
  1.0000  K+  pi-  pi-  pi+                    PHSP;
  Enddecay

  End
DECAYCARD

decay_card_Kpipi0 = <<~DECAYCARD
  Decay psi(3770)
  1.0000  D0  anti-D0                          PHSP;
  Enddecay

  # Signal side: D0 -> K- pi+ pi0
  Decay D0
  1.0000  K-  pi+  pi0                         PHSP;
  Enddecay

  Decay anti-D0
  1.0000  K+  pi-  pi0                         PHSP;
  Enddecay

  Decay pi0
  1.0000  gamma  gamma                         PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC ###
exMC_K3pi   = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psipp_D0_K3pi_DT"
  c.related_dataset = psipp_data
  c.events          = 1_500_000
  c.decay_card      = decay_card_K3pi
  c.cross_section   = :default
end
exMC_Kpipi0 = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psipp_D0_Kpipi0_DT"
  c.related_dataset = psipp_data
  c.events          = 1_500_000
  c.decay_card      = decay_card_Kpipi0
  c.cross_section   = :default
end

######################################################################
# Common tag-mode groups (declared once, re-used across sig algorithms)
######################################################################
FLAVOUR_TAGS = [:D0toKPi, :D0toKPiPi0, :D0toKPiPiPi]
CP_EVEN_TAGS = [:D0toKK, :D0toPiPi, :D0toPiPiPi0, :D0toKsPi0Pi0, :D0toKlPi0, :D0toKlOmega]
CP_ODD_TAGS  = [:D0toKsPi0, :D0toKsEta, :D0toKsOmega, :D0toKsEtaP, :D0toKsPhi, :D0toKlPi0Pi0]
SELF_CONJ_TAGS = [:D0toKsPiPi, :D0toKlPiPi]

######################################################################
# Algorithm 1: signal mode D0 → K- π+ π+ π-  (double-tag)
######################################################################
alg_k3pi = TagAnalysis.new("D0K3piVsTag")
alg_k3pi.set_header(["D0K3piVsTagAlg/D0K3piVsTag.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .with_decay_card(decay_card_K3pi)

# Signal-side tag: reconstruct D0 → K- π+ π+ π- via DTagAlg
alg_k3pi.tag_side(:D0) do |t|
  t.modes :D0toKPiPiPi
end

# Other side: any of the declared tag modes (flavour / CP / self-conjugate)
alg_k3pi.tag_side(:D0) do |t|
  t.modes(*(FLAVOUR_TAGS + CP_EVEN_TAGS + CP_ODD_TAGS + SELF_CONJ_TAGS))
  t.rank_by :inv
end

alg_k3pi.fit do |f|
  f.constrain_four_momentum
  # 6C tag-side D0 mass constraints improve phase-space resolution.
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_k3pi.note(:tag_mode_table,
              "Tag categories: flavour (like-/opposite-sign) K-π+π+π-, K-π+π0, K-π+; " \
              "CP-even K+K-, π+π-, π+π-π0, KS0π0π0, KL0π0, KL0ω; " \
              "CP-odd KS0π0, KS0η, KS0ω, KS0η', KS0φ, KL0π0π0; " \
              "self-conjugate KS0π+π-, KL0π+π-.")
        .note(:kl0_tag_missing_mass,
              "Tag modes with a K_L0 use a missing-mass technique: " \
              "M_miss^2 = (√s/2 - E_X)^2 - |p_S + p_X|^2 with p_S the " \
              "signal-side momentum and (E_X, p_X) the sum over the other-side " \
              "charged tracks and π0s. Extra π0 / extra charged tracks / any " \
              "η→γγ candidate / π0-π0 shower-sharing all veto the event.")
        .note(:ks0_flight_veto,
              "K_S0 flight-distance veto applied to D → π+π-π0 to suppress " \
              "D → K_S0 π0, and to D → K-π+π+π- to suppress D → K_S0 K-π+.")
        .note(:like_sign_tight_pid,
              "Like-sign tags apply a tightened PID: track labelled K (π) if " \
              "P_K > 100 P_π (P_π > 100 P_K).")
        .note(:best_candidate_criterion,
              "For events with multiple DT candidates (~10%), the pair whose " \
              "average reconstructed invariant mass is closest to nominal m_D0 " \
              "is retained.")
        .note(:signal_mass_constraint,
              "For the binned K3π analysis, a mass-constrained fit to nominal " \
              "m_D0 is applied to the D → K-π+π+π- side to sharpen the " \
              "phase-space position resolution.")

alg_k3pi.apply
alg_k3pi.execute_on([psipp_data, psipp_incMC, exMC_K3pi])

######################################################################
# Algorithm 2: signal mode D0 → K- π+ π0  (double-tag)
######################################################################
alg_kpipi0 = TagAnalysis.new("D0Kpipi0VsTag")
alg_kpipi0.set_header(["D0Kpipi0VsTagAlg/D0Kpipi0VsTag.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .with_decay_card(decay_card_Kpipi0)

# Signal-side tag: reconstruct D0 → K- π+ π0 via DTagAlg
alg_kpipi0.tag_side(:D0) do |t|
  t.modes :D0toKPiPi0
end

# Other side: same set of tag modes (K-π+π+π- excluded as it is the K3π mode)
alg_kpipi0.tag_side(:D0) do |t|
  t.modes(*([:D0toKPi, :D0toKPiPi0] + CP_EVEN_TAGS + CP_ODD_TAGS + SELF_CONJ_TAGS))
  t.rank_by :inv
end

alg_kpipi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_kpipi0.note(:tag_mode_table,
                "Same tag categories as the K3π analysis, except the like-sign " \
                "K-π+π+π- combination is dropped for K-π+π0 signal.")
          .note(:kl0_tag_missing_mass,
                "K_L0-containing tag modes use the missing-mass technique.")
          .note(:like_sign_tight_pid,
                "Like-sign tags apply tightened PID (P_K > 100 P_π, P_π > 100 P_K).")

alg_kpipi0.apply
alg_kpipi0.execute_on([psipp_data, psipp_incMC, exMC_Kpipi0])
