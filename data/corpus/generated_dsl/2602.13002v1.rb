# ============================================================================
# ψ(3770) → D0 D̄0 quantum-correlated measurement
# Coherence factors + strong-phase differences for
#   D0 → K−π+π+π−  (K3π)   and   D0 → K−π+π0
# ============================================================================

### Datasets ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # 7.93 fb−1 real data at √s = 3.773 GeV
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # matching inclusive MC

### Decay cards (signal MC: ψ(3770) → D0 D̄0, signal D0 by phase space) ###
decay_card_k3pi = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- pi- pi+ PHSP;
    Enddecay

    End
DECAYCARD

decay_card_kpipi0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC: 1.5M events for each signal mode ###
exMC_k3pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3770_sigK3pi"
  config.related_dataset = psi3770_data
  config.events          = 1_500_000
  config.decay_card      = decay_card_k3pi
  config.cross_section   = :default
end

exMC_kpipi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3770_sigKPiPi0"
  config.related_dataset = psi3770_data
  config.events          = 1_500_000
  config.decay_card      = decay_card_kpipi0
  config.cross_section   = :default
end

# ============================================================================
# Double-tag reconstruction — signal mode D0 → K−π+π+π−
# ============================================================================
alg_name_k3pi = "D0DTagK3pi"
alg_k3pi = TagAnalysis.new(alg_name_k3pi)
alg_k3pi.set_header(["#{alg_name_k3pi}Alg/#{alg_name_k3pi}.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .with_decay_card(decay_card_k3pi)

# Side 1 — signal side: D0 → K−π+π+π−
alg_k3pi.tag_side(:D0) do |t|
  t.modes :D0toKPiPiPi
  t.charm -1
  t.rank_by :inv   # keep the DT pair whose average invariant mass is closest to m_D0
end

# Side 2 — tag side: flavour / CP-even / CP-odd / self-conjugate D0 tags
alg_k3pi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi,        # flavour tags
          :D0toKK, :D0toPiPi, :D0toPiPiPi0,           # CP-even tags
          :D0toKsPi0Pi0, :D0toKlPi0, :D0toKlOmega,
          :D0toKsPi0, :D0toKsEta, :D0toKsOmega,       # CP-odd tags
          :D0toKsEtaP, :D0toKsPhi, :D0toKlPi0Pi0,
          :D0toKsPiPi, :D0toKlPiPi                    # self-conjugate tags
  t.charm 1
  t.rank_by :inv
end

# 6C kinematic fit: total four-momentum + both D0 candidate masses to nominal m_D0
alg_k3pi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# --- BOSS-side procedures with no formal DSL construct ---
alg_k3pi
  .note(:kl0_missing_mass_reconstruction,
        "Tag modes containing a K_L0 (K_L0 pi0, K_L0 omega, K_L0 pi0 pi0, K_L0 pi+ pi-) " \
        "are reconstructed with a missing-mass technique against the recoil D; events with " \
        "extra pi0 or charged tracks, any eta -> gamma gamma candidate, or pi0-pi0 shower " \
        "sharing are vetoed.")
  .note(:background_veto,
        "A K_S0 flight-distance veto is applied to D -> pi+ pi- pi0 and D -> K- pi+ pi+ pi- " \
        "to suppress K_S0 backgrounds.")
  .note(:pid_correction_method,
        "Like-sign tags use tightened PID: a track is called K if P_K > 100 P_pi and called " \
        "pi if P_pi > 100 P_K.")
  .note(:k3pi_mass_constraint_binned,
        "For the binned K3pi (K- pi+ pi+ pi-) analysis an additional mass constraint to " \
        "nominal m_D0 is applied on the K- pi+ pi+ pi- side.")

alg_k3pi.apply
alg_k3pi.execute_on([psi3770_data, psi3770_incMC, exMC_k3pi])

# ============================================================================
# Double-tag reconstruction — signal mode D0 → K−π+π0
# ============================================================================
alg_name_kpipi0 = "D0DTagKPiPi0"
alg_kpipi0 = TagAnalysis.new(alg_name_kpipi0)
alg_kpipi0.set_header(["#{alg_name_kpipi0}Alg/#{alg_name_kpipi0}.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .with_decay_card(decay_card_kpipi0)

# Side 1 — signal side: D0 → K−π+π0
alg_kpipi0.tag_side(:D0) do |t|
  t.modes :D0toKPiPi0
  t.charm -1
  t.rank_by :inv   # keep the DT pair whose average invariant mass is closest to m_D0
end

# Side 2 — tag side: same tag set as K3pi, but with the like-sign K−π+π+π− tag dropped
alg_kpipi0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0,                      # flavour tags (K3pi tag removed)
          :D0toKK, :D0toPiPi, :D0toPiPiPi0,           # CP-even tags
          :D0toKsPi0Pi0, :D0toKlPi0, :D0toKlOmega,
          :D0toKsPi0, :D0toKsEta, :D0toKsOmega,       # CP-odd tags
          :D0toKsEtaP, :D0toKsPhi, :D0toKlPi0Pi0,
          :D0toKsPiPi, :D0toKlPiPi                    # self-conjugate tags
  t.charm 1
  t.rank_by :inv
end

# 6C kinematic fit: total four-momentum + both D0 candidate masses to nominal m_D0
alg_kpipi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# --- BOSS-side procedures with no formal DSL construct ---
alg_kpipi0
  .note(:kl0_missing_mass_reconstruction,
        "Tag modes containing a K_L0 (K_L0 pi0, K_L0 omega, K_L0 pi0 pi0, K_L0 pi+ pi-) " \
        "are reconstructed with a missing-mass technique against the recoil D; events with " \
        "extra pi0 or charged tracks, any eta -> gamma gamma candidate, or pi0-pi0 shower " \
        "sharing are vetoed.")
  .note(:background_veto,
        "The like-sign K- pi+ pi+ pi- tag is dropped in the K- pi+ pi0 signal analysis. " \
        "A K_S0 flight-distance veto is applied to D -> pi+ pi- pi0 and D -> K- pi+ pi+ pi-.")
  .note(:pid_correction_method,
        "Like-sign tags use tightened PID: a track is called K if P_K > 100 P_pi and called " \
        "pi if P_pi > 100 P_K.")

alg_kpipi0.apply
alg_kpipi0.execute_on([psi3770_data, psi3770_incMC, exMC_kpipi0])