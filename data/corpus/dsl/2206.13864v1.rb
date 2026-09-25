# Paper: 2206.13864v1 — Absolute Measurements of BFs of Cabibbo-Suppressed Hadronic D^{0(+)} Decays Involving Multiple Pions
# Type: TagAnalysis (DT method at ψ(3770))
# Dataset: 712_3773 (ψ(3770), 2.93 fb⁻¹)
# 20 signal modes: 9 D^0 and 11 D^+ decays into multiple pions (+ some η)

# ============================================================
# D^0 signal modes with D0bar tag
# ============================================================
# Signal modes: π+π-π0, π+π-2π0, π+π-2η, 4π0, 3π0η,
#   2π+2π-π0, 2π+2π-η, π+π-3π0, 2π+2π-2π0

alg_D0 = TagAnalysis.new("D0ToMultiPion_DT")

alg_D0.set_header([
  "EventModel/Event.h",
  "EvtRecEvent/EvtRecTrack.h",
  "EvtRecEvent/EvtRecDTag.h"
])

alg_D0.set_constant(ECMS: 3.773)

# Tag side: D0bar with hadronic tag modes
alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKsPiPi,
          :D0toKsPiPiPi0, :D0toKsPiPiPi, :D0toKKPi
  t.charm(-1)  # tag D0bar
end

# Signal side: signal D^0 from remaining tracks/showers
# Representative: D^0 → π+π-π0 (most abundant mode)
alg_D0.signal_side do |s|
  s.photons 2                    # at least one π0 → γγ
  s.charged(pip: 1, pim: 1)      # π+π-
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# No kinematic fit on signal side — yield extraction via 2D M_BC fits in ROOT
alg_D0.note(:no_signal_kinematic_fit,
  "No BOSS-side kinematic fit on signal side. Signal D mesons are identified " \
  "using ΔE_sig and M_BC_sig. Yields extracted from 2D unbinned ML fits to " \
  "M_BC^tag vs M_BC^sig in ROOT stage.")

alg_D0.note(:signal_modes_D0,
  "9 D^0 signal modes: π+π-π0, π+π-2π0, π+π-2η, 4π0, 3π0η, " \
  "2π+2π-π0, 2π+2π-η, π+π-3π0, 2π+2π-2π0. " \
  "Each mode has different charged/neutral multiplicities handled in ROOT. " \
  "The signal_side shown is representative (π+π-π0). " \
  "K_S^0 veto: M(π+π-) not in [0.468,0.528] GeV/c², M(π0π0) not in [0.428,0.548] GeV/c².")

alg_D0.note(:tag_deltaE_windows,
  "ΔE_tag window: (-55, 40) MeV for tag modes with π0, (-25, 25) MeV otherwise. " \
  "Best |ΔE_tag| selected per tag mode.")

alg_D0.note(:signal_deltaE_windows,
  "ΔE_sig windows vary per signal mode (see paper Table 1). Applied in ROOT stage.")

alg_D0.note(:quantum_correlation,
  "Quantum correlation correction factors applied for neutral D decays. " \
  "CP-even fractions determined for each D^0 mode. Correction factors and " \
  "residual uncertainties assigned in ROOT (see paper Table 3).")

# ============================================================
# D^+ signal modes with D^- tag
# ============================================================
# Signal modes: 2π+π-, π+2π0, 2π+π-π0, π+3π0, 3π+2π-,
#   2π+π-2π0, 2π+π-π0η, π+4π0, π+3π0η, 3π+2π-π0, 2π+π-3π0

alg_Dp = TagAnalysis.new("DpToMultiPion_DT")

alg_Dp.set_header([
  "EventModel/Event.h",
  "EvtRecEvent/EvtRecTrack.h",
  "EvtRecEvent/EvtRecDTag.h"
])

alg_Dp.set_constant(ECMS: 3.773)

# Tag side: D^- with hadronic tag modes
alg_Dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKsPiPiPi, :DptoKKPi
  t.charm(-1)  # tag D^-
end

# Signal side: signal D^+ from remaining tracks/showers
# Representative: D^+ → 2π+π- (highest-yield 3-body mode via BODY3)
alg_Dp.signal_side do |s|
  s.photons 0
  s.charged(pip: 2, pim: 1)     # 2π+π-
  s.require_charge 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

alg_Dp.note(:no_signal_kinematic_fit,
  "No BOSS-side kinematic fit on signal side. Signal D mesons are identified " \
  "using ΔE_sig and M_BC_sig. Yields extracted from 2D unbinned ML fits to " \
  "M_BC^tag vs M_BC^sig in ROOT stage.")

alg_Dp.note(:signal_modes_Dp,
  "11 D^+ signal modes: 2π+π-, π+2π0, 2π+π-π0, π+3π0, 3π+2π-, " \
  "2π+π-2π0, 2π+π-π0η, π+4π0, π+3π0η, 3π+2π-π0, 2π+π-3π0. " \
  "Each mode has different charged/neutral multiplicities handled in ROOT. " \
  "The signal_side shown is representative (2π+π-).")

alg_Dp.note(:kS0_veto,
  "K_S^0 veto: π+π- and π0π0 invariant masses outside K_S^0 mass windows " \
  "(4σ resolution). Applied in ROOT stage.")

# ============================================================
# Datasets and execution
# ============================================================

psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Exclusive MC — phase space for multi-body, BODY3 for 3-body
# Note: actual signal MC uses mixed phase-space models + BODY3 generator

alg_D0.execute_on([psi3770_data, psi3770_incMC])
alg_Dp.execute_on([psi3770_data, psi3770_incMC])