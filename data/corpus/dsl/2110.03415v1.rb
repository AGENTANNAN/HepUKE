# Paper: 2110.03415v1
# Title: D -> K pi omega absolute branching fractions
# Dataset: psi(3770) (712_3773), 2.93 fb^{-1}
#
# Double-tag method at psi(3770):
#   Tag side: D0bar (anti-D0) via 3 hadronic modes, or D- via 6 modes
#   Signal side: D0/D+ -> final state with omega -> pi+ pi- pi0
#
# Three independent signal modes (Rule T1):
#   1. D0 -> K- pi+ omega
#   2. D0 -> K_S0 pi0 omega
#   3. D+ -> K_S0 pi+ omega

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_3773  = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================
# Algorithm 1: D0 -> K- pi+ omega, omega -> pi+ pi- pi0
# Single tag on anti-D0 (tag_side :D0, 3 modes)
# Signal side: K- pi+ pi+ pi- pi0 (4 charged, 2 photons)
# ============================================================

alg_d0_kpi_omega = TagAnalysis.new("D0toKPiOmega")
alg_d0_kpi_omega.set_header(["D0toKPiOmegaAlg/D0toKPiOmega.h"])
                 .set_constant({ "ECMS" => [:double, 3.773] })

alg_d0_kpi_omega.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toK3Pi
  t.window :deltaE, abs: 0.025       # |DeltaE| < 25 MeV for non-pi0 tag modes
end

alg_d0_kpi_omega.signal_side do |s|
  s.photons 2                          # pi0 -> gamma gamma
  s.charged(km: 1, pip: 2, pim: 1)    # K-, pi+(D0), pi+(omega), pi-(omega)
end

alg_d0_kpi_omega.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_kpi_omega.note(:omega_reconstruction, "omega -> pi+ pi- pi0. Omega mass window: |M_3pi - M_omega| < 25 MeV/c^2. " \
                                              "Two pi+pi-pi0 combinations possible (2 pi+ in final state); " \
                                              "both invariant masses denoted M_pi1+pi-pi0 and M_pi2+pi-pi0. " \
                                              "Events in 2D omega signal region kept. " \
                                              "Best candidate selected by smallest |DeltaE_sig|.")
alg_d0_kpi_omega.note(:tag_deltaE, "Tag deltaE window: (-55, +40) MeV for modes with pi0, (-25, +25) MeV otherwise.")
alg_d0_kpi_omega.note(:dt_yield, "Double-tag yield from 2D unbinned ML fit to M_BC^tag vs M_BC^sig. " \
                                  "Peaking backgrounds from K_S0/omega sideband regions subtracted.")
alg_d0_kpi_omega.note(:ks_veto, "D0 -> K_S0 K+ pi- veto: M(pi+pi-) outside (0.478, 0.518) GeV/c^2 for D0->K3Pi tag mode.")

alg_d0_kpi_omega.apply
alg_d0_kpi_omega.execute_on([data_3773, incMC_3773])

# ============================================================
# Algorithm 2: D0 -> K_S0 pi0 omega, omega -> pi+ pi- pi0
# Single tag on anti-D0 (tag_side :D0, 3 modes)
# K_S0 -> pi+ pi-, signal pi0 -> gamma gamma, omega -> pi+ pi- pi0
# Final state: pi+ pi- (K_S0) + gamma gamma (pi0) + pi+ pi- pi0 (omega)
# Total: 4 charged, 4 photons
# ============================================================

alg_d0_kspi0_omega = TagAnalysis.new("D0toKSPi0Omega")
alg_d0_kspi0_omega.set_header(["D0toKSPi0OmegaAlg/D0toKSPi0Omega.h"])
                   .set_constant({ "ECMS" => [:double, 3.773] })

alg_d0_kspi0_omega.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toK3Pi
  t.window :deltaE, abs: 0.025
end

alg_d0_kspi0_omega.signal_side do |s|
  s.photons 4                          # 2 from pi0 + 2 from omega->pi0
  s.charged(pip: 2, pim: 2)           # pi+pi- from K_S0 + pi+pi- from omega
end

alg_d0_kspi0_omega.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_kspi0_omega.note(:ks_reconstruction, "K_S0 -> pi+ pi-. K_S0 mass window: |M_pipi - M_Ks| < 12 MeV/c^2. " \
                                             "K_S0 flight distance > 2*sigma. No PID or Vr requirement on K_S0 pions. " \
                                             "Vz requirement: < 20 cm for K_S0 daughters.")
alg_d0_kspi0_omega.note(:omega_reconstruction, "omega -> pi+ pi- pi0. " \
                                                "Two pi0 in final state, so two possible pi+pi-pi0 combinations. " \
                                                "Omega candidates with M_3pi within +/-25 MeV/c^2 of nominal omega mass. " \
                                                "2D omega signal region used.")
alg_d0_kspi0_omega.note(:tag_deltaE, "Tag deltaE: (-55,+40) MeV for modes with pi0, (-25,+25) MeV otherwise. " \
                                      "Signal deltaE: (-42,+57) MeV (asymmetric due to photon energy losses).")
alg_d0_kspi0_omega.note(:peaking_bg, "Peaking backgrounds A,B,C from K_S0/omega sideband regions studied.")

alg_d0_kspi0_omega.apply
alg_d0_kspi0_omega.execute_on([data_3773, incMC_3773])

# ============================================================
# Algorithm 3: D+ -> K_S0 pi+ omega, omega -> pi+ pi- pi0
# Single tag on D- (tag_side :Dplus, 6 modes)
# Final state: K_S0 (pi+pi-) + pi+ + pi+pi-pi0 (omega)
# Total: 4 charged, 2 photons
# ============================================================

alg_dp_kspi_omega = TagAnalysis.new("DptoKSPiOmega")
alg_dp_kspi_omega.set_header(["DptoKSPiOmegaAlg/DptoKSPiOmega.h"])
                 .set_constant({ "ECMS" => [:double, 3.773] })

alg_dp_kspi_omega.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKSPi, :DptoKPiPiPi0,
          :DptoKSPiPi0, :DptoKS3Pi, :DptoKKPi
  t.window :deltaE, abs: 0.025
end

alg_dp_kspi_omega.signal_side do |s|
  s.photons 2                          # pi0 -> gamma gamma (from omega)
  s.charged(pip: 3, pim: 1)           # pi+(D+), pi+pi-(K_S0), pi+pi-(omega)
  s.require_charge 1                   # D+ charge
end

alg_dp_kspi_omega.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_dp_kspi_omega.note(:ks_reconstruction, "K_S0 -> pi+ pi-. K_S0 mass window: |M_pipi - M_Ks| < 12 MeV/c^2. " \
                                            "K_S0 flight distance > 2*sigma.")
alg_dp_kspi_omega.note(:omega_reconstruction, "omega -> pi+ pi- pi0. Omega mass window +/-25 MeV/c^2. " \
                                               "Two pi0 in final state, two pi+pi-pi0 combinations. 2D omega signal region.")
alg_dp_kspi_omega.note(:tag_deltaE, "Tag deltaE: (-55,+40) MeV for modes with pi0, (-25,+25) MeV otherwise. " \
                                     "Signal deltaE: (-37,+44) MeV.")
alg_dp_kspi_omega.note(:dt_yield, "DT yield from 2D ML fit to M_BC^tag vs M_BC^sig.")

alg_dp_kspi_omega.apply
alg_dp_kspi_omega.execute_on([data_3773, incMC_3773])