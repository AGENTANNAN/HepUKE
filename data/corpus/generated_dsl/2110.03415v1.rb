# ============================================================
# psi(3770) 3.773 GeV : real data + inclusive MC
# ============================================================
data_3773  = DatasetManager.real_data.find("712_3773")     # 2.93 fb^-1 psi(3770) data set
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773") # corresponding inclusive MC

# ------------------------------------------------------------
# Signal mode (I): D0 -> K- pi+ omega, omega -> pi+ pi- pi0
# Double tag: recoiling anti-D0 is tagged from the DTag list
# ------------------------------------------------------------
alg_d0_kpi_omega = TagAnalysis.new("D0ToKPiOmegaDTag")
alg_d0_kpi_omega.set_header(["D0ToKPiOmegaDTagAlg/D0ToKPiOmegaDTag.h"])
                 .set_constant({"ECMS" => [:double, 3.773]})
                 .note(:tag_delta_e_window,
                       "tag DeltaE restricted to |DeltaE| < 25 MeV for non-pi0 tag modes and
                        to (-55, +40) MeV for tag modes containing a pi0; mode-dependent, so
                        windowed in ROOT (store-not-cut) instead of a single tag_side window")
                 .note(:best_candidate_selection,
                       "D0 -> K- pi+ omega candidate chosen by smallest |DeltaE_sig|")
                 .note(:ks_veto,
                       "D0 -> K- 3pi tag mode: K_S0 veto excluding M(pi+ pi-) in
                        (0.478, 0.518) GeV/c2")
                 .note(:omega_selection,
                       "omega -> pi+ pi- pi0 uses |M(pi+ pi- pi0) - M(omega)| < 25 MeV/c2,
                        both pi+ pi- pi0 combinations allowed, events kept in the 2D omega
                        signal region")

alg_d0_kpi_omega.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi   # D0 -> K- pi+, K- pi+ pi0, K- 3pi
  t.charm -1                                     # tag the recoiling anti-D0
end

alg_d0_kpi_omega.signal_side do |s|
  s.charged km: 1, pip: 2, pim: 1   # K- pi+ (from D0) + pi+ pi- (from omega)
  s.photons 2                        # pi0 -> gamma gamma
end

alg_d0_kpi_omega.fit do |f|
  f.constrain_four_momentum          # 4C energy-momentum constraint
  f.chi2_cut 200                     # loose chi2 < 200
end

alg_d0_kpi_omega.apply
alg_d0_kpi_omega.execute_on([data_3773, incMC_3773])

# ------------------------------------------------------------
# Signal mode (II): D0 -> K_S0 pi0 omega, omega -> pi+ pi- pi0
# Double tag: recoiling anti-D0 tagged from the DTag list
# ------------------------------------------------------------
alg_d0_ks_pi0_omega = TagAnalysis.new("D0ToKsPi0OmegaDTag")
alg_d0_ks_pi0_omega.set_header(["D0ToKsPi0OmegaDTagAlg/D0ToKsPi0OmegaDTag.h"])
                    .set_constant({"ECMS" => [:double, 3.773]})
                    .note(:tag_delta_e_window,
                          "tag DeltaE restricted to |DeltaE| < 25 MeV for non-pi0 tag modes
                           and to (-55, +40) MeV for tag modes containing a pi0; mode-dependent,
                           windowed in ROOT (store-not-cut)")
                    .note(:signal_delta_e_window,
                          "signal DeltaE window (-42, +57) MeV applied in ROOT (store-not-cut)")
                    .note(:ks_selection,
                          "K_S0 -> pi+ pi- selected with |M(pi+ pi-) - M(K_S0)| < 12 MeV/c2,
                           flight distance > 2 sigma and K_S0-daughter |Vz| < 20 cm; no PID or
                           Vr requirement on these pions")
                    .note(:ks_veto,
                          "D0 -> K- 3pi tag mode: K_S0 veto excluding M(pi+ pi-) in
                           (0.478, 0.518) GeV/c2")
                    .note(:omega_selection,
                          "omega -> pi+ pi- pi0 uses |M(pi+ pi- pi0) - M(omega)| < 25 MeV/c2,
                           both pi+ pi- pi0 combinations allowed, events kept in the 2D omega
                           signal region")
                    .note(:sideband_subtraction,
                          "peaking backgrounds subtracted using K_S0 / omega sidebands; yields
                           from a 2D unbinned ML fit to M_BC^tag vs M_BC^sig")

alg_d0_ks_pi0_omega.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi   # D0 -> K- pi+, K- pi+ pi0, K- 3pi
  t.charm -1                                     # tag the recoiling anti-D0
end

alg_d0_ks_pi0_omega.signal_side do |s|
  s.charged pip: 2, pim: 2   # pi+ pi- (from K_S0) + pi+ pi- (from omega)
  s.photons 4                 # pi0 (from D0) + pi0 (from omega) -> 4 gamma
end

alg_d0_ks_pi0_omega.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_ks_pi0_omega.apply
alg_d0_ks_pi0_omega.execute_on([data_3773, incMC_3773])

# ------------------------------------------------------------
# Signal mode (III): D+ -> K_S0 pi+ omega, omega -> pi+ pi- pi0
# Double tag: recoiling D- tagged from the DTag list
# ------------------------------------------------------------
alg_dp_ks_pi_omega = TagAnalysis.new("DpToKsPiOmegaDTag")
alg_dp_ks_pi_omega.set_header(["DpToKsPiOmegaDTagAlg/DpToKsPiOmegaDTag.h"])
                  .set_constant({"ECMS" => [:double, 3.773]})
                  .note(:tag_delta_e_window,
                        "tag DeltaE restricted to |DeltaE| < 25 MeV for non-pi0 tag modes and
                         to (-55, +40) MeV for tag modes containing a pi0; mode-dependent,
                         windowed in ROOT (store-not-cut)")
                  .note(:signal_delta_e_window,
                        "signal DeltaE window (-37, +44) MeV applied in ROOT (store-not-cut)")
                  .note(:ks_selection,
                        "K_S0 -> pi+ pi- selected with |M(pi+ pi-) - M(K_S0)| < 12 MeV/c2 and
                         flight distance > 2 sigma; no PID or Vr requirement on these pions")
                  .note(:omega_selection,
                        "omega -> pi+ pi- pi0 uses |M(pi+ pi- pi0) - M(omega)| < 25 MeV/c2,
                         both pi+ pi- pi0 combinations allowed, events kept in the 2D omega
                         signal region")
                  .note(:sideband_subtraction,
                        "peaking backgrounds subtracted using K_S0 / omega sidebands; yields
                         from a 2D unbinned ML fit to M_BC^tag vs M_BC^sig")

alg_dp_ks_pi_omega.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0,
          :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi   # D- tag modes
  t.charm -1                                        # tag the recoiling D-
end

alg_dp_ks_pi_omega.signal_side do |s|
  s.charged pip: 3, pim: 2   # pi+ pi- (from K_S0) + pi+ (from D+) + pi+ pi- (from omega)
  s.require_charge 1         # D+ signal side has net charge +1
  s.photons 2                # pi0 (from omega) -> gamma gamma
end

alg_dp_ks_pi_omega.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_dp_ks_pi_omega.apply
alg_dp_ks_pi_omega.execute_on([data_3773, incMC_3773])