# BESIII Analysis: Search for chi_cJ -> e+ e- phi via psi(3686) -> gamma chi_cJ
# Paper: arXiv:2509.23386v1
# Dataset: psi(3686) at 3.686 GeV (709_3686)
#
# Three chi_cJ states (J=0,1,2) share the same final state: e+e- K+K- + gamma
# -> single Algorithm (they are distinguished by the radiative photon energy in ROOT).
#
# Per Rule T3: nominal kinematic fit uses loose chi2_cut 200;
# the paper's tight cut (40) is applied in ROOT.

dm = DatasetManager.new
dm.real_data.find("709_3686")
dm.inclusive_mc.find("709_3686")

alg = Algorithm.new("psip_gamma_chicJ_eePhi")
alg.set_header("Search for chi_cJ -> e+e- phi via psi(3686) -> gamma chi_cJ",
  journal: "arXiv:2509.23386v1",
  dataset: "psi(3686) 709_3686")
alg.set_constant(:ECMS, 3686.1)

sel = Selection.new

sel.select_track(4,
  vxy: [0, 1.0],
  vz: [0, 10.0],
  cos_theta: [0, 0.93])

sel.pid do
  # Electrons: Prob(e) > Prob(pi) and Prob(e) > Prob(K)
  identify_high_momentum_leptons
  # Kaons: Prob(K) > Prob(pi) and Prob(K) > Prob(e)
  identify :kaon, using: :prob, within: [:kp, :km]
end

sel.select_photon(1,
  angle_to_track: 20.0,
  tdc_emc_start: 0,
  tdc_emc_end: 14,
  energy: [0.025, 1.4])

# Post-selection filter: E/p cut for electrons (0.8 < E/p < 1.1)
# and gamma-conversion veto (Rxy < 2 cm)
sel.for_each(:ep) do |ep|
  ep.where("ep.energy / ep.momentum > 0.8 && ep.energy / ep.momentum < 1.1")
end

sel.for_each(:em) do |em|
  em.where("em.energy / em.momentum > 0.8 && em.energy / em.momentum < 1.1")
end

sel.note(:gamma_conversion_veto, "Reject e+e- pairs with Rxy < 2 cm; applied in ROOT")
sel.note(:vertex_fit_chi2, "Vertex fit chi2 < 100 applied in ROOT")

sel.kalman_kinematic_fit(:psip_chiCJ_eePhi,
  constrain: [[:ep, :em, :kp, :km, :gamma0], :four_momentum],
  nominal: true,
  chi2_cut: 200)

alg.with_decay_card(%{
  Decay psi(2S)
  1.0 gamma chi_c0 PHSP;
  1.0 gamma chi_c1 PHSP;
  1.0 gamma chi_c2 PHSP;
  CDecay chi_c0
  1.0 e+ e- phi PHSP;
  EndCDecay
  CDecay chi_c1
  1.0 e+ e- phi PHSP;
  EndCDecay
  CDecay chi_c2
  1.0 e+ e- phi PHSP;
  EndCDecay
  CDecay phi
  1.0 K+ K- VSS;
  EndCDecay
  EndDecay
}).apply(sel)

alg.execute_on(["709_3686"])