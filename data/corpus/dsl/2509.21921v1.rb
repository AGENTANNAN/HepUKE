# BESIII Analysis: Search for LNV decay eta -> pi+ pi- e+ e- via J/psi -> phi eta
# Paper: arXiv:2509.21921v1
# Dataset: J/psi at 3.097 GeV (708_3097)
#
# Two independent decay modes (Rule T1):
#   Signal:   J/psi -> phi(->K+K-) eta(->pi+pi-e+e-)
#   Reference: J/psi -> phi(->K+K-) eta(->gamma gamma)
#
# Per Rule T3: nominal kinematic fit uses loose chi2_cut 200;
# the paper's tight cuts (20, 40) are applied in ROOT.

dm = DatasetManager.new
dm.real_data.find("708_3097")
dm.inclusive_mc.find("708_3097")

# ---- Signal Algorithm: eta -> pi+ pi- e+ e- ----

alg_signal = Algorithm.new("Jpsi_phi_eta_pipiee")
alg_signal.set_header("Search for LNV decay eta->pi+pi-e+e- via J/psi->phi eta (signal channel)",
  journal: "arXiv:2509.21921v1",
  dataset: "J/psi 708_3097")
alg_signal.set_constant(:ECMS, 3096.9)

sel_signal = Selection.new
sel_signal.select_track(6,
  vxy: [0, 1.0],
  vz: [0, 10.0],
  cos_theta: [0, 0.93])

sel_signal.pid do
  identify_high_momentum_leptons
  identify :kaon, using: :prob, within: [:kp, :km]
  identify :pion, using: :prob
end

sel_signal.kalman_kinematic_fit(:Jpsiphipiee,
  constrain: [[:pip, :pim, :ep, :em, :kp, :km], :four_momentum],
  nominal: true,
  chi2_cut: 200)

alg_signal.with_decay_card(%{
  Decay J/psi
  1.0 phi eta PHSP;
  CDecay phi
  1.0 K+ K- VSS;
  EndCDecay
  CDecay eta
  1.0 pi+ pi- e+ e- PHSP;
  EndCDecay
  EndDecay
}).apply(sel_signal)

alg_signal.execute_on(["708_3097"])

# ---- Reference Algorithm: eta -> gamma gamma ----

alg_ref = Algorithm.new("Jpsi_phi_eta_gg")
alg_ref.set_header("Search for LNV decay eta->pi+pi-e+e- via J/psi->phi eta (reference channel)",
  journal: "arXiv:2509.21921v1",
  dataset: "J/psi 708_3097")
alg_ref.set_constant(:ECMS, 3096.9)

sel_ref = Selection.new
sel_ref.select_track(2,
  vxy: [0, 1.0],
  vz: [0, 10.0],
  cos_theta: [0, 0.93])

sel_ref.pid do
  identify :kaon, using: :prob, within: [:kp, :km]
end

sel_ref.select_photon(2,
  angle_to_track: 20.0,
  tdc_emc_start: 0,
  tdc_emc_end: 14,
  energy: [0.025, 1.4])

sel_ref.kalman_kinematic_fit(:Jpsiphippgg,
  constrain: [[:kp, :km, :gamma0, :gamma1], :four_momentum],
  nominal: true,
  chi2_cut: 200)

alg_ref.with_decay_card(%{
  Decay J/psi
  1.0 phi eta PHSP;
  CDecay phi
  1.0 K+ K- VSS;
  EndCDecay
  CDecay eta
  1.0 gamma gamma PHSP;
  EndCDecay
  EndDecay
}).apply(sel_ref)

alg_ref.execute_on(["708_3097"])