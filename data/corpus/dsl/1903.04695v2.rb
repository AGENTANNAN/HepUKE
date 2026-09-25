#
# 1903.04695v2: e+e- → gamma omega J/psi
#   at sqrt(s) = 4.008-4.600 GeV (~11.6 fb⁻¹)
#   J/psi → l+l- (l = e, mu), omega → pi+pi-pi0, pi0 → gamma gamma
# Ordinary analysis with energy scan
#

# ── Dataset: energy scan ──────────────────────────────────────────────
# The paper covers 4.008-4.600 GeV.  Representative datasets from that range:
data_sample_names = %w[
  703_4009 703_4180 703_4190 703_4200 703_4210 703_4220
  703_4230 703_4237 703_4246 703_4260 703_4270 703_4280
  703_4310 703_4360 703_4390 703_4420 703_4470 703_4530 703_4600
]

datasets = data_sample_names.map { |s| DatasetManager.load_real_data.find(s) }
inc_mc   = data_sample_names.map { |s| DatasetManager.load_inclusive_mc.find(s) }

# ── Decay card ─────────────────────────────────────────────────────────
# KKMC: continuum → gamma omega J/psi
# J/psi → l+l- (mixed e/mu), omega → pi+ pi- pi0, pi0 → gamma gamma

decay_card = <<~DECAY
  Decay vpho
  1.0 gamma omega J/psi  PHSP;
  Decay omega
  1.0 pi+ pi- pi0  VVPIPI;
  Decay J/psi
  0.5 e+ e-  VLL;
  0.5 mu+ mu-  VLL;
  Decay pi0
  1.0 gamma gamma  PHSP;
  Enddecay
DECAY

exclusive_mc = DatasetManager.create_exclusive_mc_for(datasets) do |c|
  c.decay_card(decay_card, generator: :kkmc)
end

# ── Algorithm ──────────────────────────────────────────────────────────
alg = Algorithm.new("gamma_omega_jpsi", version: "00-00-01")
alg.set_header(%w[KinematicFit/KinematicFit.h])
alg.set_constant(ECMS: "MeasuredEcmsSvc::get_ecms()")
alg.with_decay_card(decay_card)
  .note(:generator, "KKMC continuum production — vpho as top mother with ISR")
  .note(:beam_energy, "Energy scan from 4.008 to 4.600 GeV, ECMS obtained run-by-run from MeasuredEcmsSvc")

# ── Selection ──────────────────────────────────────────────────────────
sel = Selection.new("sel_gamma_omega_jpsi")

sel.select_track do |t|
  t.cos_theta Range.new(-0.93, 0.93)
  t.Vr 1.0
  t.Vz 10.0
  t.nChrp 4
  t.nChrn 0
  t.nNet 0
end

# Momentum-based lepton/pion separation:
# tracks with p > 1.0 GeV/c → lepton (e/mu from J/psi)
# tracks with p < 1.0 GeV/c → pion (from omega)
# E/p selection separates e from mu in EMC
sel.note(:momentum_based_separation,
  "Tracks with p > 1 GeV/c assigned as lepton candidates (J/psi daughters); " \
  "tracks with p < 1 GeV/c assigned as pion candidates (omega daughters).")
sel.note(:em_separation,
  "EMC energy deposition used to separate e/mu: " \
  "muon: E_EMC < 0.35 GeV; electron: E_EMC > 1.1 GeV. " \
  "Implementation requires custom C++ code in the ROOT stage.")

sel.assign(:lep_p, from: :all)
sel.assign(:lep_m, from: :all)
sel.assign(:pim1, from: :all)
sel.assign(:pip1, from: :all)

sel.select_photon do |p|
  p.nGam 3..20
  p.energyThreshold_b 0.025
  p.energyThreshold_e 0.050
  p.angle_to_track 20.0
end

# Reconstruct pi0 → gamma gamma
sel.kalman_kinetic_fit(:gamma, :gamma) do |kf|
  kf.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
end

# 5C kinematic fit: 4-momentum conservation + pi0 mass constraint
sel.kinematic_fit(:pi0, :pim1, :pip1, :lep_p, :lep_m) do |kf|
  kf.constrain_four_momentum
  kf.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  kf.chi2_cut 100
  kf.nominal
end

# The radiative photon (the "gamma" in gamma omega J/psi) is not explicitly fitted
# here — it is one of the extra photons not used in pi0 reconstruction.
sel.note(:radiative_photon,
  "Radiative photon from e+e- → gamma omega J/psi is among extra photons. " \
  "Combination with smallest chi2 of the 5C fit selected in case of ambiguity.")

# ── Background vetoes ──────────────────────────────────────────────────
# These are applied in ROOT, noted for reference.
alg.note(:psi3686_veto,
  "Events rejected if |RM(pi+pi-) - m[psi(3686)]| < 8 MeV/c^2 " \
  "or |M(pi+pi-J/psi) - m[psi(3686)]| < 7 MeV/c^2.")
alg.note(:etap_jpsi_veto,
  "Events rejected if 0.93 < M(gamma omega) < 0.97 GeV/c^2 and " \
  "M(omega J/psi) > 3.9 GeV/c^2 (eta' J/psi veto).")
alg.note(:omega_chi_c0_bkg,
  "e+e- → omega chi_c0 is irreducible background; estimated from MC " \
  "and included as a component in the M(omega J/psi) fit.")
alg.note(:mass_windows,
  "J/psi mass window: 3.07 < M(l+l-) < 3.14 GeV/c^2. " \
  "omega mass window: 0.72 < M(pi+pi-pi0) < 0.81 GeV/c^2 (asymmetric " \
  "to accommodate X(3872) decay kinematics).")
alg.note(:omega_jpsi_mass,
  "M(omega J/psi) = M(pi+pi-pi0 l+l-) - M(l+l-) + m(J/psi) " \
  "used to partially cancel lepton-pair mass resolution.")

alg.apply(sel)

# ── Execute ───────────────────────────────────────────────────────────
alg.execute_on(datasets + inc_mc + exclusive_mc)