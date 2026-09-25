#
# 1903.02359v3: Cross section measurement of e+e- → omega chi_c0
#   at sqrt(s) = 4.178-4.278 GeV (9 energy points)
#   chi_c0 → pi+pi-/K+K-, omega → pi+pi-pi0, pi0 → gamma gamma
# Ordinary analysis with energy scan
#

# ── Dataset: 9 energy points ──────────────────────────────────────────
data_samples = %w[703_4180 703_4190 703_4200 703_4210 703_4220 703_4237 703_4246 703_4270 703_4280]

datasets = data_samples.map { |s| DatasetManager.load_real_data.find(s) }
inc_mc    = data_samples.map { |s| DatasetManager.load_inclusive_mc.find(s) }

# ── Exclusive MC ──────────────────────────────────────────────────────
# Mode I: e+e- → omega chi_c0, chi_c0 → pi+pi-
# Mode II: e+e- → omega chi_c0, chi_c0 → K+K-
# KKMC generator for continuum production

decay_card_mode1 = <<~DECAY
  Decay vpho
  1.0 omega chi_c0  VSS;
  Decay omega
  1.0 pi+ pi- pi0  VVPIPI;
  Decay chi_c0
  1.0 pi+ pi-  SVS;
  Decay pi0
  1.0 gamma gamma  PHSP;
  Enddecay
DECAY

decay_card_mode2 = <<~DECAY
  Decay vpho
  1.0 omega chi_c0  VSS;
  Decay omega
  1.0 pi+ pi- pi0  VVPIPI;
  Decay chi_c0
  1.0 K+ K-  SVS;
  Decay pi0
  1.0 gamma gamma  PHSP;
  Enddecay
DECAY

exclusive_mc_mode1 = DatasetManager.create_exclusive_mc_for(datasets) do |c|
  c.decay_card(decay_card_mode1, generator: :kkmc)
end

exclusive_mc_mode2 = DatasetManager.create_exclusive_mc_for(datasets) do |c|
  c.decay_card(decay_card_mode2, generator: :kkmc)
end

# ── Mode I: chi_c0 → pi+pi- ────────────────────────────────────────
alg_mode1 = Algorithm.new("omega_chi_c0_to_pipi", version: "00-00-01")
alg_mode1.set_header(%w[KinematicFit/KinematicFit.h])
alg_mode1.set_constant(ECMS: "MeasuredEcmsSvc::get_ecms()")
alg_mode1.with_decay_card(decay_card_mode1)
  .note(:generator, "KKMC continuum production — vpho as top mother with ISR")
  .note(:beam_energy, "Energy scan from 4.178 to 4.278 GeV, ECMS obtained run-by-run from MeasuredEcmsSvc")

sel_mode1 = Selection.new("sel_omega_chi_c0_pipi")
sel_mode1.select_track do |t|
  t.cos_theta Range.new(-0.93, 0.93)
  t.Vr 1.0
  t.Vz 10.0
  t.nChrp 4
  t.nChrn 0
  t.nNet 0
end

# Assign high-momentum tracks (>1 GeV/c) as chi_c0 daughters, low-momentum as omega daughters
# NOTE: momentum-based assignment not directly expressible in DSL — ROOT-level selection
sel_mode1.note(:momentum_based_assignment,
  "Tracks with p > 1 GeV/c assigned as chi_c0 daughters (pi+pi- or K+K-); " \
  "tracks with p < 1 GeV/c assigned as omega daughters (pi+pi-). " \
  "Mis-identification rate negligible due to clear momentum separation.")

sel_mode1.assign(:chi_c0_daughter_1, from: :all)
sel_mode1.assign(:chi_c0_daughter_2, from: :all)
sel_mode1.assign(:omega_daughter_1, from: :all)
sel_mode1.assign(:omega_daughter_2, from: :all)

sel_mode1.select_photon do |p|
  p.nGam 2..20
  p.energyThreshold_b 0.025
  p.energyThreshold_e 0.050
  p.tdc_emc_start 0
  p.tdc_emc_end 700
  p.angle_to_track 10.0
end

# Reconstruct pi0 → gamma gamma
sel_mode1.kalman_kinetic_fit(:gamma, :gamma) do |kf|
  kf.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
end

# 5C kinematic fit: 4-momentum conservation + pi0 mass constraint
sel_mode1.kinematic_fit(:pi0, :chi_c0_daughter_1, :chi_c0_daughter_2,
                         :omega_daughter_1, :omega_daughter_2) do |kf|
  kf.constrain_four_momentum
  kf.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  kf.chi2_cut 100
  kf.nominal
end

alg_mode1.apply(sel_mode1)

# ── Mode II: chi_c0 → K+K- ─────────────────────────────────────────
alg_mode2 = Algorithm.new("omega_chi_c0_to_KK", version: "00-00-01")
alg_mode2.set_header(%w[KinematicFit/KinematicFit.h])
alg_mode2.set_constant(ECMS: "MeasuredEcmsSvc::get_ecms()")
alg_mode2.with_decay_card(decay_card_mode2)
  .note(:generator, "KKMC continuum production — vpho as top mother with ISR")
  .note(:beam_energy, "Energy scan from 4.178 to 4.278 GeV, ECMS obtained run-by-run from MeasuredEcmsSvc")

sel_mode2 = Selection.new("sel_omega_chi_c0_KK")
sel_mode2.select_track do |t|
  t.cos_theta Range.new(-0.93, 0.93)
  t.Vr 1.0
  t.Vz 10.0
  t.nChrp 4
  t.nChrn 0
  t.nNet 0
end

sel_mode2.note(:momentum_based_assignment,
  "Same momentum-based assignment: p > 1 GeV/c → chi_c0 daughters, p < 1 GeV/c → omega daughters")

sel_mode2.assign(:chi_c0_daughter_1, from: :all)
sel_mode2.assign(:chi_c0_daughter_2, from: :all)
sel_mode2.assign(:omega_daughter_1, from: :all)
sel_mode2.assign(:omega_daughter_2, from: :all)

sel_mode2.select_photon do |p|
  p.nGam 2..20
  p.energyThreshold_b 0.025
  p.energyThreshold_e 0.050
  p.tdc_emc_start 0
  p.tdc_emc_end 700
  p.angle_to_track 10.0
end

sel_mode2.kalman_kinetic_fit(:gamma, :gamma) do |kf|
  kf.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
end

sel_mode2.kinematic_fit(:pi0, :chi_c0_daughter_1, :chi_c0_daughter_2,
                         :omega_daughter_1, :omega_daughter_2) do |kf|
  kf.constrain_four_momentum
  kf.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  kf.chi2_cut 100
  kf.nominal
end

alg_mode2.apply(sel_mode2)

# ── Hypothesis comparison ──────────────────────────────────────────
# In the BOSS stage both modes are reconstructed; the mode with smaller
# chi2_5C is selected.  This comparison is performed in the ROOT analysis.
alg_mode1.note(:hypothesis_selection,
  "If chi2_5C(pi+pi-) < chi2_5C(K+K-), event assigned to Mode I; " \
  "otherwise assigned to Mode II.  Decision made in ROOT stage by comparing " \
  "the two chi2 values stored in the NTuple.")

# ── Execute ───────────────────────────────────────────────────────────
alg_mode1.execute_on(datasets + inc_mc + exclusive_mc_mode1)
alg_mode2.execute_on(datasets + inc_mc + exclusive_mc_mode2)