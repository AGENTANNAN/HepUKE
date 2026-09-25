#
# 1903.05375v1: eta_c exclusive decays via e+e- → pi+pi- h_c, h_c → gamma eta_c
#   at sqrt(s) = 4.23, 4.26, 4.36, 4.42 GeV
#   4 exclusive modes + inclusive mode
# Ordinary analysis — 4 exclusive decay modes each get own Algorithm (Rule T1)
#

# ── Dataset: 4 energy points ──────────────────────────────────────────
data_sample_names = %w[703_4230 703_4260 703_4360 703_4420]

datasets = data_sample_names.map { |s| DatasetManager.load_real_data.find(s) }
inc_mc   = data_sample_names.map { |s| DatasetManager.load_inclusive_mc.find(s) }

# ── Common note: h_c tagging ───────────────────────────────────────────
# The analysis tags h_c via RM(pi+pi-) → h_c signal region [3.515, 3.535] GeV/c^2,
# then measures eta_c via RM(pi+pi-gamma).  This tagging logic is in ROOT stage.
COMMON_TAG_NOTE = "h_c tagged by RM(pi+pi-) in [3.515, 3.535] GeV/c^2 signal region. " \
  "eta_c identified via RM(pi+pi-gamma) after selecting h_c candidates. " \
  "Simultaneous unbinned maximum likelihood fit to RM(pi+pi-gamma) spectra " \
  "across all 4 energy points for both inclusive and exclusive modes."

# ═══════════════════════════════════════════════════════════════════════
# Mode I: eta_c → K+ K- pi0
# ═══════════════════════════════════════════════════════════════════════
decay_card_kkpi0 = <<~DECAY
  Decay vpho
  1.0 pi+ pi- h_c  PHSP;
  Decay h_c
  1.0 gamma eta_c  VSP_PWAVE;
  Decay eta_c
  1.0 K+ K- pi0  PHSP;
  Decay pi0
  1.0 gamma gamma  PHSP;
  Enddecay
DECAY

exclusive_mc_kkpi0 = DatasetManager.create_exclusive_mc_for(datasets) do |c|
  c.decay_card(decay_card_kkpi0, generator: :kkmc)
end

alg_kkpi0 = Algorithm.new("eta_c_to_KKpi0", version: "00-00-01")
alg_kkpi0.set_header(%w[KinematicFit/KinematicFit.h])
alg_kkpi0.set_constant(ECMS: "MeasuredEcmsSvc::get_ecms()")
alg_kkpi0.with_decay_card(decay_card_kkpi0)
  .note(:h_c_tag, COMMON_TAG_NOTE)

sel_kkpi0 = Selection.new("sel_eta_c_KKpi0")

sel_kkpi0.select_track do |t|
  t.cos_theta Range.new(-0.93, 0.93)
  t.Vr 1.0
  t.Vz 10.0
  t.nChrp 2
  t.nChrn 0
  t.nNet 0
end

sel_kkpi0.pid(method: :probability) do |pid|
  pid.identify(:kp, :km, against: :pip)
  pid.prob_cut 0.001
end

sel_kkpi0.select_photon do |p|
  p.nGam 3..20
  p.energyThreshold_b 0.025
  p.energyThreshold_e 0.050
end

# pi0 → gamma gamma
sel_kkpi0.kalman_kinetic_fit(:gamma, :gamma) do |kf|
  kf.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
end

# 4C kinematic fit: 4-momentum conservation
sel_kkpi0.kinematic_fit(:pi0, :kp, :km) do |kf|
  kf.constrain_four_momentum
  kf.chi2_cut 25   # optimized per energy, ~15-25 (Table II)
  kf.nominal
end

sel_kkpi0.note(:combined_chi2,
  "Combination with min chi2 = chi2_4C + sum(chi2_1C_pi0) + sum(chi2_PID) selected. " \
  "The two pi+pi- tagging pions for h_c are already reconstructed.")

alg_kkpi0.apply(sel_kkpi0)

# ═══════════════════════════════════════════════════════════════════════
# Mode II: eta_c → K_S0 K± pi∓
# ═══════════════════════════════════════════════════════════════════════
decay_card_kskpi = <<~DECAY
  Decay vpho
  1.0 pi+ pi- h_c  PHSP;
  Decay h_c
  1.0 gamma eta_c  VSP_PWAVE;
  Decay eta_c
  1.0 K_S0 K+ pi-  PHSP;
  Decay K_S0
  1.0 pi+ pi-  PHSP;
  Enddecay
DECAY

exclusive_mc_kskpi = DatasetManager.create_exclusive_mc_for(datasets) do |c|
  c.decay_card(decay_card_kskpi, generator: :kkmc)
end

alg_kskpi = Algorithm.new("eta_c_to_KS0_K_pi", version: "00-00-01")
alg_kskpi.set_header(%w[KinematicFit/KinematicFit.h VertexFit/VertexFit.h])
alg_kskpi.set_constant(ECMS: "MeasuredEcmsSvc::get_ecms()")
alg_kskpi.with_decay_card(decay_card_kskpi)
  .note(:h_c_tag, COMMON_TAG_NOTE)

sel_kskpi = Selection.new("sel_eta_c_KS0_K_pi")

sel_kskpi.select_track do |t|
  t.cos_theta Range.new(-0.93, 0.93)
  t.Vr 1.0
  t.Vz 10.0
  t.nChrp 4
  t.nChrn 0
  t.nNet 0
end

# PID: K± identifed against pi, remaining pi track(s)
sel_kskpi.pid(method: :probability) do |pid|
  pid.identify(:kp, against: :pip)
  pid.prob_cut 0.001
end

# K_S0 → pi+ pi- secondary vertex fit
sel_kskpi.secondary_vertex_fit(:K_S0, daughters: [:pip, :pim]) do |sv|
  sv.remove_used_particle_from_candidate_list true
end

sel_kskpi.select_photon do |p|
  p.nGam 1..20
  p.energyThreshold_b 0.025
  p.energyThreshold_e 0.050
end

# 4C kinematic fit
sel_kskpi.kinematic_fit(:K_S0, :kp, :pim) do |kf|
  kf.constrain_four_momentum
  kf.chi2_cut 45   # optimized per energy (Table II)
  kf.nominal
end

sel_kskpi.note(:combined_chi2,
  "Combination with min chi2 = chi2_4C + chi2_vertex(KS0) + sum(chi2_PID) selected.")

alg_kskpi.apply(sel_kskpi)

# ═══════════════════════════════════════════════════════════════════════
# Mode III: eta_c → 2(pi+ pi- pi0)
# ═══════════════════════════════════════════════════════════════════════
decay_card_2pipipi0 = <<~DECAY
  Decay vpho
  1.0 pi+ pi- h_c  PHSP;
  Decay h_c
  1.0 gamma eta_c  VSP_PWAVE;
  Decay eta_c
  1.0 pi+ pi- pi0 pi+ pi- pi0  PHSP;
  Decay pi0
  1.0 gamma gamma  PHSP;
  Enddecay
DECAY

exclusive_mc_2pipipi0 = DatasetManager.create_exclusive_mc_for(datasets) do |c|
  c.decay_card(decay_card_2pipipi0, generator: :kkmc)
end

alg_2pipipi0 = Algorithm.new("eta_c_to_2pipipi0", version: "00-00-01")
alg_2pipipi0.set_header(%w[KinematicFit/KinematicFit.h])
alg_2pipipi0.set_constant(ECMS: "MeasuredEcmsSvc::get_ecms()")
alg_2pipipi0.with_decay_card(decay_card_2pipipi0)
  .note(:h_c_tag, COMMON_TAG_NOTE)

sel_2pipipi0 = Selection.new("sel_eta_c_2pipipi0")

sel_2pipipi0.select_track do |t|
  t.cos_theta Range.new(-0.93, 0.93)
  t.Vr 1.0
  t.Vz 10.0
  t.nChrp 4
  t.nChrn 0
  t.nNet 0
end

# All charged tracks assumed to be pions
sel_2pipipi0.pid(method: :probability) do |pid|
  pid.identify(:pip, against: :kp)
  pid.prob_cut 0.001
end

sel_2pipipi0.select_photon do |p|
  p.nGam 5..20
  p.energyThreshold_b 0.025
  p.energyThreshold_e 0.050
end

# Two pi0 → gamma gamma each
sel_2pipipi0.kalman_kinetic_fit(:gamma, :gamma) do |kf|
  kf.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
end

# 4C kinematic fit with two pi0
sel_2pipipi0.kinematic_fit(:pi0_1, :pi0_2, :pip1, :pim1, :pip2, :pim2) do |kf|
  kf.constrain_four_momentum
  kf.chi2_cut 35   # optimized per energy (Table II)
  kf.nominal
end

sel_2pipipi0.note(:two_pi0,
  "Two pi0 candidates required, each from gamma-gamma pairs. " \
  "DSL v1 kalman_kinematic_fit handles one pi0; second requires custom code.")
sel_2pipipi0.note(:combined_chi2,
  "Combination with min chi2 = chi2_4C + sum(chi2_1C_pi0) + sum(chi2_PID) selected.")

alg_2pipipi0.apply(sel_2pipipi0)

# ═══════════════════════════════════════════════════════════════════════
# Mode IV: eta_c → p pbar
# ═══════════════════════════════════════════════════════════════════════
decay_card_ppbar = <<~DECAY
  Decay vpho
  1.0 pi+ pi- h_c  PHSP;
  Decay h_c
  1.0 gamma eta_c  VSP_PWAVE;
  Decay eta_c
  1.0 p+ anti-p-  SVS;
  Enddecay
DECAY

exclusive_mc_ppbar = DatasetManager.create_exclusive_mc_for(datasets) do |c|
  c.decay_card(decay_card_ppbar, generator: :kkmc)
end

alg_ppbar = Algorithm.new("eta_c_to_ppbar", version: "00-00-01")
alg_ppbar.set_header(%w[KinematicFit/KinematicFit.h])
alg_ppbar.set_constant(ECMS: "MeasuredEcmsSvc::get_ecms()")
alg_ppbar.with_decay_card(decay_card_ppbar)
  .note(:h_c_tag, COMMON_TAG_NOTE)

sel_ppbar = Selection.new("sel_eta_c_ppbar")

sel_ppbar.select_track do |t|
  t.cos_theta Range.new(-0.93, 0.93)
  t.Vr 1.0
  t.Vz 10.0
  t.nChrp 2
  t.nChrn 0
  t.nNet 0
end

# PID: proton and anti-proton
sel_ppbar.pid(method: :probability) do |pid|
  pid.identify(:prp, :prm, against: :pip)
  pid.prob_cut 0.001
end

sel_ppbar.select_photon do |p|
  p.nGam 1..20
  p.energyThreshold_b 0.025
  p.energyThreshold_e 0.050
end

# 4C kinematic fit
sel_ppbar.kinematic_fit(:prp, :prm) do |kf|
  kf.constrain_four_momentum
  kf.chi2_cut 40
  kf.nominal
end

sel_ppbar.note(:combined_chi2,
  "Combination with min chi2 = chi2_4C + sum(chi2_PID) selected.")

alg_ppbar.apply(sel_ppbar)

# ── Inclusive mode ─────────────────────────────────────────────────────
# The inclusive eta_c decay mode requires at least 2 charged tracks and 1 photon,
# with RM(pi+pi-) in h_c signal region and RM(pi+pi-gamma) in [2.52, 3.4] GeV/c^2.
# Not expressible as a full selection chain in DSL v1.
# Tagged via algorithm note.
alg_incl = Algorithm.new("eta_c_inclusive", version: "00-00-01")
alg_incl.set_header(%w[])
alg_incl.set_constant(ECMS: "MeasuredEcmsSvc::get_ecms()")
alg_incl.note(:inclusive_mode,
  "Inclusive eta_c selection: ≥2 charged tracks (all assumed pions), ≥1 photon. " \
  "RM(pi+pi-) in h_c signal region [3.515, 3.535] GeV/c^2. " \
  "RM(pi+pi-gamma) in [2.52, 3.4] GeV/c^2. " \
  "Simultaneous fit across all 4 energies with exclusive modes to extract BFs. " \
  "Charged track multiplicity at production level measured via unfolding method.")

# No Selection for inclusive — too broad to express in DSL

# ── Execute ───────────────────────────────────────────────────────────
alg_kkpi0.execute_on(datasets + inc_mc + exclusive_mc_kkpi0)
alg_kskpi.execute_on(datasets + inc_mc + exclusive_mc_kskpi)
alg_2pipipi0.execute_on(datasets + inc_mc + exclusive_mc_2pipipi0)
alg_ppbar.execute_on(datasets + inc_mc + exclusive_mc_ppbar)
alg_incl.execute_on(datasets + inc_mc)