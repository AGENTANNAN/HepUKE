# DSL for 2311.07043v2: Partial Wave Analysis of J/ψ → φπ0η
# Ordinary analysis: Algorithm + Selection
# J/ψ → φ(→K+K-) π0(→γγ) η(→γγ)
# 10B J/ψ events at 3.097 GeV

# ============================================================
# Datasets
# ============================================================
data_jpsi = DatasetManager.load_real_data.find("708_3097")
inc_jpsi  = DatasetManager.load_inclusive_mc.find("708_3097")

# Signal MC: J/ψ → φπ0η, φ→K+K-, π0→γγ, η→γγ
sig_mc = DatasetManager.create_exclusive_mc do |c|
  c.decay_card = <<~DECAY
    Decay J/psi
      1.0 phi pi0 eta PHSP;
    Enddecay
    Decay phi
      1.0 K+ K- VSS;
    Enddecay
    Decay pi0
      1.0 gamma gamma PHSP;
    Enddecay
    Decay eta
      1.0 gamma gamma PHSP;
    Enddecay
  DECAY
  c.related_dataset [data_jpsi]
end

# ============================================================
# Algorithm + Selection
# ============================================================
alg = Algorithm.new("Jpsi2PhiPi0Eta")
alg.set_header(["Jpsi2PhiPi0Eta/Jpsi2PhiPi0Eta.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })

event_selection = Selection.new("Jpsi2PhiPi0EtaSel")

# Step 1: Charged tracks — exactly 2 kaons (φ→K+K-)
# |cosθ| < 0.93, Vz < 10 cm, Vxy < 1 cm
event_selection.select_track do |t|
  t.nChrp 1
  t.nChrn 1
  t.nTot 2
  t.cos_theta_range(-0.93, 0.93)
  t.vz_range(-10.0, 10.0)
  t.vxy_range(0.0, 1.0)
end

# Step 2: Photons — at least 4 (π0→γγ + η→γγ)
# E > 25 MeV barrel, > 50 MeV endcap, TDC within 700 ns
event_selection.select_photon do |ph|
  ph.n_min 4
  ph.min_energy_barrel 0.025
  ph.min_energy_endcap 0.050
end

# Step 3: PID — identify kaons
event_selection.pid(method: :probability) do |pid|
  pid.identify :kp, :km, against: :pim
  pid.prob_cut 0.001
end

# Step 4: 4C kinematic fit to J/ψ → K+K-γγγγ hypothesis
# χ²_4C < 45; best combination kept for events with >4 photons
event_selection.kinematic_fit([:kp, :km, :gamma, :gamma, :gamma, :gamma]) do |kf|
  kf.chi2_cut 45
  kf.nominal
end

# Step 5: Reconstruct π0 and η via 1C mass-constrained fits
event_selection.kalman_kinematic_fit([:gamma, :gamma]) do |kf|
  kf.mass_constraint :pi0
  kf.name :pi0_fit
end
event_selection.kalman_kinematic_fit([:gamma, :gamma]) do |kf|
  kf.mass_constraint :eta
  kf.name :eta_fit
end

# Step 6: 6C kinematic fit with π0 and η mass constraints
# Enforces energy-momentum conservation + nominal π0 and η masses
event_selection.kinematic_fit([:kp, :km, :gamma, :gamma, :gamma, :gamma]) do |kf|
  kf.mass_constraint :pi0, :eta
  kf.nominal
end

# Post-selection: apply mass windows and select best candidate
event_selection.for_each(:phi_candidate) do |fe|
  fe.define(:m_kk) { (p4(:kp) + p4(:km)).m }
  fe.where { m_kk > 1.01 && m_kk < 1.03 }
end

alg.note(:photon_combination,
  "For events with >4 photons, all γγ combinations are considered. " \
  "The combination minimizing χ²_{π0η} = (M_{γ1γ2}-m_π0)²/σ²_π0 + (M_{γ3γ4}-m_η)²/σ²_η is kept. " \
  "Not expressible in DSL v1 — implemented via manual combination loop.")
alg.note(:pi0_eta_mass_windows,
  "π0: |M_{γ1γ2} - m_π0| < 15 MeV/c². η: |M_{γ3γ4} - m_η| < 25 MeV/c². " \
  "Applied after the 4C fit, before the 6C fit.")
alg.note(:pi0pi0_veto,
  "Veto J/ψ→K+K-π0π0 background: require χ²_{π0π0} > 90, " \
  "where χ²_{π0π0} = (M_{γ1γ2}-m_π0)²/σ²_π0 + (M_{γ3γ4}-m_π0)²/σ²_π0.")
alg.note(:etaeta_veto,
  "Veto J/ψ→K+K-ηη background: require χ²_{ηη} > 8.")
alg.note(:qfactor,
  "Q-factor method applied for non-π0/non-η background subtraction. " \
  "Each event assigned a weight Q_i from 2D fit to M_{γ1γ2} vs M_{γ3γ4}.")
alg.note(:pwa,
  "Partial Wave Analysis performed with GPUPWA framework on the selected events. " \
  "PWA fit includes φ(1680), h1(1900), X(2000), a0(980)-f0(980) mixing, and " \
  "non-φ background modeled from sideband PWA. Beyond DSL scope.")

alg.with_decay_card(<<~DECAY)
  Decay J/psi
    1.0 phi pi0 eta PHSP;
  Enddecay
  Decay phi
    1.0 K+ K- VSS;
  Enddecay
  Decay pi0
    1.0 gamma gamma PHSP;
  Enddecay
  Decay eta
    1.0 gamma gamma PHSP;
  Enddecay
DECAY
alg.apply(event_selection)
alg.execute_on([data_jpsi, inc_jpsi, sig_mc])