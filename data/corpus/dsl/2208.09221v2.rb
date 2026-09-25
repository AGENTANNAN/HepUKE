# 2208.09221v2: Search for hyperon ΔS=ΔQ violating decay Ξ⁰ → Σ⁻ e⁺ ν_e
# J/ψ data at √s = 3.097 GeV, double-tag technique
# ST: anti-Ξ⁰ → anti-Λ π⁰ → anti-p π⁺ γγ
# DT: Ξ⁰ → Σ⁻ e⁺ ν_e, Σ⁻ → n π⁻ (n and ν_e missing)

decay_card = <<~DECAYCARD
  Decay J/psi
  1.000 Xi0 anti-Xi0 PHSP;
  Enddecay
  Decay anti-Xi0
  1.000 anti-Lambda pi0 PHSP;
  Enddecay
  Decay anti-Lambda
  1.000 anti-p- pi+ PHSP;
  Enddecay
  Decay Xi0
  1.000 Sigma- e+ nu_e PHSP;
  Enddecay
  Decay Sigma-
  1.000 n pi- PHSP;
  Enddecay
  End
DECAYCARD

algorithm = Algorithm.new("Xi0toSigmaENuESearch", version: '00-00-01')
algorithm.set_header(["Xi0toSigmaENuESearch/Xi0toSigmaENuESearch.h"])
algorithm.set_constant({ "ECMS" => [:double, 3.097] })
algorithm.with_decay_card(decay_card)

# Dataset: J/ψ at 3.097 GeV
data_jpsi = DatasetManager.load_real_data.find("708_3097")
inc_mc_jpsi = DatasetManager.load_inclusive_mc.find("708_3097")

# Exclusive signal MC
signal_mc = DatasetManager.create_exclusive_mc do |c|
  c.name = "Xi0toSigmaENuE_signal"
  c.decay_card = decay_card
  c.related_dataset = data_jpsi
end

# Event selection
event_selection = Selection.new

# Track selection: |cosθ| < 0.93
event_selection.select_track do |t|
  t.cos_theta 0.93
end

# Photon selection: E > 25 MeV barrel, > 50 MeV endcap
event_selection.select_photon do |p|
  p.energyThreshold_b 0.025
  p.energyThreshold_e 0.050
  p.angle_to_track 10.0
end

# Note: EMC time window [0, 700] ns
algorithm.note(:emc_time_window, "EMC shower time required within [0, 700] ns of event start time")

# Reconstruct Λ → p π⁻ (anti-Λ → anti-p π⁺)
# Secondary vertex fit for Λ candidates
event_selection.secondary_vertex_fit(:Lambda) do |v|
  v.build_virtual_particle(:Lambda).from([:prp, :pim])
  v.build_virtual_particle(:anti_Lambda).from([:prm, :pip])
end

# Note: Λ decay length L/σ_L > 2, |M_pπ - m_Λ| < 5 MeV/c²
algorithm.note(:lambda_selection, "Λ: L/σ_L > 2, |M_pπ - m_Λ| < 5 MeV/c²")

# Reconstruct π⁰ → γγ with Kalman kinematic fit
event_selection.kalman_kinematic_fit(:pi0, [:gamma, :gamma]) do |k|
  k.chi2_cut 25
end
event_selection.build_virtual_particle(:pi0, from: [:gamma, :gamma])

# Note: M_γγ ∈ (115, 150) MeV/c², reject both photons in endcap
algorithm.note(:pi0_selection, "M_γγ ∈ (115,150) MeV/c², reject events with both photons from endcap EMC")

# Build Ξ⁰ → Λ π⁰
event_selection.build_virtual_particle(:Xi0_bar, from: [:anti_Lambda, :pi0])

# Note: |M_Λπ⁰ - m_Ξ⁰| < 20 MeV/c², pick best by min|ΔM|
algorithm.note(:xi0_candidate_selection, "|M_Λπ⁰ - m_Ξ⁰| < 20 MeV/c², pick candidate with min|ΔM|")

# DT side: exactly 2 remaining tracks, one positive one negative
event_selection.select_track do |t|
  t.nTot 2
end

# Note: DT tracks - exactly one positive and one negative
algorithm.note(:dt_charge_requirement, "Exactly 2 remaining tracks: one positive, one negative")

# PID for pion (negative track) and positron (positive track)
event_selection.pid(method: :probability) do |pid|
  pid.identify(:pim, against: [:kp, :km, :prp, :prm])
  pid.identify(:ep, against: [:pip, :pim, :kp, :km, :prp, :prm])
  pid.prob_cut 0.001
end

# Note: Pion PID: L(π) > L(K) and L(π) > L(p) and L(π) > 0.001
# Positron PID: L'(e) > 0.001, L'(e)/(L'(e)+L'(π)+L'(K)) > 0.8
algorithm.note(:pid_requirements, "pion: L(π) > L(K) and L(π) > L(p) and L(π) > 0.001; e⁺: L'(e) > 0.001 and L'(e)/(L'(e)+L'(π)+L'(K)) > 0.8")

# Note: χ_dE/dx cuts and momentum cuts
algorithm.note(:chi_dedx_and_momentum, "χ_dE/dx(e as π) < -4.5, χ_dE/dx(π as e) < -2.5; p(π⁻) ∈ (0.20,0.38) GeV/c, p(e⁺) < 0.20 GeV/c")

# Note: Both neutron and neutrino treated as missing particles
# DT yield extracted from M_BC fit on ST side
algorithm.note(:missing_particles, "Both neutron and ν_e treated as missing; DT yield from M_BC fit on ST side")

algorithm.apply(event_selection)
algorithm.execute_on([data_jpsi, inc_mc_jpsi, signal_mc])