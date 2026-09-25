# DSL for: Search for J/ψ → D⁻ e⁺ νₑ + c.c.
# Paper: 2104.06628v3
# Classification: Ordinary (Algorithm + Selection), partial_rec with missing neutrino
# Data: √s = 3.097 GeV, J/ψ
# Process: J/ψ → D⁻ e⁺ νₑ, D⁻ → K⁺ π⁻ π⁻

data_3097 = DatasetManager.real_data.find("708_3097")
incMC_3097 = DatasetManager.inclusive_mc.find("708_3097")

sig_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Jpsi_to_Dm_ep_nue"
  config.related_dataset = data_3097
  config.events          = 500_000
  config.decay_card      = <<~DECAYCARD
    Decay J/psi
    1.0000 D- e+ nu_e VLL;
    Enddecay
    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay
    End
  DECAYCARD
  config.cross_section   = :default
end

alg = Algorithm.new("JpsiToDmEpNue", version: "00-00-01")
alg.set_header(["JpsiToDmEpNueAlg/JpsiToDmEpNue.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })

sel = Selection.new("DmEpNueSelection")

# --- Track selection: exactly 4 charged tracks, zero net charge ---
sel.select_track do |t|
  t.cos_theta 0.93
  t.Vz 10.0
  t.Vr 1.0
  t.nTot 4
  t.nNet 0
end

# --- PID: kaon, pion, positron identification ---
# Kaon vs pion: L(K) > L(π) for kaon, L(π) > L(K) for pion
# Positron: combined likelihood P_e > 0.001, P_e/(P_π + P_K) > 1/4
sel.pid(method: :probability) do |p|
  p.identify :kp, against: [:pip]
  p.identify :pim, against: [:kp]
  p.identify :ep, against: [:pip, :kp]
  p.prob_cut 0.001
end

# High-momentum lepton identification for positron
sel.note(:positron_pid_detail,
  "Positron PID: combined likelihood from MDC, TOF, EMC. P_e > 0.001 AND P_e/(P_π + P_K) > 1/4 (i.e., factor > 4 per paper text). Additional E/p requirement: 0.85 < E/p < 1.05 applied in ROOT analysis.")

# --- Photon selection: suppress background with extra photons ---
sel.select_photon do |g|
  g.energyThreshold_b 0.025
  g.energyThreshold_e 0.050
  g.angle_to_track :kp, 10.0
  g.angle_to_track :pim, 10.0
  g.angle_to_track :ep, 10.0
end

sel.note(:photon_timing,
  "EMC time within [0, 700] ns of event start time required for all photon candidates.")

sel.note(:total_photon_energy_cut,
  "Total energy of all good photons E_γ_tot < 0.2 GeV to suppress backgrounds with extra photons. Applied in ROOT analysis.")

# --- D⁻ reconstruction from K⁺ π⁻ π⁻ ---
sel.build_virtual_particle(:Dm, from: [:kp, :pim, :pim]) do |v|
  v.remove_used_particle_from_candidate_list false
end

# D⁻ mass window: [1.85, 1.89] GeV/c² (±3σ around nominal D⁻ mass)
sel.note(:Dm_mass_window,
  "D⁻ candidate M_Kππ in [1.85, 1.89] GeV/c² applied in ROOT analysis.")

# --- 1C kinematic fit: constrain K⁺ π⁻ π⁻ invariant mass to D⁻ mass ---
sel.kinematic_fit([:kp, :pim, :pim]) do |k|
  k.invariant_mass_of(:kp, :pim, :pim).constrain_to_nominal_mass_of(:Dm)
  k.chi2_cut 10
end

# Note: this is a 1C fit (not 4C) — only D⁻ mass constrained
# 4C momentum constraint is NOT used here because neutrino is missing

# --- Partial reconstruction: missing neutrino ---
# Umis = Emiss - c|p_miss| peaks at 0 for signal
# |p_miss| > 50 MeV/c required
sel.partial_rec([:Dm, :ep]) do |pr|
  pr.partial_miss([:nu_e])
end

sel.note(:missing_momentum_cut,
  "|p_miss| > 50 MeV/c required to suppress hadronic backgrounds with mis-identified positrons. U_miss = E_miss - c|p_miss| used to extract signal yield in ROOT.")

sel.note(:Umis_fit,
  "Signal yield extracted via unbinned maximum likelihood fit to U_miss distribution in ROOT. Signal shape from MC, background modeled with linear function. Upper limit set via Bayesian approach at 90% CL.")

sel.note(:signal_mc_model,
  "Signal MC generated assuming weak interaction via c→d charged current, ignoring hadronization and quark spin-flip effects. Systematic uncertainty 3.0% from comparison with phase space model.")

alg.with_decay_card(nil)  # decay card embedded in create_exclusive_mc
alg.apply(sel)
alg.execute_on([data_3097, incMC_3097, sig_mc])