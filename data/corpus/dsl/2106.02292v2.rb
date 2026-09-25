# BOSS Ruby DSL: D0 → ρ- μ+ ν_μ first observation
# Paper: 2106.02292v2 (Phys. Rev. D)
# ψ(3770), 2.93 fb^-1, BESIII
# TagAnalysis ST+missing: hadronic ST D0bar (3 modes), signal D0 → ρ- μ+ ν_μ, ρ- → π-π0
# Missing particle: ν_μ (massless); signal side: π-π0 + μ+

# --- Datasets ---
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card: ψ(3770) → D0 D0bar, signal D0 → ρ- μ+ ν_μ, ρ- → π-π0
decay_card = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay
  Decay D0
  1.0000 rho- mu+ nu_mu PHSP;
  Enddecay
  Decay rho-
  1.0000 pi- pi0 PHSP;
  Enddecay
  End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_D0_rhomunu"
  config.related_dataset = psi3770_data
  config.events          = 1_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# --- TagAnalysis ---
alg = TagAnalysis.new("D0RhoMuNu")
alg.set_header(["D0RhoMuNuAlg/D0RhoMuNu.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })

# Tag side: anti-D0 (charm -1) with 3 hadronic modes
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi,       # K+π-
          :D0toKPiPi0,    # K+π-π0
          :D0toKPiPiPi    # K+π-π-π+
  t.charm -1
end

# Signal side: ρ- → π-π0 and μ+
# π- and μ+ are charged; π0 → γγ requires ≥2 photons
alg.signal_side do |s|
  s.photons 2
  s.charged(pim: 1, mup: 1)
  s.require_charge 0    # ρ-(π-π0) has charge -1, μ+ has +1 → total 0
  s.missing :nu_mu       # massless neutrino
end

# Kinematic fit: tag + signal (π-, μ+, 2γ → π0) + ν_μ = ecms_lab
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg.note(:tag_modes, "3 D0bar hadronic tag modes: K+π-, K+π-π0, K+π-π-π+. ST D0bar yield from M_BC fits: signal MC shape convolved with double-Gaussian + ARGUS background.")
   .note(:st_selection, "ST candidates: ΔE ∈ (-0.055, 0.040) for Kππ0; (-0.025, 0.025) for Kπ and Kπππ. M_BC ∈ (1.859, 1.873) GeV/c^2. Best candidate per mode per charge by min |ΔE|.")
   .note(:signal_side, "Signal: ρ- → π-π0 with M(π-π0) ∈ (0.625, 0.925) GeV/c^2. Muon ID: CL_μ > 0.001, CL_μ > CL_e, CL_μ > CL_K. No CL_μ > CL_π cut (poor μ/π separation at low p). E_μ,EMC ∈ (0.125, 0.275) GeV.")
   .note(:background_vetoes, "K_S^0 veto: M(pi-,pi+) outside (0.458,0.538) and recoil M(pi+,pi-) outside (0.468,0.528). M(ρ-μ+) < 1.5 GeV/c^2. E_extraγ_max < 0.25 GeV. K*(892) veto: |M^2_miss(π→K)| > 0.05 GeV^2/c^4. No extra charged track or π0.")
   .note(:signal_extraction, "DT yield from unbinned ML fit to M_miss^2. Signal: MC shape convolved with Gaussian (resolution params from D0→ρ-e+ν_e control sample). Peaking bkg: D0→π+π-π0π0 from data control sample. BF = N_DT / (N_ST^tot · ε_sig).")
   .note(:bf_result, "B(D0→ρ-μ+ν_μ) = (1.35 ± 0.09_stat ± 0.09_syst) × 10^-3. LFU test: R(μ/e) = 0.90 ± 0.11. Isospin test: Γ(D0→ρ-μ+ν_μ) / 2Γ(D+→ρ0μ+ν_μ) = 0.71 ± 0.14.")
   .note(:muon_pid_note, "Muon PID uses dE/dx + TOF + EMC only (no MUC, since most muons have p < 0.6 GeV/c). CL_μ > 0.001, CL_μ > CL_e, CL_μ > CL_K. E_μ,EMC requirement suppresses ~40% of hadronic background.")
   .with_decay_card(decay_card)

alg.apply
alg.execute_on([psi3770_data, psi3770_incMC, sig_mc])