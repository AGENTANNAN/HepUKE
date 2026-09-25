# DSL for 2212.09048v1: D^0 → K_L^0 π^+ π^- amplitude analysis
# D^0 DT at ψ(3770) with hadronic flavor tags
# Two signal modes: K_S^0π^+π^- (fully reconstructed) and K_L^0π^+π^- (missing K_L^0)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================
# Signal Mode 1: D^0 → K_L^0 π^+ π^-  (ST + missing K_L^0)
# Tag side: anti-D^0 → K^+π^-, K^+π^-π^+π^-, K^+π^-π^0
# Signal side: π^+π^- + missing K_L^0
# ============================================================

decay_card_kl = <<~DECAY
Decay psi(3770)
  1.0  D0  anti-D0  VSS;
Enddecay
Decay D0
  1.0  K_L0  pi+  pi-  PHSP;
Enddecay
Decay anti-D0
  1.0  K+  pi-  PHSP;
Enddecay
DECAY

sigMC_kl = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_D0_KLpipi"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_kl
  config.cross_section   = :default
end

alg_kl = TagAnalysis.new("D0TagKL0PiPi")
alg_kl.set_header(["D0TagKL0PiPiAlg/D0TagKL0PiPi.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })

alg_kl.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPiPi, :D0toKPiPi0
  t.charm -1
end

alg_kl.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.missing :K_L0
  s.require_charge 0
end

alg_kl.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_kl.note(:signal_mode, "K_L^0 reconstructed via missing-mass technique: M_miss^2 = (E_beam - E_π+ - E_π-)^2 - |p_tag + p_π+ + p_π-|^2")
alg_kl.note(:pi0_eta_veto, "π^0 and η vetoes: reject events with any γγ pair in [0.095,0.165] GeV/c^2 (π^0) or [0.48,0.58] GeV/c^2 (η)")
alg_kl.note(:dcs_tag, "DCS contamination in hadronic flavor tags accounted for in the amplitude fit via coherent amplitude addition")
alg_kl.note(:amplitude_analysis, "This is an amplitude analysis. The ROOT stage performs the full isobar-model fit with U-spin breaking parameters rho_hat. The BOSS stage performs event selection and kinematic fit.")
alg_kl.note(:ks_peaking_bkg, "K_S^0π^+π^- peaking background (~5%) modeled in the amplitude fit; K_S^0 mass window cut [0.485,0.510] GeV/c^2 on tagged Kπππ mode")
alg_kl.note(:background_model, "Non-peaking background (~6%) from sideband of M_miss^2 modeled with Gaussian kernel estimator")
alg_kl.note(:tag_selection, "Tag-side ΔE ±3σ and M_BC ±3σ windows applied. Best tag candidate selected by minimum |ΔE|")

alg_kl.apply
alg_kl.execute_on([data_3773, incMC_3773, sigMC_kl])

# ============================================================
# Signal Mode 2: D^0 → K_S^0 π^+ π^-  (fully reconstructed DT)
# Tag side: anti-D^0 → K^+π^-, K^+π^-π^+π^-, K^+π^-π^0
# Signal side: K_S^0 → π^+π^- + π^+π^- (4 charged tracks total)
# ============================================================

decay_card_ks = <<~DECAY
Decay psi(3770)
  1.0  D0  anti-D0  VSS;
Enddecay
Decay D0
  1.0  K_S0  pi+  pi-  PHSP;
Enddecay
Decay anti-D0
  1.0  K+  pi-  PHSP;
Enddecay
Decay K_S0
  1.0  pi+  pi-  PHSP;
Enddecay
DECAY

sigMC_ks = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_D0_KSpipi"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_ks
  config.cross_section   = :default
end

alg_ks = Algorithm.new("D0TagKS0PiPi", "00-00-01")
alg_ks.set_header(["D0TagKS0PiPiAlg/D0TagKS0PiPi.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .with_decay_card(decay_card_ks)

sel_ks = Selection.new

sel_ks.select_track do |t|
  t.nChrp 4
  t.nTot 4
end

sel_ks.select_photon do |p|
  p.nGam 0
end

sel_ks.secondary_vertex_fit(:K_S0, daughters: [:pip, :pim]) do |v|
  v.mass_window [0.485, 0.510]
  v.flight_significance 2.0
end

sel_ks.build_virtual_particle(:D0_sig, from: [:K_S0, :pip, :pim])

sel_ks.kinematic_fit([:K_S0, :pip, :pim, :kp, :pim]) do |fit|
  fit.constrain_four_momentum
  fit.invariant_mass_of(:K_S0, :pip, :pim).constrain_to_nominal_mass_of(:D0)
  fit.chi2_cut 200
  fit.nominal
end

sel_ks.note(:ks_reconstruction, "K_S^0 candidates: secondary vertex fit, flight significance L/σ_L > 2, mass window [0.485,0.510] GeV/c^2. Events with multiple K_S^0 candidates rejected to suppress D→K_S^0K_S^0X")
sel_ks.note(:kinematic_fit, "6C kinematic fit: total 4-momentum conservation + D^0 mass constraint + K_S^0 mass constraint. Events where fit doesn't converge are discarded")
sel_ks.note(:tag_selection, "Hadronic flavor tag modes D^0→K^+π^-, K^+π^-π^+π^-, K^+π^-π^0. Tag-side ΔE ±3σ and M_BC ±3σ windows. Best tag candidate: minimum |ΔE|")
sel_ks.note(:amplitude_analysis, "The ROOT stage performs the full amplitude fit with isobar model. BOSS stage performs event selection, K_S vertex fit, and 6C kinematic fit")
sel_ks.note(:pi0_eta_veto, "π^0 veto: reject γγ pairs in [0.095,0.165] GeV/c^2. η veto: reject γγ pairs in [0.48,0.58] GeV/c^2")

alg_ks.apply(sel_ks)
alg_ks.execute_on([data_3773, incMC_3773, sigMC_ks])