# Paper: 2307.03024v2
# Analysis: Ds+ -> K+ K- mu+ nu_mu (double-tag, semileptonic, form factors)
# Energy: 8 points from 4.128 to 4.226 GeV
# Method: Tag-based analysis (TagAnalysis) — single tag Ds- + signal Ds+ semileptonic

# --- Datasets (8 energy points) ---
data_4128 = DatasetManager.real_data.find("705_4130")
data_4157 = DatasetManager.real_data.find("705_4160")
data_4178 = DatasetManager.real_data.find("703_4180")
data_4189 = DatasetManager.real_data.find("703_4190")
data_4199 = DatasetManager.real_data.find("703_4200")
data_4209 = DatasetManager.real_data.find("703_4210")
data_4219 = DatasetManager.real_data.find("703_4220")
data_4226 = DatasetManager.real_data.find("703_4230")

all_data = [data_4128, data_4157, data_4178, data_4189, data_4199, data_4209, data_4219, data_4226]

# Inclusive MC at 4.178 GeV (40x lumi)
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")

# --- Decay card for signal MC ---
decay_card = <<~DECAYCARD
  Decay D_s+
  1.000 K+ K- mu+ nu_mu PHSP;
  Enddecay
  End
DECAYCARD

# Signal MC per energy point
sig_mc_samples = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name = "sig_Dsp_KKmunu"
  config.events = 200_000
  config.decay_card = decay_card
  config.cross_section = :default
end

# --- TagAnalysis ---
algorithm = TagAnalysis.new("DsKKMuNu")
algorithm.set_header(["DsKKMuNuAlg/DsKKMuNu.h"])
  .set_constant({ "ECMS" => [:double, 4.178] })
  .with_decay_card(decay_card)
  .note(:analysis_method, "Double-tag (DT) method via single-tag Ds- + semileptonic Ds+ signal reconstruction; Ds mesons produced via e+e- -> Ds*+ [-> gamma/pi0 Ds+] Ds-")
  .note(:ds_star_transition, "Transition gamma or pi0 from Ds*+ decay selected among unused photons; if multiple candidates, the one giving minimum |DeltaE| is kept. This is part of DTagAlg internal DT reconstruction.")
  .note(:kinematic_fit, "3C kinematic fit applied: energy-momentum conservation + both Ds masses constrained to nominal Ds mass + Ds gamma/pi0 mass constrained to nominal Ds* mass. The DSL fit block approximates with 4C only; the Ds* constraint is applied at the ROOT level.")
  .note(:umiss_signal, "Umiss = (ECM - E_Ds- - E_gamma/pi0 - E_KK - E_mu+) - |p_Ds- + p_gamma/pi0 + p_KK + p_mu+| used as signal variable; signal shape from MC convolved with Gaussian")
  .note(:muon_pid, "Muon candidates satisfy L_mu > L_K, L_mu > L_e, L_mu > 0.001 based on combined MDC dE/dx + TOF + EMC likelihoods")
  .note(:background_vetoes, "M(KKnu_mu) > 1.30 GeV/c^2 to suppress Ds+ -> KK pi+ peaking background; M(KKmu+) < 1.75 GeV/c^2 to reject remaining Ds+ -> KK pi+; E_gamma_extra_max < 0.2 GeV to suppress Ds+ -> KK pi+ pi0")
  .note(:extra_track_veto, "Zero unused charged tracks and zero unused pi0 candidates required in all DT candidate events")
  .note(:pwa, "Partial wave analysis extracts form factor ratios rV = V(0)/A1(0) and r2 = A2(0)/A1(0); phi(1020) meson is the dominant KK resonance; S-wave f0(980) component searched for")
  .note(:phi_bf, "Result: BF(Ds+ -> phi mu+ nu_mu) = (2.25 +/- 0.09 +/- 0.07) x 10^{-2}; ratio to phi e+ nu_e = 0.94 +/- 0.08")
  .note(:tag_modes_mapped, "14 ST Ds- hadronic modes mapped to DTagAlg channel names: DstoKKPi, DstoKPiPi, DstoPiPiPi, DstoKKPiPi0, DstoPiEPRhoGam, DstoPiPi0Eta, DstoKsKminusPiPi, DstoKsKplusPiPi, DstoPiEta, DstoKsKsPi, DstoPiPiPiEta, DstoPiEPPiPiEta, DstoKsKPi0, DstoKsK")

# Tag side: Ds- reconstructed in 14 hadronic modes
algorithm.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKPiPi, :DstoPiPiPi, :DstoKKPiPi0,
          :DstoPiEPRhoGam, :DstoPiPi0Eta,
          :DstoKsKminusPiPi, :DstoKsKplusPiPi,
          :DstoPiEta, :DstoKsKsPi,
          :DstoPiPiPiEta, :DstoPiEPPiPiEta,
          :DstoKsKPi0, :DstoKsK
  t.charm -1
end

# Signal side: K+ K- mu+ + missing nu_mu
algorithm.signal_side do |s|
  s.charged(kp: 1, km: 1, mup: 1)
  s.require_charge 1          # Ds+ has charge +1: K+(+1) K-(-1) mu+(+1) = +1
  s.missing :nu_mu            # massless form — muon neutrino
end

# 4C kinematic fit (D_s* constraint applied at ROOT level)
algorithm.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

algorithm.apply
algorithm.execute_on(all_data + [incMC_4178] + sig_mc_samples)