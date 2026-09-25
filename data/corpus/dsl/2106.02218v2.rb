# BOSS Ruby DSL: D_s+ → τ+ ν_τ with τ+ → e+ ν_e ν_τ
# Paper: 2106.02218v2 (Phys. Rev. D)
# 6.32 fb^-1 at 4178–4226 MeV, BESIII
# TagAnalysis ST+missing: hadronic ST D_s- (11 modes), signal D_s+ → τ+ ν_τ → e+ ν_e ν_τ
# Missing particle: ν_τ (massless); signal side: 1 positron

# --- Datasets ---
DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_4178 = DatasetManager.real_data.find("703_4178")
data_4189 = DatasetManager.real_data.find("703_4189")
data_4199 = DatasetManager.real_data.find("703_4199")
data_4209 = DatasetManager.real_data.find("703_4209")
data_4219 = DatasetManager.real_data.find("703_4219")
data_4226 = DatasetManager.real_data.find("703_4226")

all_data = [data_4178, data_4189, data_4199, data_4209, data_4219, data_4226]
all_incMC = all_data.map { |d| DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}") }

# Decay card: e+e- → D_s*+ D_s- / D_s+ D_s-
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s+ D_s- PHSP;
  Enddecay
  Decay D_s+
  1.0000 tau+ nu_tau PHSP;
  Enddecay
  Decay tau+
  1.0000 e+ nu_e anti-nu_tau PHSP;
  Enddecay
  End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Dsp_taunu"
  config.related_dataset = data_4226
  config.events          = 500_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# --- TagAnalysis ---
alg = TagAnalysis.new("DspTaunu")
alg.set_header(["DspTaunuAlg/DspTaunu.h"])

# Tag side: D_s- with 11 hadronic modes (all verified in authoritative mode vocabulary)
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK,        # K_S^0 K-
          :DstoKKPi,        # K+K-π-
          :DstoKKPiPi0,     # K+K-π-π0
          :DstoKsKminusPiPi, # K_S^0 K-π+π-
          :DstoKsKplusPiPi,  # K_S^0 K+π-π-
          :DstoPiPiPi,      # π+π-π-
          :DstoPiEta,       # π-η
          :DstoPiPi0Eta,    # π-π0η
          :DstoPiEPPiPiEta, # π-η'_{π+π-η}
          :DstoPiEPRhoGam,  # π-η'_{γρ0}
          :DstoKPiPi        # K-π+π-
  t.charm -1
end

# Signal side: exactly 1 positron, massless missing ν_τ
alg.signal_side do |s|
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_tau
end

# Kinematic fit: tag + e+ + ν_τ = ecms_lab (4C)
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.note(:tag_modes, "11 D_s- hadronic tag modes: K_S^0 K-, K+K-π-, K+K-π-π0, K_S^0 K-π+π-, K_S^0 K+π-π-, π+π-π-, π-η, π-π0η, π-η'(π+π-η), π-η'(γρ0), K-π+π-.")
   .note(:st_selection, "ST D_s- candidates selected via M_ST fits (signal region: [1.89, 2.04] GeV/c^2). Recoil mass M_rec windows vary by energy point for D_s*D_s / D_sD_s selection.")
   .note(:signal_side, "Signal: exactly 1 extra positron (charge opposite to tag). Positron ID: CL_e > 0.1%, CL_e/(CL_e+CL_π+CL_K) > 0.8, p > 0.2 GeV/c, E/p > 0.8. FSR recovery: photons within 5° of e+ added to e+ 4-momentum.")
   .note(:dt_selection, "DT yield from E_extra^tot distribution: signal region E_extra^tot < 0.4 GeV. Backgrounds: non-D_s-, D_s+ → K_L^0 e+ ν_e, D_s+ → X e+ ν_e. BF = N_DT/ε_DT / (Σ N_ST^i/ε_ST^i) / B(τ+→e+ν_eν_τ).")
   .note(:energy_points, "6 energy points: 4178, 4189, 4199, 4209, 4219, ~4226 MeV. ST yields per point; DT yield from combined sample. Recoil mass window varies: [2.050,2.195] at 4178 to [2.040,2.220] at 4226.")
   .note(:bhabha_suppression, "Relative-probability sum cut to suppress Bhabha events: Σ[(CL_e)/(CL_π+CL_K+CL_e)] < 2.0 (n>1) or < 0.9 (n=1).")
   .with_decay_card(decay_card)

alg.apply
alg.execute_on(all_data + all_incMC + [sig_mc])