# Paper: 2201.09710v2 — Partial Wave Analysis of J/ψ → γη'η'
# Type: ORDINARY, J/ψ dataset (708_3097), 10.09×10^9 J/ψ events
# Two independent η' decay modes → two separate Algorithm objects (Rule T1)
# arXiv: https://arxiv.org/abs/2201.09710

# === Dataset ===
data_jpsi = DatasetManager.real_data.find("708_3097")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")

# === Decay Cards ===
# Mode I: J/ψ → γη'η', both η' → ηπ+π-, η → γγ
# Final state: γ + 2(ηπ+π-) = 5γ + 2(π+π-)
decay_card_modeI = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma etap etap PHSP;
  Enddecay
  Decay etap
  1.0000 eta pi+ pi- PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode II: J/ψ → γη'η', η'_1 → γπ+π-, η'_2 → ηπ+π-, η → γγ
# Final state: γ + (γπ+π-) + (ηπ+π-) = 4γ + 2(π+π-)
decay_card_modeII = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma etap1 etap2 PHSP;
  Enddecay
  Decay etap1
  1.0000 gamma pi+ pi- PHSP;
  Enddecay
  Decay etap2
  1.0000 eta pi+ pi- PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# === Exclusive MC ===
sig_mc_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_etap_etap_modeI"
  config.related_dataset = data_jpsi
  config.events          = 500_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

sig_mc_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_etap_etap_modeII"
  config.related_dataset = data_jpsi
  config.events          = 500_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# ============================================================
# === Mode I: both η' → ηπ+π-, η → γγ (5γ + 2(π+π-)) ===
# ============================================================
alg_modeI = Algorithm.new("JpsiGammaEtapEtap_ModeI")
alg_modeI.set_header(["JpsiGammaEtapEtapAlg_ModeI/JpsiGammaEtapEtap_ModeI.h"])
          .set_constant({ "ECMS" => [:double, 3.097] })
          .note(:eta_prime_selection, "Two η' candidates reconstructed from ηπ+π- combinations;
            best pair selected by minimizing |M(ηπ+π-)_1 - m_η'|^2 + |M(ηπ+π-)_2 - m_η'|^2.
            Selected with |M(ηπ+π-) - m_η'| < 0.01 GeV/c^2.")
          .note(:photon_number_veto, "4C fits under 4γ/5γ/6γ hypotheses for photon-number veto:
            events with χ²_5γ < χ²_4γ and χ²_5γ < χ²_6γ accepted.")

sel_modeI = Selection.new

# Charged tracks: 4 charged pions, zero net charge
sel_modeI.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==2"
  nNet "==0"
end

# Photons: ≥5 good photons
sel_modeI.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=5"
end

# PID: all 4 charged tracks identified as pions
sel_modeI.pid(method: :probability) do
  prob_cut 0.001
  identify :pip, against: [:kaon, :proton]
  identify :pim, against: [:kaon, :proton]
  npip "==2"
  npim "==2"
end

# Reconstruct η from γγ pairs (1C Kalman fit)
sel_modeI.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 25
  neta ">=2"
end

# Nominal 6C kinematic fit: 4C + 2 η mass constraints
# Participant list: η(from J/ψ→γηη': both η→γγ) + 4π + radiative γ
# Wait: the fit participants are the FINAL state particles BEFORE η reconstruction
# Since η is reconstructed via kalman, use :eta composite in kinematic fit
# The final state after both η' decays: 5 photons + 2(π+π-)
# Two photons form each η → both η composites → 2:eta + 1 radiative γ + 4π
sel_modeI.kinematic_fit([:gamma, :pip, :pim, :pip, :pim, :eta, :eta]) do
  nominal
  constrain_four_momentum
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 85
end

# Photon-number veto: 4C fits under competing hypotheses (Rule T2)
# 4γ hypothesis (non-nominal)
sel_modeI.kinematic_fit([:gamma, :pip, :pim, :pip, :pim, :eta, :eta]) do
  constrain_four_momentum
end

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on([data_jpsi, incMC_jpsi, sig_mc_modeI])

# ============================================================
# === Mode II: one η' → γπ+π-, other η' → ηπ+π-, η → γγ (4γ + 2(π+π-)) ===
# ============================================================
alg_modeII = Algorithm.new("JpsiGammaEtapEtap_ModeII")
alg_modeII.set_header(["JpsiGammaEtapEtapAlg_ModeII/JpsiGammaEtapEtap_ModeII.h"])
           .set_constant({ "ECMS" => [:double, 3.097] })
           .note(:eta_prime_selection, "η' candidates: one from γπ+π-, one from ηπ+π-.
             Best combination by minimizing weighted |M - m_η'|^2/σ^2 sum.
             Selected: |M(γπ+π-) - m_η'| < 0.02, |M(ηπ+π-) - m_η'| < 0.01 GeV/c^2.")
           .note(:pi0_veto, "|M(γγ) - m_π0| > 0.02 GeV/c^2 for all γγ pairs
             excluding the η pair. Suppresses J/ψ→π0π+π-η' background.")
           .note(:rho_region, "M(π+π-) in η'→γπ+π- required within [0.4, 0.85] GeV/c^2
             to select ρ mass region.")
           .note(:photon_number_veto, "4C fits under 3γ/4γ/5γ hypotheses for photon-number veto.")
           .note(:mass_cut, "M(η'η') < 3.0 GeV/c^2 to suppress fake-photon backgrounds.")

sel_modeII = Selection.new

# Charged tracks: 4 charged pions, zero net charge
sel_modeII.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==2"
  nChrn "==2"
  nNet "==0"
end

# Photons: ≥4 good photons
sel_modeII.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track 10.0
  nGam ">=4"
end

# PID: all 4 charged tracks identified as pions
sel_modeII.pid(method: :probability) do
  prob_cut 0.001
  identify :pip, against: [:kaon, :proton]
  identify :pim, against: [:kaon, :proton]
  npip "==2"
  npim "==2"
end

# Reconstruct η from γγ pairs (1C Kalman fit)
sel_modeII.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 25
  neta ">=1"
end

# Nominal 5C kinematic fit: 4C + 1 η mass constraint
# Final state: 1 radiative γ + 1 η (from η'→ηπ+π- decay) + 1γ (from η'→γπ+π- decay) + 4π
# In the fit: participants = radiative γ + η composite + γ (from η'→γπ+π-) + 4π
# Wait: the two η' decay modes give different particle content:
# η'_1 → γπ+π-: gives 1γ + π+π-
# η'_2 → ηπ+π- (η→γγ): gives (composite η) + π+π-
# Total: 1 radiative γ + (1γ from η'_1 decay) + (η composite from η'_2 decay) + 4π
# Participants: gamma, gamma, eta, pip, pim, pip, pim
sel_modeII.kinematic_fit([:gamma, :gamma, :eta, :pip, :pim, :pip, :pim]) do
  nominal
  constrain_four_momentum
  chi2_cut 75
end

# Photon-number veto: 4C fits under competing hypotheses (Rule T2)
# 4γ hypothesis (non-nominal)
sel_modeII.kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :pip, :pim]) do
  constrain_four_momentum
end

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on([data_jpsi, incMC_jpsi, sig_mc_modeII])