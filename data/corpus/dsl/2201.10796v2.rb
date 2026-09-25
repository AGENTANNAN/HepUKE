# Paper: 2201.10796v2 — Observation of X(2600) in J/ψ → γπ+π-η'
# Type: ORDINARY, J/ψ dataset (708_3097), 10.09×10^9 J/ψ events
# Two independent η' decay modes → two separate Algorithm objects (Rule T1)
# arXiv: https://arxiv.org/abs/2201.10796

# === Dataset ===
data_jpsi = DatasetManager.real_data.find("708_3097")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")

# === Decay Cards ===
# Mode I: J/ψ → γπ+π-η', η' → γπ+π-
# Final state: 2γ + 2(π+π-) = γγπ+π-π+π-
decay_card_modeI = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma pi+ pi- etap PHSP;
  Enddecay
  Decay etap
  1.0000 gamma pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Mode II: J/ψ → γπ+π-η', η' → π+π-η, η → γγ
# Final state: 3γ + 3(π+π-)? No: γ + π+π- + (π+π-η) = γ + 2(π+π-) + η → γ + 2(π+π-) + γγ
# = 3γ + 2(π+π-)
decay_card_modeII = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma pi+ pi- etap PHSP;
  Enddecay
  Decay etap
  1.0000 pi+ pi- eta PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# === Exclusive MC ===
sig_mc_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_pipi_etap_modeI"
  config.related_dataset = data_jpsi
  config.events          = 200_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

sig_mc_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_pipi_etap_modeII"
  config.related_dataset = data_jpsi
  config.events          = 200_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# ============================================================
# === Mode I: η' → γπ+π- (2γ + 2(π+π-)) ===
# ============================================================
alg_modeI = Algorithm.new("JpsiGammaPiPiEtap_ModeI")
alg_modeI.set_header(["JpsiGammaPiPiEtapAlg_ModeI/JpsiGammaPiPiEtap_ModeI.h"])
          .set_constant({ "ECMS" => [:double, 3.097] })
          .note(:high_energy_photons, "two photons required with E > 100 MeV each.")
          .note(:eta_prime_selection, "η' candidate selected from γπ+π- with
            |M(γπ+π-) - m_η'| < 15 MeV/c^2. Best combination by minimum |M-m_η'|.")
          .note(:pi0_eta_veto, "radiative photon paired with all additional photons:
            events with |M(γγ) - m_π0| < 40 MeV/c^2, |M(γγ) - m_η| < 30 MeV/c^2,
            or 720 < M(γγ) < 820 MeV/c^2 (ω→γπ0 veto) rejected.")
          .note(:eta_veto_pipi, "400 < M(γπ+π-) < 563 MeV/c^2 veto suppresses
            J/ψ→γη(η→γπ+π-)π+π- and J/ψ→γη(η→π0π+π-)π+π-.")
          .note(:pi0_etap_veto, "radiative photon paired with extra photons:
            events with |M(γγ) - m_π0| < 15 MeV/c^2 rejected to suppress
            J/ψ→π0π+π-η' background.")

sel_modeI = Selection.new

# Charged tracks: 4 tracks, ≥3 identified as pions
sel_modeI.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
  nNet "==0"
  nTot ">=4"
end

# Photons: ≥2 photons with E>100 MeV
sel_modeI.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.100
  energyThreshold_e 0.100
  angle_to_track 10.0
  nGam ">=2"
end

# PID: ≥3 charged tracks identified as pions
sel_modeI.pid(method: :probability) do
  prob_cut 0.001
  identify :pip, against: [:kaon, :proton]
  identify :pim, against: [:kaon, :proton]
  npip ">=1"
  npim ">=1"
end

# 4C kinematic fit: γγπ+π-π+π- (radiative photon + photon from η' decay + 4π)
# Loose chi2_cut 40 as specified in paper
sel_modeI.kinematic_fit([:gamma, :gamma, :pip, :pim, :pip, :pim]) do
  nominal
  constrain_four_momentum
  chi2_cut 40
end

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on([data_jpsi, incMC_jpsi, sig_mc_modeI])

# ============================================================
# === Mode II: η' → π+π-η, η → γγ (3γ + 2(π+π-)) ===
# ============================================================
alg_modeII = Algorithm.new("JpsiGammaPiPiEtap_ModeII")
alg_modeII.set_header(["JpsiGammaPiPiEtapAlg_ModeII/JpsiGammaPiPiEtap_ModeII.h"])
           .set_constant({ "ECMS" => [:double, 3.097] })
           .note(:high_energy_photons, "three photons required with E > 100 MeV each.")
           .note(:eta_reconstruction, "η reconstructed with |M(γγ) - m_η| < 30 MeV/c^2.")
           .note(:pi0_veto, "|M(γγ) - m_π0| > 40 MeV/c^2 for all photon pairs to suppress π0 background.")
           .note(:eta_prime_selection, "η' candidate from π+π-η with |M(π+π-η) - m_η'| < 10 MeV/c^2.
             Best combination by minimum |M - m_η'|.")
           .note(:pi0_etap_veto, "radiative photon paired with extra photons:
             events with |M(γγ) - m_π0| < 15 MeV/c^2 rejected to suppress
             J/ψ→π0π+π-η' background.")

sel_modeII = Selection.new

# Charged tracks: 4 tracks, ≥3 identified as pions
sel_modeII.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
  nNet "==0"
  nTot ">=4"
end

# Photons: ≥3 photons
sel_modeII.select_photon do
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.100
  energyThreshold_e 0.100
  angle_to_track 10.0
  nGam ">=3"
end

# PID: ≥3 charged tracks identified as pions
sel_modeII.pid(method: :probability) do
  prob_cut 0.001
  identify :pip, against: [:kaon, :proton]
  identify :pim, against: [:kaon, :proton]
  npip ">=1"
  npim ">=1"
end

# Reconstruct η from γγ (1C Kalman fit)
sel_modeII.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 25
  neta ">=1"
end

# Nominal 5C kinematic fit: 4C + η mass constraint
# Final state: 1 radiative γ + 2γ (from η→γγ) + 4π → use eta composite + 1 radiative γ + 4π
# Participants: gamma, gamma, gamma, pip, pim, pip, pim
# But: we reconstructed η via kalman → use :eta composite
# Remaining particles: 1 radiative γ + :eta + 4π
sel_modeII.kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim, :pip, :pim]) do
  nominal
  constrain_four_momentum
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 40
end

# Competing 4C fit (without η mass constraint) for background suppression (Rule T2)
sel_modeII.kinematic_fit([:gamma, :gamma, :gamma, :pip, :pim, :pip, :pim]) do
  constrain_four_momentum
end

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on([data_jpsi, incMC_jpsi, sig_mc_modeII])