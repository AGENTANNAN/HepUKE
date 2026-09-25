# Paper: 2202.00621v4 — Observation of η1(1855) with exotic JPC=1-+ in J/ψ → γηη'
# Type: ORDINARY, J/ψ dataset (708_3097), 10.09×10^9 J/ψ events
# Two η' decay modes → two Algorithm objects (Rule T1)
# Selection from companion paper Phys. Rev. D 106, 072012 (2022)
# arXiv: https://arxiv.org/abs/2202.00621

# === Dataset ===
data_jpsi = DatasetManager.real_data.find("708_3097")
incMC_jpsi = DatasetManager.inclusive_mc.find("708_3097")

# === Decay Cards ===
# Mode I: J/ψ → γηη', η → γγ, η' → ηπ+π- (η→γγ)
# Final state: γ + (γγ) + ((γγ)π+π-) = 5γ + π+π-
decay_card_modeI = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta etap PHSP;
  Enddecay
  Decay etap
  1.0000 eta pi+ pi- PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode II: J/ψ → γηη', η → γγ, η' → γπ+π-
# Final state: γ + (γγ) + (γπ+π-) = 4γ + π+π-
# Use Alias to distinguish the two η' decay daughters
decay_card_modeII = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta etap PHSP;
  Enddecay
  Decay etap
  1.0000 gamma pi+ pi- PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# === Exclusive MC ===
sig_mc_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_eta_etap_modeI"
  config.related_dataset = data_jpsi
  config.events          = 500_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

sig_mc_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_eta_etap_modeII"
  config.related_dataset = data_jpsi
  config.events          = 500_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# ============================================================
# === Mode I: η' → ηπ+π- (both η→γγ) — 5γ + π+π- ===
# ============================================================
# Selection from companion paper PRD 106, 072012 (Ref. [27])
# Final state: γ_radiative + 2γ(η) + 2γ(η from η') + π+π-
alg_modeI = Algorithm.new("JpsiGammaEtaEtap_ModeI")
alg_modeI.set_header(["JpsiGammaEtaEtapAlg_ModeI/JpsiGammaEtaEtap_ModeI.h"])
          .set_constant({ "ECMS" => [:double, 3.097] })
          .note(:selection_reference, "Selection criteria from companion paper
            Phys. Rev. D 106, 072012 (2022) [Ref. 27]. 4788 events in ηπ+π- mode.")
          .note(:eta_prime_selection, "η' candidates selected from ηπ+π- with
            invariant mass window. Best combination selection applied.")
          .note(:sideband, "Non-ηη' background estimated from η' sideband events in data.")
          .note(:phi_veto, "J/ψ → φη', φ → γη rejected; events around M(γη) ~ 1.02 GeV/c^2
            depleted. Applied in PWA not at BOSS selection level.")

sel_modeI = Selection.new

# Charged tracks: 2 pions with opposite charge
sel_modeI.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==1"
  nChrn "==1"
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

# PID: pions identified
sel_modeI.pid(method: :probability) do
  prob_cut 0.001
  identify :pip, against: [:kaon, :proton]
  identify :pim, against: [:kaon, :proton]
  npip "==1"
  npim "==1"
end

# Reconstruct η from γγ pairs (1C Kalman fit) — need 2 η candidates
sel_modeI.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 25
  neta ">=2"
end

# 5C kinematic fit: 4C + η mass constraint
# Final state: 5γ + π+π- with two η composites from kalman fit
# Participants: :gamma (radiative) + :eta (from J/ψ) + :eta (from η' decay) + π+π-
sel_modeI.kinematic_fit([:gamma, :eta, :eta, :pip, :pim]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on([data_jpsi, incMC_jpsi, sig_mc_modeI])

# ============================================================
# === Mode II: η' → γπ+π- — 4γ + π+π- ===
# ============================================================
# Selection from companion paper PRD 106, 072012 (Ref. [27])
# Final state: γ_radiative + 2γ(η→γγ) + γ(η'→γπ+π-) + π+π- = 4γ + π+π-
alg_modeII = Algorithm.new("JpsiGammaEtaEtap_ModeII")
alg_modeII.set_header(["JpsiGammaEtaEtapAlg_ModeII/JpsiGammaEtaEtap_ModeII.h"])
           .set_constant({ "ECMS" => [:double, 3.097] })
           .note(:selection_reference, "Selection criteria from companion paper
             Phys. Rev. D 106, 072012 (2022) [Ref. 27]. 10544 events in γπ+π- mode.")
           .note(:eta_prime_selection, "η' candidates selected from γπ+π- with
             invariant mass window. Best combination selection applied.")
           .note(:phi_veto, "J/ψ → φη', φ → γη rejected; events around M(γη) ~ 1.02 GeV/c^2
             depleted. Applied in PWA not at BOSS selection level.")

sel_modeII = Selection.new

# Charged tracks: 2 pions with opposite charge
sel_modeII.select_track do
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp "==1"
  nChrn "==1"
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

# PID: pions identified
sel_modeII.pid(method: :probability) do
  prob_cut 0.001
  identify :pip, against: [:kaon, :proton]
  identify :pim, against: [:kaon, :proton]
  npip "==1"
  npim "==1"
end

# Reconstruct η from γγ (1C Kalman fit)
sel_modeII.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 25
  neta ">=1"
end

# 4C kinematic fit: 4γ + π+π- with η as composite
# Participants: γ_radiative + γ(η') + η_composite + π+π-
sel_modeII.kinematic_fit([:gamma, :gamma, :eta, :pip, :pim]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on([data_jpsi, incMC_jpsi, sig_mc_modeII])