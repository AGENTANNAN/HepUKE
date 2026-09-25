# =============================================================================
# BOSS DSL — ψ(3770) → D⁰ D̄⁰ : CP-eigenstate tag + semileptonic signal
# (y_CP measurement from double tags, K⁻e⁺ν_e and K⁻μ⁺ν_μ)
# =============================================================================

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("712_3773")     # 2.92 fb⁻¹ real data at 3.773 GeV
psip_incMC = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

# -----------------------------------------------------------------------------
# Decay cards: ψ(3770) → D⁰ D̄⁰ ; the D⁰ decays semileptonically (signal),
# the D̄⁰ decays to the CP eigenstate (tag). Ten combinations:
# 5 CP modes { K⁺K⁻, π⁺π⁻, K_S⁰π⁰π⁰, K_S⁰π⁰, K_S⁰ω } × { e, μ }.
# -----------------------------------------------------------------------------

# --- CP = K⁺K⁻ , signal = K⁻e⁺ν_e ---
dc_KK_e = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.000 K- e+ nu_e PHSP;
  Enddecay

  Decay anti-D0
  1.000 K+ K- PHSP;
  Enddecay

  End
DECAYCARD

# --- CP = K⁺K⁻ , signal = K⁻μ⁺ν_μ ---
dc_KK_mu = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.000 K- mu+ nu_mu PHSP;
  Enddecay

  Decay anti-D0
  1.000 K+ K- PHSP;
  Enddecay

  End
DECAYCARD

# --- CP = π⁺π⁻ , signal = K⁻e⁺ν_e ---
dc_PiPi_e = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.000 K- e+ nu_e PHSP;
  Enddecay

  Decay anti-D0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# --- CP = π⁺π⁻ , signal = K⁻μ⁺ν_μ ---
dc_PiPi_mu = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.000 K- mu+ nu_mu PHSP;
  Enddecay

  Decay anti-D0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# --- CP = K_S⁰π⁰π⁰ , signal = K⁻e⁺ν_e ---
dc_KsPi0Pi0_e = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.000 K- e+ nu_e PHSP;
  Enddecay

  Decay anti-D0
  1.000 K_S0 pi0 pi0 PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- CP = K_S⁰π⁰π⁰ , signal = K⁻μ⁺ν_μ ---
dc_KsPi0Pi0_mu = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.000 K- mu+ nu_mu PHSP;
  Enddecay

  Decay anti-D0
  1.000 K_S0 pi0 pi0 PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- CP = K_S⁰π⁰ , signal = K⁻e⁺ν_e ---
dc_KsPi0_e = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.000 K- e+ nu_e PHSP;
  Enddecay

  Decay anti-D0
  1.000 K_S0 pi0 PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- CP = K_S⁰π⁰ , signal = K⁻μ⁺ν_μ ---
dc_KsPi0_mu = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.000 K- mu+ nu_mu PHSP;
  Enddecay

  Decay anti-D0
  1.000 K_S0 pi0 PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- CP = K_S⁰ω , signal = K⁻e⁺ν_e ---
dc_KsOmega_e = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.000 K- e+ nu_e PHSP;
  Enddecay

  Decay anti-D0
  1.000 K_S0 omega PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- CP = K_S⁰ω , signal = K⁻μ⁺ν_μ ---
dc_KsOmega_mu = <<~DECAYCARD
  Decay psi(3770)
  1.000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.000 K- mu+ nu_mu PHSP;
  Enddecay

  Decay anti-D0
  1.000 K_S0 omega PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# -----------------------------------------------------------------------------
# Exclusive MC: 200k events for each of the ten CP-tag × semileptonic combos.
# -----------------------------------------------------------------------------
exMC_KK_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0CP_KK_enu"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = dc_KK_e
  config.cross_section   = :default
end

exMC_KK_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0CP_KK_munu"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = dc_KK_mu
  config.cross_section   = :default
end

exMC_PiPi_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0CP_PiPi_enu"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = dc_PiPi_e
  config.cross_section   = :default
end

exMC_PiPi_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0CP_PiPi_munu"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = dc_PiPi_mu
  config.cross_section   = :default
end

exMC_KsPi0Pi0_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0CP_KsPi0Pi0_enu"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = dc_KsPi0Pi0_e
  config.cross_section   = :default
end

exMC_KsPi0Pi0_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0CP_KsPi0Pi0_munu"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = dc_KsPi0Pi0_mu
  config.cross_section   = :default
end

exMC_KsPi0_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0CP_KsPi0_enu"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = dc_KsPi0_e
  config.cross_section   = :default
end

exMC_KsPi0_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0CP_KsPi0_munu"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = dc_KsPi0_mu
  config.cross_section   = :default
end

exMC_KsOmega_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0CP_KsOmega_enu"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = dc_KsOmega_e
  config.cross_section   = :default
end

exMC_KsOmega_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0CP_KsOmega_munu"
  config.related_dataset = psip_data
  config.events          = 200000
  config.decay_card      = dc_KsOmega_mu
  config.cross_section   = :default
end

### Tag analysis — electron channel (K⁻e⁺ν_e) ###
alg_e = TagAnalysis.new("D0CPTagKenu")
alg_e.set_header(["D0CPTagKenuAlg/D0CPTagKenu.h"])
     .set_constant({"ECMS" => [:double, 3.773]})
     .note(:tag_track_selection, "tag tracks: |cosθ|<0.93, |Vz|<10 cm, Vr<1 cm; "
                                 "hadron PID from combined MDC dE/dx and TOF likelihoods "
                                 "(K: L(K)>L(π), π: otherwise); K_S0 daughters exempt from PID "
                                 "and use a looser 20 cm longitudinal IP cut")
     .note(:deltaE_window, "mode-dependent ΔE windows applied in ROOT (stored, not cut at BOSS): "
                           "K+K− ±0.020, π+π− ±0.030, K_S0π0π0 [−0.080,0.045], "
                           "K_S0π0 [−0.070,0.040], K_S0ω [−0.050,0.030] GeV; "
                           "only the smallest-|ΔE| candidate is retained")
     .note(:photon_selection, "photons: EMC showers separated by >10σ from any charged track, "
                              "E>25 MeV in the barrel (|cosθ|<0.80) and E>50 MeV in the endcap "
                              "(0.84<|cosθ|<0.92)")
     .note(:pi0_eta_reconstruction, "π0 and η built from photon pairs with 1C mass constraints: "
                                    "π0 in 0.115–0.150 GeV/c² and η in 0.505–0.570 GeV/c²")
     .note(:omega_selection, "ω → π+π−π0 invariant mass window 0.7600–0.8050 GeV/c² "
                             "with sideband subtraction")
     .note(:background_veto, "cosmic/Bhabha veto applied to the K+K− and π+π− tags")
     .note(:pid_correction_method, "signal-side electron: L'(e)>0.001 and R'(e)>0.8")

alg_e.tag_side(:D0) do |t|
  # CP-eigenstate tag modes (both D⁰ / D̄⁰ scanned)
  t.modes :D0toKK, :D0toPiPi, :D0toKsPi0Pi0, :D0toKsPi0, :D0toKsOmega
  # explicit opt-in tag-side window (ΔE windows are mode-dependent → handled in ROOT)
  t.window :mBC, min: 1.855, max: 1.875
end

alg_e.signal_side do |s|
  # exactly two oppositely charged unused tracks (K⁻ , e⁺), net charge 0, massless ν_e
  s.charged(km: 1, ep: 1)
  s.require_charge 0
  s.missing :nu_e
end

alg_e.fit do |f|
  f.constrain_four_momentum                                  # 4C energy–momentum constraint
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)  # tag D⁰ mass to nominal
  f.chi2_cut 200
end

alg_e.apply
alg_e.execute_on([psip_data, psip_incMC,
                  exMC_KK_e, exMC_PiPi_e, exMC_KsPi0Pi0_e, exMC_KsPi0_e, exMC_KsOmega_e])

### Tag analysis — muon channel (K⁻μ⁺ν_μ) ###
alg_mu = TagAnalysis.new("D0CPTagKmunu")
alg_mu.set_header(["D0CPTagKmunuAlg/D0CPTagKmunu.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .note(:tag_track_selection, "tag tracks: |cosθ|<0.93, |Vz|<10 cm, Vr<1 cm; "
                                  "hadron PID from combined MDC dE/dx and TOF likelihoods "
                                  "(K: L(K)>L(π), π: otherwise); K_S0 daughters exempt from PID "
                                  "and use a looser 20 cm longitudinal IP cut")
      .note(:deltaE_window, "mode-dependent ΔE windows applied in ROOT (stored, not cut at BOSS): "
                            "K+K− ±0.020, π+π− ±0.030, K_S0π0π0 [−0.080,0.045], "
                            "K_S0π0 [−0.070,0.040], K_S0ω [−0.050,0.030] GeV; "
                            "only the smallest-|ΔE| candidate is retained")
      .note(:photon_selection, "photons: EMC showers separated by >10σ from any charged track, "
                               "E>25 MeV in the barrel (|cosθ|<0.80) and E>50 MeV in the endcap "
                               "(0.84<|cosθ|<0.92)")
      .note(:pi0_eta_reconstruction, "π0 and η built from photon pairs with 1C mass constraints: "
                                     "π0 in 0.115–0.150 GeV/c² and η in 0.505–0.570 GeV/c²")
      .note(:omega_selection, "ω → π+π−π0 invariant mass window 0.7600–0.8050 GeV/c² "
                              "with sideband subtraction")
      .note(:background_veto, "cosmic/Bhabha veto applied to the K+K− and π+π− tags")
      .note(:pid_correction_method, "signal-side muon: EMC deposit <0.3 GeV and R_l'(e)<0.8")
      .note(:kmu_selection, "Kμ channel: M(Kμ)<1.65 GeV/c² and unmatched EMC energy <0.2 GeV")

alg_mu.tag_side(:D0) do |t|
  t.modes :D0toKK, :D0toPiPi, :D0toKsPi0Pi0, :D0toKsPi0, :D0toKsOmega
  t.window :mBC, min: 1.855, max: 1.875
end

alg_mu.signal_side do |s|
  # exactly two oppositely charged unused tracks (K⁻ , μ⁺), net charge 0, massless ν_μ
  s.charged(km: 1, mup: 1)
  s.require_charge 0
  s.missing :nu_mu
end

alg_mu.fit do |f|
  f.constrain_four_momentum                                     # 4C energy–momentum constraint
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)  # tag D⁰ mass to nominal
  f.chi2_cut 200
end

alg_mu.apply
alg_mu.execute_on([psip_data, psip_incMC,
                   exMC_KK_mu, exMC_PiPi_mu, exMC_KsPi0Pi0_mu, exMC_KsPi0_mu, exMC_KsOmega_mu])