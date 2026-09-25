# 2112.14369v1: Observation of EM Dalitz decays J/ψ → e+e-π+π-η'
# with η' → γπ+π- (Mode I) and η' → π+π-η, η → γγ (Mode II).
# √s = 3.097 GeV, ~10 billion J/ψ events.
# Two intermediate states: X(1835), X(2120), X(2370) → π+π-η'.
# Transition form factor (TFF) measurement for J/ψ → e+e-X(1835).
# Two independent decay modes → two Algorithm objects (Rule T1).

# ── Datasets ──────────────────────────────────────────────────────
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card: J/ψ → e+e-π+π-η' (the full EM Dalitz decay chain)
# The η' is decayed generically; specific modes selected in reconstruction.
decay_card = <<~DECAYCARD
    Decay J/psi
    1.000 e+ e- pi+ pi- etaprime PHSP;
    Enddecay
    Decay etaprime
    0.289 pi+ pi- gamma PHSP;
    0.411 pi+ pi- eta PHSP;
    Enddecay
    Decay eta
    0.394 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_ee_pipi_etap"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# ═══════════════════════════════════════════════════════════════════
# Mode I: η' → γπ+π-
# 6 charged tracks (e+e- + 4 pions) + 1 photon → 4C kinematic fit
# ═══════════════════════════════════════════════════════════════════
alg_I = Algorithm.new("JPsiEEPiPiEtap_ModeI")
alg_I.set_header(["JPsiEEPiPiEtapAlg/JPsiEEPiPiEtap.h"])
      .set_constant({ "ECMS" => [:double, 3.097] })
      .with_decay_card(decay_card)

sel_I = Selection.new
sel_I.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nTot  "==6"        # exactly 6 charged tracks
  nNet  "==0"        # net charge zero
}
.select_photon {
  energyThreshold_b 0.025      # barrel (|cosθ| < 0.80): 25 MeV
  energyThreshold_e 0.050      # endcap (0.86 < |cosθ| < 0.92): 50 MeV
  nGam ">=1"                   # at least 1 photon from η' → γπ+π-
}
.pid(method: :probability) {
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.6
  identify :pion, against: [:kaon, :proton]
  nlp  "==1"          # exactly 1 e+ from lepton list
  nlm  "==1"          # exactly 1 e- from lepton list
  npip "==2"          # 2 π+ (from X decay + η' decay)
  npim "==2"          # 2 π-
}
# 4C kinematic fit: e+e- + 4π + γ
.kinematic_fit([:lp, :lm, :pip, :pim, :pip, :pim, :gamma]) {
  constrain_four_momentum
  chi2_cut 200
  nominal
}

alg_I
  .note(:electron_pid, "e+ PID: combined dE/dx + TOF + EMC → probabilities Prob(i)_{i=e;π;K}. Track is electron when Prob(e) > Prob(π) and Prob(K). Remaining tracks identified as pions. E/p > 0.85c equivalent via high-momentum lepton path.")
  .note(:photon_angle, "Photons must be separated from charged tracks by more than 10° (min angle).")
  .note(:etaprime_modeI, "η' → γπ+π- mode. η' candidates selected with |M(γπ+π-) − m_η'| < 15 MeV/c² (ROOT-level cut). ρ0 mass window: 575 < M_π+π- < 920 MeV/c² for the π+π- from η' decay (ROOT-level cut).")
  .note(:gamma_conversion, "Reject γ conversion background (J/ψ→γπ+π-η') via R_xy < 2 cm based on γ conversion finder algorithm (ROOT-level cut).")
  .note(:ee_veto, "Remove events with |M_e+e- − m_ω| < 25 MeV/c² and |M_e+e- − m_φ| < 30 MeV/c² to suppress J/ψ→ωπ+π-η' and J/ψ→φπ+π-η' backgrounds (ROOT-level cuts).")
  .note(:gamma_ee_veto, "Remove events with M_γe+e- < 210 MeV/c², 490 < M_γe+e- < 600 MeV/c², and 700 < M_γe+e- < 820 MeV/c² to suppress π0, η, ω → γe+e- Dalitz backgrounds (ROOT-level cuts).")
  .note(:chi2_cut_actual, "Paper applies χ²_4C < 60 at BOSS level (optimized for max S/√(S+B)). Per Rule T3, DSL uses loose chi2_cut 200; tight cut applied in ROOT.")
  .note(:signal_extraction, "Simultaneous unbinned ML fit to M_π+π-η' spectra from both η' decay modes (1.36–2.80 GeV/c²). Six resonances: f1(1510), X(1835)+X(1870) coherent sum (Flatté-like line shape distortion near ppbar threshold), X(2120), X(2370), X(2600). Non-η' background from η' mass sideband; J/ψ→π0π+π-η' from PHSP MC.")
  .note(:tff_measurement, "Transition form factor |F(q²)|² measured by comparing M_e+e- spectrum in 5 intervals with point-like QED prediction. Simple pole approximation: F(q²) = 1/(1 − q²/Λ²), Λ = 1.75±0.29±0.05 GeV/c².")
  .apply(sel_I)

# ═══════════════════════════════════════════════════════════════════
# Mode II: η' → π+π-η, η → γγ
# 6 charged tracks (e+e- + 4 pions) + 2 photons → 5C kinematic fit
# ═══════════════════════════════════════════════════════════════════
alg_II = Algorithm.new("JPsiEEPiPiEtap_ModeII")
alg_II.set_header(["JPsiEEPiPiEtapAlg/JPsiEEPiPiEtap.h"])
       .set_constant({ "ECMS" => [:double, 3.097] })
       .with_decay_card(decay_card)

sel_II = Selection.new
sel_II.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nTot  "==6"
  nNet  "==0"
}
.select_photon {
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  nGam ">=2"          # at least 2 photons for η → γγ
}
.pid(method: :probability) {
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                 treat_as_electron_if_energy_above: 0.6
  identify :pion, against: [:kaon, :proton]
  nlp  "==1"
  nlm  "==1"
  npip "==2"
  npim "==2"
}
# Reconstruct η → γγ via Kalman fit
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 25
  neta ">=1"
}
# 5C kinematic fit: e+e- + 4π + η (4C + η mass constraint)
.kinematic_fit([:lp, :lm, :pip, :pim, :pip, :pim, :eta]) {
  constrain_four_momentum
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 200
  nominal
}

alg_II
  .note(:electron_pid, "e+ PID: combined dE/dx + TOF + EMC → probabilities. Electron: Prob(e) > Prob(π) and Prob(K). Remaining = pions.")
  .note(:photon_angle, "Photons separated from charged tracks by > 10°.")
  .note(:etaprime_modeII, "η' → π+π-η, η → γγ mode. η' candidates: |M(π+π-η) − m_η'| < 8.1 MeV/c² (ROOT-level cut). All η' combinations in an event kept.")
  .note(:gamma_conversion, "γ conversion rejection: R_xy < 2 cm (ROOT-level cut).")
  .note(:ee_veto, "Mass veto: |M_e+e- − m_ω| < 25 MeV/c² and |M_e+e- − m_φ| < 30 MeV/c² (ROOT-level cuts).")
  .note(:chi2_cut_actual, "Paper applies χ²_5C < 60 (optimized). Per Rule T3, DSL uses loose chi2_cut 200.")
  .note(:signal_extraction, "Same simultaneous unbinned ML fit as Mode I. Two solutions (constructive/destructive interference between X(1835) and X(1870)). Branching fractions reported for both solutions.")
  .apply(sel_II)

# Execute both algorithms
alg_I.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
alg_II.execute_on([jpsi_data, jpsi_incMC, exMC_signal])