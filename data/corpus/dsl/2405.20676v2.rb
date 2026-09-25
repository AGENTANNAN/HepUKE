# DSL for paper 2405.20676v2: Search for e+e- → η'ψ(2S) at √s from 4.66 to 4.95 GeV
# ψ(2S) → π+π-J/ψ, J/ψ → e+e-/μ+μ-, η' → γπ+π- (Mode I) / ηπ+π- (Mode II)
# 9 energy points, 4.67 fb-1 total; no significant signal, upper limits at 90% CL

# Data: BOSS 706/707 scan points 4660-4950 MeV
scan_data = DatasetManager.real_data.where(cms_energy: {value: 4660..4950})
scan_incMC = DatasetManager.inclusive_mc.where(cms_energy: {value: 4660..4950})

# ============================================================
# Mode I: η' → γπ+π-, ψ(2S) → π+π-J/ψ, J/ψ → l+l-
# ============================================================
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.000 etap psi(2S) PHSP;
    Enddecay

    Decay etap
    1.000 gamma pi+ pi- PHSP;
    Enddecay

    Decay psi(2S)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHSP;
    Enddecay

    End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name = "sig_etap_psip_modeI"
  config.events = 100000
  config.decay_card = decay_card_modeI
  config.cross_section = :default
end

alg_modeI = Algorithm.new("EtapPsi2SModeI")
alg_modeI.set_header(["EtapPsi2SModeIAlg/EtapPsi2SModeI.h"])
         .set_constant({"ECMS" => [:double, 4.800]})
         .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 4C/1C kinematic fit")
         .note(:lepton_pid, "charged tracks with momentum > 1.0 GeV/c categorized as leptons; electrons if E/p > 0.7, muons if E/p < 0.7 and E < 0.45 GeV; tracks < 0.8 GeV/c as pions; ROOT level")
         .note(:missing_pion_events, "5-track events (net charge ±1) account for ~32% of Mode I; 1C kinematic fit constraining missing particle mass to pion nominal mass applied; ROOT level")
         .note(:jpsi_mass_constraint, "J/ψ invariant mass resolution improved using M(π+π-l+l-) - M(l+l-) + m(J/ψ)")
         .note(:best_candidate, "best candidate selected by minimizing |M(γπ+π-)-m(η')|²/σ²_η' + |M(π+π-J/ψ)-m(ψ(2S))|²/σ²_ψ(2S); ROOT level")
         .note(:competing_hypothesis_veto, "χ²_4C,nγ < χ²_4C,(n-1)γ required for Mode I to suppress background with unexpected photon count; ROOT level")
         .note(:jpsi_mass_window, "J/ψ mass window [3.083, 3.111] GeV/c² for Mode I; ROOT level")
         .note(:psi2s_mass_window, "ψ(2S) mass window [3.680, 3.693] GeV/c²; ROOT level")
         .note(:cross_section_formula, "Born cross section σ_B = N_sig / (L_int × (1+δ) × (1/|1-Π|²) × ε × B_intermediate); upper limits at 90% CL via profile likelihood; simultaneous fit to M(γπ+π-) and M(ηπ+π-) sharing Born cross section")

sel_modeI = Selection.new

sel_modeI.select_track do
  cos_theta  0.93
  Vz         10.0
  Vr         1.0
  nChrp      ">=2"
  nChrn      ">=2"
end

sel_modeI.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=1"
end

sel_modeI.pid(method: :probability) do
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                  treat_as_electron_if_energy_over_momentum_above: 0.7
  nlp ">=1"
  nlm ">=1"
end

sel_modeI.remove([:lp <= :chrgp, :lm <= :chrgn])
         .assign({chrgp: :pip, chrgn: :pim})

# 4C kinematic fit: 6-track events (γ + π+π-π+π- + l+l-)
sel_modeI.kinematic_fit([:gamma, :pip, :pip, :pim, :pim, :lp, :lm]) do
  nominal
  constrain_four_momentum
  chi2_cut 200  # paper cut < 30 for Mode I; optimal in ROOT
end

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)

# ============================================================
# Mode II: η' → ηπ+π-, η → γγ, ψ(2S) → π+π-J/ψ, J/ψ → l+l-
# ============================================================
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.000 etap psi(2S) PHSP;
    Enddecay

    Decay etap
    1.000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay psi(2S)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHSP;
    Enddecay

    End
DECAYCARD

exMC_modeII = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name = "sig_etap_psip_modeII"
  config.events = 100000
  config.decay_card = decay_card_modeII
  config.cross_section = :default
end

alg_modeII = Algorithm.new("EtapPsi2SModeII")
alg_modeII.set_header(["EtapPsi2SModeIIAlg/EtapPsi2SModeII.h"])
          .set_constant({"ECMS" => [:double, 4.800]})
          .note(:helix_correction, "helix parameter correction applied to all charged tracks before the 4C/1C kinematic fit")
          .note(:lepton_pid, "charged tracks with momentum > 1.0 GeV/c categorized as leptons; electrons if E/p > 0.7, muons if E/p < 0.7 and E < 0.45 GeV; tracks < 0.8 GeV/c as pions; ROOT level")
          .note(:missing_pion_events, "5-track events (net charge ±1) account for ~37% of Mode II; 1C kinematic fit constraining missing particle mass to pion nominal mass applied; ROOT level")
          .note(:jpsi_mass_constraint, "J/ψ invariant mass resolution improved using M(π+π-l+l-) - M(l+l-) + m(J/ψ)")
          .note(:best_candidate, "best candidate selected by minimizing |M(ηπ+π-)-m(η')|²/σ²_η' + |M(π+π-J/ψ)-m(ψ(2S))|²/σ²_ψ(2S); ROOT level")
          .note(:jpsi_mass_window, "J/ψ mass window [3.073, 3.121] GeV/c² for Mode II; ROOT level")
          .note(:eta_mass_window, "η mass window [0.482, 0.604] GeV/c² for Mode II; ROOT level")
          .note(:psi2s_mass_window, "ψ(2S) mass window [3.680, 3.693] GeV/c²; ROOT level")
          .note(:cross_section_formula, "Born cross section σ_B = N_sig / (L_int × (1+δ) × (1/|1-Π|²) × ε × B_intermediate); simultaneous fit shares Born cross section with Mode I")

sel_modeII = Selection.new

sel_modeII.select_track do
  cos_theta  0.93
  Vz         10.0
  Vr         1.0
  nChrp      ">=2"
  nChrn      ">=2"
end

sel_modeII.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=2"
end

sel_modeII.pid(method: :probability) do
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                  treat_as_electron_if_energy_over_momentum_above: 0.7
  nlp ">=1"
  nlm ">=1"
end

sel_modeII.remove([:lp <= :chrgp, :lm <= :chrgn])
          .assign({chrgp: :pip, chrgn: :pim})

# Reconstruct η → γγ via 1C Kalman fit
sel_modeII.kalman_kinematic_fit([:gamma, :gamma]) do
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 200
  neta ">=1"
end

# 4C kinematic fit: 6-track events (η + π+π-π+π- + l+l-)
sel_modeII.kinematic_fit([:eta, :pip, :pip, :pim, :pim, :lp, :lm]) do
  nominal
  constrain_four_momentum
  chi2_cut 200  # paper cut < 60 for Mode II; optimal in ROOT
end

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

# ============================================================
# Execute all algorithms
# ============================================================
alg_modeI.execute_on([scan_data, scan_incMC, exMC_modeI])
alg_modeII.execute_on([scan_data, scan_incMC, exMC_modeII])