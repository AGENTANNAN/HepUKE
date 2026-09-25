# Paper: 2212.03693v2 — Study of e+e- → Omega- Omega+ at 3.49-3.67 GeV
# Single hyperon tag technique:
#   Reconstruct Omega- → Lambda K-, Lambda → p pi-; Omega+ inferred via recoil mass
#   Also c.c. mode: Omega+ reconstructed, Omega- missed (separate Algorithm, Rule T1)
# 8 energy points near chi_c1 threshold

### Dataset preparation ###
scan_data = [
  DatasetManager.real_data.find("703_chic1_scan_1"),   # 3.4900
  DatasetManager.real_data.find("703_chic1_scan_3"),   # 3.5080
  DatasetManager.real_data.find("703_chic1_scan_2"),   # 3.5097
  DatasetManager.real_data.find("703_chic1_scan_4"),   # 3.5104
  DatasetManager.real_data.find("703_chic1_scan_5"),   # 3.5146
  DatasetManager.real_data.find("704_psip_scan_1"),    # 3.5815
  DatasetManager.real_data.find("709_3650"),            # 3.6500
  DatasetManager.real_data.find("704_psip_scan_2"),    # 3.6702
]

scan_incMC = scan_data.map do |ds|
  DatasetManager.inclusive_mc.find(ds.sample_name)
end

# Decay card: e+e- → Omega- Omega+
decay_card = <<~DECAYCARD
    Decay vpho
    1.0 Omega- anti-Omega+ PHSP;
    Enddecay

    Decay Omega-
    1.0 Lambda K- PHSP;
    Enddecay

    Decay Lambda
    1.0 p+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name = "sig_Omega_Omega"
  config.events = 200000
  config.decay_card = decay_card
  config.cross_section = :default
end

# ============================================================
# Algorithm I: Reconstruct Omega-, miss Omega+ via recoil mass
# DecayCard recIDs:
#   0: vpho        [skip]
#   1: Omega-      [reconstruct via SV fit]
#   2: anti-Omega+ [miss — inferred from recoil mass]
#   3: Lambda      [reconstruct via SV fit]
#   4: K-          [from Omega- — reconstruct]
#   5: p+          [from Lambda — reconstruct]
#   6: pi-         [from Lambda — reconstruct]
# Reconstruct recIDs: 1 (Omega-, expands to 3,4,5,6)
# ============================================================
alg_Omega_tag = Algorithm.new("OmegaMinusTag")
alg_Omega_tag.set_header(["OmegaMinusTagAlg/OmegaMinusTag.h"])
               .set_constant({ "ECMS" => [:double, 3.510] })
               .note(:single_hyperon_tag,
                 "Single hyperon tag: reconstruct Omega- → Lambda K-, Lambda → p pi-.
                  Secondary vertex fits for Lambda and Omega-.
                  Best Omega- candidate by δ = (M_pπ - m_Λ)² + (M_KΛ - m_Ω)² minimization.
                  Omega+ inferred from RM(K-Lambda) corrected: RM_corr = RM + M_KΛ - m_Ω.
                  Lambda mass window: |M_pπ - m_Λ| < 4 MeV. Both decay lengths > 0.
                  Signal region: M_ΛK in [1.6665, 1.6785], RM_ΛK in [1.65, 1.69] GeV/c².")
               .note(:pid_highest_confidence_level,
                 "Tracks identified as pi/K/p by highest PID probability from dE/dx + TOF.
                  Event requires ≥1 p+, ≥1 pi-, ≥1 K- after PID assignment.
                  Cannot be directly expressed as standard PID blocks — requires custom logic.")

sel_Omega_tag = Selection.new
sel_Omega_tag
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=1"
    nChrn ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npim ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkm ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
    nprp ">=1"
  }
  # Lambda secondary vertex fit: p+ pi- → Lambda
  .secondary_vertex_fit(:Lambda, daughters: [:prp, :pim]) {
    require_mass_window :Lambda, 0.004
    decay_length_significance_cut 0.0
  }
  # Omega- secondary vertex fit: Lambda K-
  .secondary_vertex_fit(:Omega_minus, daughters: [:Lambda, :km]) {
    decay_length_significance_cut 0.0
  }
  # Partial reconstruction: reconstruct Omega-, miss Omega+
  .partial_rec([1, 3, 4, 5, 6]) {
    best_combination_by_mass :Omega, 1.67245
    require_recoil_mass 1.59, 1.75
  }

alg_Omega_tag.with_decay_card(decay_card).apply(sel_Omega_tag)

# ============================================================
# Algorithm II: Reconstruct Omega+, miss Omega- via recoil mass
# (c.c. mode — identical logic with charge conjugation, separate Algorithm per Rule T1)
# ============================================================
alg_OmegaPlus_tag = Algorithm.new("OmegaPlusTag")
alg_OmegaPlus_tag.set_header(["OmegaPlusTagAlg/OmegaPlusTag.h"])
                   .set_constant({ "ECMS" => [:double, 3.510] })
                   .note(:single_hyperon_tag,
                     "c.c. of Omega- tag: reconstruct Omega+ → anti-Lambda K+,
                      anti-Lambda → pbar pi+; Omega- via recoil mass.")

sel_OmegaPlus_tag = Selection.new
sel_OmegaPlus_tag
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:pion, :kaon]
    nprm ">=1"
  }
  .secondary_vertex_fit(:Lambda_bar, daughters: [:prm, :pip]) {
    require_mass_window :Lambda_bar, 0.004
    decay_length_significance_cut 0.0
  }
  .secondary_vertex_fit(:Omega_plus, daughters: [:Lambda_bar, :kp]) {
    decay_length_significance_cut 0.0
  }
  .partial_rec([1, 3, 4, 5, 6]) {
    best_combination_by_mass :Omega, 1.67245
    require_recoil_mass 1.59, 1.75
  }

alg_OmegaPlus_tag.with_decay_card(decay_card).apply(sel_OmegaPlus_tag)

# ============================================================
# Execute all algorithms on the scan datasets
# ============================================================
all_datasets = scan_data + scan_incMC + exMC_signal
alg_Omega_tag.execute_on(all_datasets)
alg_OmegaPlus_tag.execute_on(all_datasets)