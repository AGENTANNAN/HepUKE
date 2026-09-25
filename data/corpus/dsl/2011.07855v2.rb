# ============================================================
# Paper: Observation of a Near-Threshold Structure in the K+ Recoil-Mass Spectra
#        in e+e- → K+(Ds-D*0 + Ds*-D0)
# arXiv: 2011.07855v2
# ============================================================

###
### Dataset preparation: 5 energy points at BOSS 706
###
data_4628 = DatasetManager.real_data.find("706_4620")   # 4.628 GeV, 511.1 pb-1
data_4641 = DatasetManager.real_data.find("706_4640")   # 4.641 GeV, 541.4 pb-1
data_4661 = DatasetManager.real_data.find("706_4660")   # 4.661 GeV, 523.6 pb-1
data_4681 = DatasetManager.real_data.find("706_4680")   # 4.681 GeV, 1643.4 pb-1
data_4698 = DatasetManager.real_data.find("706_4700")   # 4.698 GeV, 526.2 pb-1

data_points = [data_4628, data_4641, data_4661, data_4681, data_4698]

inc_mc_points = data_points.map { |d| DatasetManager.inclusive_mc.find("#{d.boss}_#{d.sample_name}") }

# Decay card for Mode I: K+ D_s-(→K+K-π-) anti-D*0
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.000 K+ D_s- anti-D*0  PHSP;
    Enddecay

    Decay anti-D*0
    1.000 anti-D0 pi0       VSS;
    Enddecay

    Decay anti-D0
    1.000 K+ pi-            PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma        PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi-         PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for Mode II: K+ D_s-(→K_S0 K-) anti-D*0
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.000 K+ D_s- anti-D*0  PHSP;
    Enddecay

    Decay anti-D*0
    1.000 anti-D0 pi0       VSS;
    Enddecay

    Decay anti-D0
    1.000 K+ pi-            PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma        PHSP;
    Enddecay

    Decay D_s-
    1.000 K_S0 K-            PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi-            PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for both modes (multi-energy scan)
exMC_modeI = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_KpDsDst0_KpKmPim"
  config.events        = 500_000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_KpDsDst0_KsKm"
  config.events        = 500_000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

# ============================================================
# Mode I: D_s- → K+ K- π-  (partial reconstruction)
# Reconstruct bachelor K+ and D_s-(→K+K-π-); D*0/D0 is missing
# ============================================================
alg_modeI = Algorithm.new("KpDsDst0_KpKmPim")
alg_modeI.set_header(["KpDsDst0_KpKmPimAlg/KpDsDst0_KpKmPim.h"])

# Multi-energy scan: no ECMS constant (energy injected at runtime)
alg_modeI.set_constant({})
          .set_alias({})

sel_modeI = Selection.new

# --- Charged track selection ---
# Topology: 2 K+ (bachelor + Ds- daughter), 1 K-, 1 π- → 2 positive, 2 negative
sel_modeI.select_track do
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       ">=2"
  nChrn       "==2"
end

# --- Particle identification: kaons ---
# Identify kaons against pion and proton hypotheses
sel_modeI.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp ">=2"
  nkm "==1"
end

# Remove identified kaons from generic charged lists; remaining are pions
sel_modeI.remove([:kp <= :chrgp, :km <= :chrgn])
          .assign({chrgn: :pim})

# --- Partial reconstruction ---
# Tag: K+ bachelor (recID 1), D_s- (recID 2, expands to K+,K-,π-)
# Miss: anti-D*0 (recID 6, auto-expands to anti-D0, pi0, K+, π-, gamma, gamma)
# The recID mapping assumes DFS pre-order traversal of the decay card:
#   0: psi(4260) [skip]
#   1: K+ (bachelor)
#   2: D_s-
#   3: K+ (from D_s-)
#   4: K- (from D_s-)
#   5: pi- (from D_s-)
#   6: anti-D*0 [miss — auto-expands to 7..12]
#   7: anti-D0
#   8: K+ (from anti-D0)
#   9: pi- (from anti-D0)
#  10: pi0
#  11: gamma (from pi0)
#  12: gamma (from pi0)
sel_modeI.partial_miss([6]) do
  best_combination_by_mass :"D_s-", 1.968
  # Signal region: RM(K+Ds-) consistent with D*0 peak
  # Paper uses RM(K+Ds-) + M(Ds-) - m(Ds-) in (1.990, 2.027) GeV/c^2
  # This corresponds approximately to the D*0 mass window in raw RM
  require_recoil_mass 1.990, 2.027
end

alg_modeI.note(:ds_mass_window, "Ds- mass window 1.955 < M(K+K-pi-) < 1.980 GeV/c2 applied; candidate with best Ds- mass selected")
alg_modeI.note(:dalitz_plot_selection, "Ds- candidates retained only in Dalitz regions: M(K+K-) < 1.05 GeV/c2 (phi region) OR 0.850 < M(K+pi-) < 0.930 GeV/c2 (K*(892) region)")
alg_modeI.note(:partial_reconstruction, "Partial reconstruction: only bachelor K+ and Ds-(→K+K-pi-) reconstructed; D*0/D0 inferred from recoil; RM(K+) and RM(K+Ds-) stored for ROOT-level fit")
alg_modeI.note(:wrong_sign_background, "Comb. background modeled with wrong-sign (WS) K-Ds- combinations; WS data used for data-driven background estimation in ROOT")

alg_modeI.with_decay_card(decay_card_modeI)
          .apply(sel_modeI)

# ============================================================
# Mode II: D_s- → K_S0 K-  (partial reconstruction)
# Reconstruct bachelor K+, D_s-(→K_S0 K-); D*0/D0 is missing
# ============================================================
alg_modeII = Algorithm.new("KpDsDst0_KsKm")
alg_modeII.set_header(["KpDsDst0_KsKmAlg/KpDsDst0_KsKm.h"])

# Multi-energy scan: no ECMS constant
alg_modeII.set_constant({})
           .set_alias({})

sel_modeII = Selection.new

# --- Charged track selection ---
# Topology: 1 K+ (bachelor), 1 π+ (from K_S0), 1 π- (from K_S0), 1 K- (from Ds-)
# → 2 positive, 2 negative
sel_modeII.select_track do
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       "==2"
  nChrn       "==2"
end

# --- Particle identification: kaons ---
sel_modeII.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp "==1"
  nkm "==1"
end

# Remove kaons from generic lists; remaining tracks assigned as pions for K_S0 reconstruction
sel_modeII.remove([:kp <= :chrgp, :km <= :chrgn])
           .assign({chrgp: :pip, chrgn: :pim})

# --- K_S0 reconstruction via secondary vertex fit ---
sel_modeII.secondary_vertex_fit([:pip, :pim]) do
  build_virtual_particle(:K_S0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end

# --- Partial reconstruction ---
# Tag: K+ bachelor (recID 1), D_s- (recID 2), K_S0 (recID 3, pre-built), K- (recID 6)
# Miss: anti-D*0 (recID 7, auto-expands to anti-D0, K+, pi-, pi0, gamma, gamma)
# recID mapping (DFS pre-order of the Mode II decay card):
#   0: psi(4260) [skip]
#   1: K+ (bachelor)
#   2: D_s-
#   3: K_S0 (daughter of D_s-, pre-built by secondary_vertex_fit)
#   4: pi+ (from K_S0)
#   5: pi- (from K_S0)
#   6: K- (from D_s-)
#   7: anti-D*0 [miss — auto-expands to 8..13]
#   8: anti-D0
#   9: K+ (from anti-D0)
#  10: pi- (from anti-D0)
#  11: pi0
#  12: gamma (from pi0)
#  13: gamma (from pi0)
sel_modeII.partial_miss([7]) do
  best_combination_by_mass :"D_s-", 1.968
  require_recoil_mass 1.990, 2.027
end

alg_modeII.note(:ds_mass_window, "Ds- mass window 1.955 < M(K_S0 K-) < 1.985 GeV/c2 applied; candidate selected closest to Ds- nominal mass")
alg_modeII.note(:partial_reconstruction, "Partial reconstruction: bachelor K+ and Ds-(→K_S0 K-) reconstructed; D*0/D0 inferred from recoil")
alg_modeII.note(:k_s0_mass_window, "K_S0 mass window 0.485 < M(pi+pi-) < 0.511 GeV/c2 applied via secondary vertex fit")
alg_modeII.note(:wrong_sign_background, "Comb. background modeled with wrong-sign (WS) K-Ds- combinations in ROOT")

alg_modeII.with_decay_card(decay_card_modeII)
           .apply(sel_modeII)

# ============================================================
# Execute both algorithms
# ============================================================
all_real = data_points
all_inc  = inc_mc_points

alg_modeI.execute_on(all_real + all_inc + exMC_modeI)
alg_modeII.execute_on(all_real + all_inc + exMC_modeII)