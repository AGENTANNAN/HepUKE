# BESIII DSL for arXiv:2403.13437v2
# Search for ΔS=2 nonleptonic hyperon decays Ω- → Σ0π- and Ω- → nK-
# ψ(3686) → Ω- Ω+, Ω+ → Λ K+ (ST), Ω- signal side
# Double-tag method with Ω hyperons

# ============================================================
# Datasets
# ============================================================
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# ============================================================
# Decay card for ψ(3686) → Ω- Ω+ (shared)
# ============================================================
decay_card_omega_pair = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Omega- anti-Omega+ PHSP;
    Enddecay

    End
DECAYCARD

# ============================================================
# Signal 1: Ω- → Σ0 π-
# Ω+ → Λ K+, Λ → pbar π+, Σ0 → X (inclusive)
# Only bachelor π- reconstructed on signal side
# ============================================================
decay_card_omega_sigma_pi = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Omega- anti-Omega+ PHSP;
    Enddecay

    Decay anti-Omega+
    1.0000 Lambda K+ PHSP;
    Enddecay

    Decay Omega-
    1.0000 Sigma0 pi- PHSP;
    Enddecay

    Decay Lambda
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda PHSP;
    Enddecay

    End
DECAYCARD

exMC_sigma_pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "omega_to_sigma0_pi"
  config.related_dataset = psip_data
  config.events = 1_270_000
  config.decay_card = decay_card_omega_sigma_pi
  config.cross_section = :default
end

# Algorithm 1: Ω- → Σ0 π-
# ST: Ω+ → Λ(→ pbar π+) K+
# Signal: bachelor π- only, Σ0 → X inclusive (not reconstructed)
alg_omega_sigma_pi = Algorithm.new("OmegaToSigma0Pi")
alg_omega_sigma_pi.set_header(["OmegaToSigma0PiAlg/OmegaToSigma0Pi.h"])
                   .set_constant({"ECMS" => [:double, 3.686]})

sel_omega_sigma_pi = Selection.new

# Charged track selection: |cosθ| < 0.93, Vz 10cm, Vr 1cm
# Tag side: pbar, π+ (from Λ), K+; Signal side: π-
# Total on tag side: 3 tracks minimum
sel_omega_sigma_pi.select_track do
  cos_theta  0.93
  Vz         10.0
  Vr         1.0
  nChrp      ">=2"    # K+, π+ (from Λ)
  nChrn      ">=2"    # pbar (from Λ) + π- (signal side bachelor)
  nNet       "==0"
end

# PID: proton (pbar), kaon (K+), pion (π+ from Λ, π- bachelor)
sel_omega_sigma_pi.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprm ">=1"     # anti-proton from Lambda decay
end

sel_omega_sigma_pi.remove([:prm <= :chrgn])

sel_omega_sigma_pi.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp ">=1"       # K+ from Omega+ decay
end

sel_omega_sigma_pi.remove([:kp <= :chrgp])
               .assign({chrgp: :pip, chrgn: :pim})

# Reconstruct Lambda → pbar π+
sel_omega_sigma_pi.secondary_vertex_fit([:prm, :pip]) do
  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end

# Reconstruct ST Omega+: Lambda_bar K+ vertex fit + Omega mass
# The Omega+ signal region: M_Lambda_K+ in [1.664, 1.680] and RM_Omega+ in [1.652, 1.695]
# Note: these mass windows are applied at ROOT level; only the kinematic fit chain is expressed here

# ST Omega+ candidate selection: vertex fit Lambda K+
# Partial reconstruction: only the bachelor π- from Ω- decay is reconstructed
# The Σ0 → X is NOT reconstructed (inclusive tag on Ω- via recoil mass)
# RMS_Omega+_pi- distribution used to extract DT yield at ROOT level

sel_omega_sigma_pi.kinematic_fit([:Lambda_bar, :kp, :pim]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

alg_omega_sigma_pi.with_decay_card(decay_card_omega_sigma_pi)
                  .note(:omega_st_tag, "ST Omega+ selection: vertex fit Lambda_bar K+, Omega mass window M_Lambda_K+ in [1.664, 1.680] GeV/c2; ST yield from fit to RM_Omega+ spectrum; ROOT level")
                  .note(:partial_signal_side, "signal side: only bachelor pi- reconstructed (highest momentum pi-); Sigma0 → X is inclusive (not reconstructed); Omega- → Sigma0 pi- yield extracted from RM_Omega+pi- distribution; ROOT level")
                  .note(:dt_yield_extraction, "DT yield from fit to RM_Omega+pi- distribution; signal shape from MC, background from Chebychev polynomial")
                  .apply(sel_omega_sigma_pi)

# ============================================================
# Signal 2: Ω- → n K-
# Ω+ → Λ K+ (ST), signal side: bachelor K-
# ============================================================
decay_card_omega_nK = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Omega- anti-Omega+ PHSP;
    Enddecay

    Decay anti-Omega+
    1.0000 Lambda K+ PHSP;
    Enddecay

    Decay Omega-
    1.0000 anti-n- K- PHSP;
    Enddecay

    Decay Lambda
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

exMC_nK = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "omega_to_nK"
  config.related_dataset = psip_data
  config.events = 1_270_000
  config.decay_card = decay_card_omega_nK
  config.cross_section = :default
end

alg_omega_nK = Algorithm.new("OmegaToNK")
alg_omega_nK.set_header(["OmegaToNKAlg/OmegaToNK.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

sel_omega_nK = Selection.new

# Charged track selection
# Tag: pbar, π+ (Lambda), K+; Signal: K-
# Total 4 charged tracks required (paper: N_tracks = 4 on both sides)
sel_omega_nK.select_track do
  cos_theta  0.93
  Vz         10.0
  Vr         1.0
  nChrp      ">=2"     # K+ (tag), π+ (Lambda)
  nChrn      ">=2"     # pbar (Lambda), K- (signal)
  nTot       "==4"     # exactly 4 charged tracks total
  nNet       "==0"
end

# PID: proton (pbar), kaon (K+, K-)
sel_omega_nK.pid(method: :probability) do
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  nprm ">=1"
end

sel_omega_nK.remove([:prm <= :chrgn])

sel_omega_nK.pid(method: :probability) do
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp ">=1"      # K+ from Omega+
  nkm ">=1"      # K- from Omega-
end

sel_omega_nK.remove([:kp <= :chrgp, :km <= :chrgn])
            .assign({chrgp: :pip, chrgn: :pim})

# Reconstruct Lambda → pbar π+
sel_omega_nK.secondary_vertex_fit([:prm, :pip]) do
  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
end

# ST Omega+ with Lambda_bar K+; signal K- from Omega- → n K-
# Neutron is not reconstructed
sel_omega_nK.kinematic_fit([:Lambda_bar, :kp, :km]) do
  nominal
  constrain_four_momentum
  chi2_cut 200
end

alg_omega_nK.with_decay_card(decay_card_omega_nK)
            .note(:omega_st_tag, "ST Omega+ selection: vertex fit Lambda_bar K+, Omega mass window M_Lambda_K+ in [1.664, 1.680] GeV/c2; ST yield from fit to RM_Omega+; ROOT level")
            .note(:partial_signal_side, "signal side: only bachelor K- reconstructed; neutron not reconstructed; DT yield from fit to RM_Omega+K-; ROOT level")
            .note(:four_track_requirement, "N_tracks = 4 (tag + signal) required on both sides to suppress background")
            .apply(sel_omega_nK)

# ============================================================
# Execute
# ============================================================
alg_omega_sigma_pi.execute_on([psip_data, psip_incMC, exMC_sigma_pi])
alg_omega_nK.execute_on([psip_data, psip_incMC, exMC_nK])