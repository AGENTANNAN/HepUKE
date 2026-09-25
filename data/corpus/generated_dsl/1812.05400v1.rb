# =============================================================================
# ψ(3770) single-tag D⁻ recoil analysis:  D⁺ → K_{S,L}⁰ K⁺ (π⁰)  @ √s = 3.773 GeV
# Tag side: D⁻ reconstructed from pre-stored DTag candidates (six hadronic modes)
# Signal side: recoil D⁺  (four independent signal modes → four TagAnalysis specs)
# =============================================================================

### ---------------------------------- Dataset preparation ---------------------------------- ###
psip_data  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data
psip_incMC = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

# --- Decay cards, one per signal mode (EvtGen format) ---
decay_card_ks_k = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-
  Enddecay
  Decay D+
  1.0000 K_S0 K+
  Enddecay
  Decay K_S0
  1.0000 pi+ pi-
  Enddecay
  End
DECAYCARD

decay_card_ks_k_pi0 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-
  Enddecay
  Decay D+
  1.0000 K_S0 K+ pi0
  Enddecay
  Decay K_S0
  1.0000 pi+ pi-
  Enddecay
  Decay pi0
  1.0000 gamma gamma
  Enddecay
  End
DECAYCARD

decay_card_kl_k = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-
  Enddecay
  Decay D+
  1.0000 K_L0 K+
  Enddecay
  End
DECAYCARD

decay_card_kl_k_pi0 = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-
  Enddecay
  Decay D+
  1.0000 K_L0 K+ pi0
  Enddecay
  Decay pi0
  1.0000 gamma gamma
  Enddecay
  End
DECAYCARD

# --- Exclusive MC: 100k events for each of the four signal modes ---
exMC_ks_k = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToKsK"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_ks_k
  config.cross_section   = :default
end

exMC_ks_k_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToKsKPi0"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_ks_k_pi0
  config.cross_section   = :default
end

exMC_kl_k = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToKlK"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_kl_k
  config.cross_section   = :default
end

exMC_kl_k_pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToKlKPi0"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_kl_k_pi0
  config.cross_section   = :default
end

### ------------------------------ Event selection (tag-based) ------------------------------ ###
# Single-tag pattern (one tag_side call) is inferred from the declarations.
# Store-not-cut: tag mBC / ΔE / missing mass are stored, windowed later in ROOT.
# Charged tracks are constrained only by species counts + net charge (no |cosθ|/Vz/Vr/PID block).

# ============================ Mode 1: D⁺ → K_S0 K⁺ ============================
alg_ks_k = TagAnalysis.new("DpToKsK_ST")
alg_ks_k.set_header(["DpToKsK_STAlg/DpToKsK_ST.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .with_decay_card(decay_card_ks_k)

alg_ks_k.tag_side(:Dplus) do |t|
  # six D⁻ hadronic tag modes (DTagAlg channel names)
  t.modes :DplusToKPiPi,       # K+ π- π-
          :DplusToKPiPiPi0,    # K+ π- π- π0
          :DplusToKsPi,        # K_S0 π-
          :DplusToKsPiPi0,     # K_S0 π- π0
          :DplusToKsPiPiPi,    # K_S0 π+ π- π-
          :DplusToKKPi         # K+ K- π-
  t.charm -1                   # pin the tagged side to D⁻
end

alg_ks_k.signal_side do |s|
  s.charged(kp: 1, pip: 1, pim: 1)  # K⁺ and the K_S0 → π⁺π⁻ daughters
  s.require_charge(1)               # signal net charge = +1
end

alg_ks_k.fit do |f|
  f.constrain_four_momentum                                              # 4-momentum conservation
  f.invariant_mass_of(:pip, :pim).between(0.4856, 0.5096)                # M(π⁺π⁻) within ±12 MeV of m(K_S0)
  f.chi2_cut 200
end

alg_ks_k
  .note(:secondary_vertex, "K_S0 is built from a π+π- secondary vertex; decay-length significance > 2σ required. The tag-based signal side consumes the tag's unused tracks and has no secondary-vertex primitive, so the vertex fit and L/σ_L cut are emitted inside the generated BOSS selection rather than declared in the DSL.")
  .note(:tag_deltae_window, "Tag-mode-dependent ΔE requirements (Table I) applied per tag mode. Mode-dependent windows are not expressible as a single tag_side window, so ΔE_tag is stored unconditionally and cut in ROOT.")
  .apply

alg_ks_k.execute_on([psip_data, psip_incMC, exMC_ks_k])

# ========================= Mode 2: D⁺ → K_S0 K⁺ π⁰ ==========================
alg_ks_k_pi0 = TagAnalysis.new("DpToKsKPi0_ST")
alg_ks_k_pi0.set_header(["DpToKsKPi0_STAlg/DpToKsKPi0_ST.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            .with_decay_card(decay_card_ks_k_pi0)

alg_ks_k_pi0.tag_side(:Dplus) do |t|
  t.modes :DplusToKPiPi, :DplusToKPiPiPi0, :DplusToKsPi,
          :DplusToKsPiPi0, :DplusToKsPiPiPi, :DplusToKKPi
  t.charm -1
end

alg_ks_k_pi0.signal_side do |s|
  s.photons 2                       # π⁰ → γγ
  s.min_photon_energy 0.025         # E(γ) > 25 MeV
  s.min_photon_angle 10.0           # photon opening/isolating angle > 10°
  s.charged(kp: 1, pip: 1, pim: 1)  # K⁺ and the K_S0 → π⁺π⁻ daughters
  s.require_charge(1)               # signal net charge = +1
end

alg_ks_k_pi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).between(0.4856, 0.5096)                # M(π⁺π⁻) within ±12 MeV of m(K_S0)
  f.invariant_mass_of(:gamma, :gamma).between(0.110, 0.155)             # M(γγ) π⁰ window
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0) # γγ constrained to m(π⁰) in the fit
  f.chi2_cut 200
end

alg_ks_k_pi0
  .note(:tag_deltae_window, "Tag-mode-dependent ΔE requirements (Table I); applied per tag mode in ROOT (ΔE_tag stored unconditionally).")
  .note(:secondary_vertex, "K_S0 built from a π+π- secondary vertex with decay-length significance > 2σ; emitted inside the generated BOSS selection.")
  .note(:pi0_mass_chi2, "π⁰ candidates additionally required to have χ² to the nominal π⁰ mass < 20; applied on the γγ Kalman mass-constraint output, which has no dedicated DSL knob on the tag signal side.")
  .note(:signal_deltae_window, "Signal ΔE required in [-0.057, 0.040] GeV; applied in ROOT on the stored signal ΔE.")
  .apply

alg_ks_k_pi0.execute_on([psip_data, psip_incMC, exMC_ks_k_pi0])

# ============================ Mode 3: D⁺ → K_L0 K⁺ ==========================
alg_kl_k = TagAnalysis.new("DpToKlK_ST")
alg_kl_k.set_header(["DpToKlK_STAlg/DpToKlK_ST.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .with_decay_card(decay_card_kl_k)

alg_kl_k.tag_side(:Dplus) do |t|
  t.modes :DplusToKPiPi, :DplusToKPiPiPi0, :DplusToKsPi,
          :DplusToKsPiPi0, :DplusToKsPiPiPi, :DplusToKKPi
  t.charm -1
end

alg_kl_k.signal_side do |s|
  s.charged(kp: 1)                  # K⁺
  s.require_charge(1)               # signal net charge = +1
  s.missing :K_L0                   # massive missing K_L0 of known mass
  s.min_photon_energy 0.1           # EMC shower (>0.1 GeV) fixes the K_L0 direction
end

alg_kl_k.fit do |f|
  f.constrain_four_momentum         # K_L0 momentum inferred from ΔE_sig = 0
  f.chi2_cut 200
end

alg_kl_k
  .note(:tag_deltae_window, "Tag-mode-dependent ΔE requirements (Table I); applied per tag mode in ROOT (ΔE_tag stored unconditionally).")
  .note(:kl_direction, "K_L0 direction taken from the EMC shower with E > 0.1 GeV; momentum inferred from ΔE_sig = 0 through the 4-momentum constraint. Best candidate chosen by minimum fit χ² (DSL default).")
  .apply

alg_kl_k.execute_on([psip_data, psip_incMC, exMC_kl_k])

# ========================= Mode 4: D⁺ → K_L0 K⁺ π⁰ ==========================
alg_kl_k_pi0 = TagAnalysis.new("DpToKlKPi0_ST")
alg_kl_k_pi0.set_header(["DpToKlKPi0_STAlg/DpToKlKPi0_ST.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            .with_decay_card(decay_card_kl_k_pi0)

alg_kl_k_pi0.tag_side(:Dplus) do |t|
  t.modes :DplusToKPiPi, :DplusToKPiPiPi0, :DplusToKsPi,
          :DplusToKsPiPi0, :DplusToKsPiPiPi, :DplusToKKPi
  t.charm -1
end

alg_kl_k_pi0.signal_side do |s|
  s.photons 2                       # π⁰ → γγ
  s.min_photon_energy 0.025         # E(γ) > 25 MeV
  s.min_photon_angle 10.0           # photon opening/isolating angle > 10°
  s.charged(kp: 1)                  # K⁺
  s.require_charge(1)               # signal net charge = +1
  s.missing :K_L0                   # massive missing K_L0 of known mass
end

alg_kl_k_pi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).between(0.110, 0.155)              # M(γγ) π⁰ window
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0) # γγ constrained to m(π⁰) in the fit
  f.chi2_cut 200
end

alg_kl_k_pi0
  .note(:tag_deltae_window, "Tag-mode-dependent ΔE requirements (Table I); applied per tag mode in ROOT (ΔE_tag stored unconditionally).")
  .note(:pi0_mass_chi2, "π⁰ candidates additionally required to have χ² to the nominal π⁰ mass < 20; applied on the γγ Kalman mass-constraint output.")
  .note(:kl_direction, "K_L0 direction from the EMC shower with E > 0.1 GeV, momentum inferred from ΔE_sig = 0 via the 4-momentum constraint; best candidate chosen by minimum fit χ² (DSL default). The 0.1 GeV K_L0 shower floor is applied on top of the 25 MeV π⁰ photon floor inside the generated selection.")
  .apply

alg_kl_k_pi0.execute_on([psip_data, psip_incMC, exMC_kl_k_pi0])

# -----------------------------------------------------------------------------
# Out of scope (ROOT analysis part): unbinned 2D fit to M_BC(tag) vs M_BC(signal),
# K_S0 (~2%) / K_L0 (~10%) efficiency corrections, and the fixed peaking
# D⁺ → K_S0 K⁺(π⁰), K_S0 → π⁰π⁰ backgrounds (~3% / ~5%).
# -----------------------------------------------------------------------------