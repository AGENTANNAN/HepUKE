### Dataset preparation ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data at 3.773 GeV
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # Corresponding ψ(3770) inclusive MC (used at 10× data)

# ------------------------------------------------------------------
# Decay card attached to the algorithm: the signal side only.
# ψ(3770) -> D⁻ (the tag, reconstructed from the pre-stored DTag) + D⁺ -> μ⁺ν_μ
# ------------------------------------------------------------------
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 D- D+ PHSP;
    Enddecay

    Decay D+
    1.000 mu+ nu_mu PHSP;
    Enddecay

    End
DECAYCARD

# ------------------------------------------------------------------
# Nine tag-mode decay cards, one per hadronic D⁻ tag mode,
# each with the recoiling D⁺ -> μ⁺ν_μ (K_S0 -> π⁺π⁻, π⁰ -> γγ)
# ------------------------------------------------------------------
decay_card_Kpipi = <<~DECAYCARD            # D⁻ -> K⁺π⁻π⁻
    Decay psi(3770)
    1.000 D- D+ PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    Decay D+
    1.000 mu+ nu_mu PHSP;
    Enddecay

    End
DECAYCARD

decay_card_KsPi = <<~DECAYCARD             # D⁻ -> K_S⁰π⁻
    Decay psi(3770)
    1.000 D- D+ PHSP;
    Enddecay

    Decay D-
    1.000 K_S0 pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay D+
    1.000 mu+ nu_mu PHSP;
    Enddecay

    End
DECAYCARD

decay_card_KsK = <<~DECAYCARD              # D⁻ -> K_S⁰K⁻
    Decay psi(3770)
    1.000 D- D+ PHSP;
    Enddecay

    Decay D-
    1.000 K_S0 K- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay D+
    1.000 mu+ nu_mu PHSP;
    Enddecay

    End
DECAYCARD

decay_card_KKpi = <<~DECAYCARD             # D⁻ -> K⁺K⁻π⁻
    Decay psi(3770)
    1.000 D- D+ PHSP;
    Enddecay

    Decay D-
    1.000 K+ K- pi- PHSP;
    Enddecay

    Decay D+
    1.000 mu+ nu_mu PHSP;
    Enddecay

    End
DECAYCARD

decay_card_Kpipipi0 = <<~DECAYCARD         # D⁻ -> K⁺π⁻π⁻π⁰
    Decay psi(3770)
    1.000 D- D+ PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay D+
    1.000 mu+ nu_mu PHSP;
    Enddecay

    End
DECAYCARD

decay_card_pipipim = <<~DECAYCARD          # D⁻ -> π⁺π⁻π⁻
    Decay psi(3770)
    1.000 D- D+ PHSP;
    Enddecay

    Decay D-
    1.000 pi+ pi- pi- PHSP;
    Enddecay

    Decay D+
    1.000 mu+ nu_mu PHSP;
    Enddecay

    End
DECAYCARD

decay_card_KsPipi0 = <<~DECAYCARD          # D⁻ -> K_S⁰π⁻π⁰
    Decay psi(3770)
    1.000 D- D+ PHSP;
    Enddecay

    Decay D-
    1.000 K_S0 pi- pi0 PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay D+
    1.000 mu+ nu_mu PHSP;
    Enddecay

    End
DECAYCARD

decay_card_Kpipipipip = <<~DECAYCARD       # D⁻ -> K⁺π⁻π⁻π⁻π⁺
    Decay psi(3770)
    1.000 D- D+ PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- pi- pi+ PHSP;
    Enddecay

    Decay D+
    1.000 mu+ nu_mu PHSP;
    Enddecay

    End
DECAYCARD

decay_card_KsPipipim = <<~DECAYCARD        # D⁻ -> K_S⁰π⁻π⁻π⁺
    Decay psi(3770)
    1.000 D- D+ PHSP;
    Enddecay

    Decay D-
    1.000 K_S0 pi- pi- pi+ PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay D+
    1.000 mu+ nu_mu PHSP;
    Enddecay

    End
DECAYCARD

# ------------------------------------------------------------------
# Create a 100k-event exclusive MC sample for each of the nine tag modes
# ------------------------------------------------------------------
tag_mode_cards = {
  "Kpipi"        => decay_card_Kpipi,
  "KsPi"         => decay_card_KsPi,
  "KsK"          => decay_card_KsK,
  "KKpi"         => decay_card_KKpi,
  "Kpipipi0"     => decay_card_Kpipipi0,
  "pipipim"      => decay_card_pipipim,
  "KsPipi0"      => decay_card_KsPipi0,
  "Kpipipipip"   => decay_card_Kpipipipip,
  "KsPipipim"    => decay_card_KsPipipim
}

exMC_tags = tag_mode_cards.map do |tag_name, card|
  DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "dtag_munu_#{tag_name}"
    config.related_dataset = psi3770_data
    config.events          = 100000
    config.decay_card      = card
    config.cross_section   = :default
  end
end

### Event selection (BOSS) — tag analysis ###
alg = TagAnalysis.new("DpToMuNuTag")
alg.set_header(["DpToMuNuTagAlg/DpToMuNuTag.h"])
   .set_constant({"ECMS" => [:double, 3.773]})

# Single tag: one D⁻ reconstructed in nine hadronic modes (m_BC and ΔE stored, windowed later in ROOT)
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKsK, :DptoKKPi,
          :DptoKPiPiPi0, :DptoPiPiPi, :DptoKsPiPi0,
          :DptoKPiPiPiPi, :DptoKsPiPiPi
  t.charm -1   # pin the tagged side to D⁻
end

# Signal side: exactly one leftover track (μ⁺) and one undetected massless ν_μ
alg.signal_side do |s|
  s.charged(mup: 1)     # exactly one signal-side track, identified as μ⁺
  s.missing :nu_mu      # massless missing ν_μ
end

# 4C kinematic fit: tag + μ⁺ + ν_μ to the measured CMS four-momentum
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# BOSS-side procedure that has no dedicated DSL construct
alg.note(:pid_correction_method,
         "π/K separation: CL_π > CL_K (CL_K > CL_π) below 0.75 GeV/c and CL > 0.1% above; "
         "handled by DTagAlg's internal PID, not tunable from the DSL")
   .note(:track_vertex_constraint,
         "all charged tracks except K_S0 daughters required DCA < 1.0 cm transverse and "
         "< 15.0 cm along the beam, and constrained to a common vertex")
   .note(:background_veto,
         "any additional good EMC photon required E_gamma_max < 300 MeV")
   .note(:pi0_reconstruction,
         "π⁰ mass window from a 1C (Kalman) fit with χ² < 100")
   .note(:muon_id,
         "signal-side μ⁺ identified via MUC penetration depth; DTagTool signal-side lepton "
         "PID uses fixed v1 thresholds")

alg.with_decay_card(decay_card_signal).apply

# Execute on ψ(3770) real data, inclusive MC, and the nine tag-mode exclusive MC samples
root_files = alg.execute_on([psi3770_data, psi3770_incMC] + exMC_tags)