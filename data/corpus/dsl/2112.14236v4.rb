# 2112.14236v4: Search for D0 → π0 ν νbar at √s = 3.773 GeV.
# ψ(3770) → D0 D0bar, double-tag method.
# Tag: anti-D0 via 3 hadronic modes.
# Signal: D0 → π0 ν νbar (all neutral, no charged tracks on signal side).
# Data-driven K_L background modeling.
# First experimental constraint on charmed-hadron dineutrino decays.

psip3770_data  = DatasetManager.real_data.find("712_3773")
psip3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card: ψ(3770) → D0 anti-D0
decay_card = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay
    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_d0_pi0_nunu"
  config.related_dataset = psip3770_data
  config.events          = 300_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# ── TagAnalysis ───────────────────────────────────────────────────
alg = TagAnalysis.new("D0toPi0NuNubar")

alg.set_header(["D0toPi0NuNubarAlg/D0toPi0NuNubar.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })
   .with_decay_card(decay_card)

# Tag side: anti-D0 via three hadronic decay modes (Ref. [19])
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# Signal side: D0 → π0 ν νbar
# All neutral final state — no charged tracks; exactly 1 π0 from photons
# not used by the tag; missing: two neutrinos (treated as massless missing
# particle with 4-momentum inferred from recoil).
alg.signal_side do |s|
  s.photons 2                   # π0 → γγ
  s.missing :nu                 # massless missing (recoil mass method)
end

# Kinematic fit: 4C + π0 mass constraint
# The fit constrains the two photons to the π0 mass and enforces 4-momentum
# conservation with the tag side + signal side + missing neutrino.
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg
  .note(:signal_definition, "D0 → π0 ν νbar. Signal final state is all neutral. After single-tag D0_bar, require zero additional charged tracks and exactly one π0 from photon pairs not used by the tag. The π0 invariant mass must satisfy M_γγ ∈ (0.115, 0.150) GeV/c² (ROOT-level cut).")
  .note(:pi0_kinematic_fit, "π0 kinematic fit constraining M_γγ to nominal π0 mass; χ² < 20 required for background suppression (ROOT-level cut).")
  .note(:mmiss_cut, "Recoil mass squared cut: M_miss² ∈ (1.1, 1.9) GeV²/c⁴ to remove D0→K_L⁰π0 and D0→K̄*(892)⁰π0 backgrounds (ROOT-level cut). M_miss² = (p_e+e- − p_tag − p_π0)².")
  .note(:data_driven_kl_bkg, "Background dominated by D0 → π0 K_L⁰ X. K_L⁰ energy deposit modeled via data-driven method: control samples J/ψ→φK±π∓K_L⁰ and J/ψ→K±π∓K_L⁰ (purity ~99%). Showers within 10° cone along K_L⁰ flight direction summed. Unaccounted showers (<50 MeV) corrected via multi-dimensional reweighting. K_L⁰ momentum distribution reweighted to match signal simulation.")
  .note(:data_driven_x_bkg, "X component (soft π/K) modeled with D0 → π0 K_S⁰ X, K_S⁰ → π+π- control sample. Showers not associated with tag tracks, π0, or K_S⁰ daughters summed as E_EMC^X.")
  .note(:eemc_fit, "Summed EMC energy E_EMC (all showers excluding those used for tag and signal π0) used as discriminating variable. Signal peaks near zero. Extended ML fit with signal shape from MC, K_L⁰X background from data-driven sampling, wrong-tag background from M_BC sideband (factor 0.611±0.001), and other D0 decays from ψ(3770)→D0D0 simulation.")
  .note(:wrong_tag_bkg, "Wrong-tag (non-D0D0bar) events estimated from M_BC sideband region (1.830, 1.855) GeV/c². Signal region: (1.858, 1.874) GeV/c². Scaling factor 0.611±0.001. Fixed to 1919±34 in fit.")
  .note(:upper_limit, "Branching fraction upper limit: 2.1×10⁻⁴ at 90% CL. Ensemble of 100k toy samples incorporating systematic uncertainties. No significant signal observed (N_sig = 14±30).")
  .note(:branching_formula, "B_sig = N_sig / [B(π0→γγ) · Σ_α N_tag^α · (ε_tag,sig^α / ε_tag^α)]. α = {Kπ, Kππ0, K3π} tag modes. Single-tag yields from M_BC fits; tag efficiencies from generic MC.")
  .apply

alg.execute_on([psip3770_data, psip3770_incMC, exMC_signal])