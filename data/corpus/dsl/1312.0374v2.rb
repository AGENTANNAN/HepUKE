# =====================================================================
# BESIII: Precision measurement of B(D+ -> mu+ nu_mu), f_D+, |V_cd|
# Data: 2.92 fb^-1 at sqrt(s) = 3.773 GeV (psi(3770) peak)
#
# Tag-based analysis: a D- meson is reconstructed from the pre-stored
# EvtRecDTag collection (nine hadronic tag modes); the recoiling D+ is
# then searched for in the purely leptonic decay D+ -> mu+ nu_mu.
# This is a TagAnalysis (single tag + missing neutrino), NOT an
# Algorithm + Selection.
# =====================================================================

### Dataset description ###
psipp_data  = DatasetManager.real_data.find("712_3773")    # psi(3770) data, 2.92 fb^-1
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

# ---------------------------------------------------------------------
# Decay cards (EvtGen)
#
# Combined card: psi(3770) -> D+ D-, the tagged D- proceeds through the
# nine hadronic tag modes (equal fractions used as a generation
# placeholder; the per-mode efficiency weighting is applied in ROOT),
# while the recoiling D+ decays to mu+ nu_mu.
# ---------------------------------------------------------------------
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D-
    0.11111 K+ pi- pi- PHSP;
    0.11111 K_S0 pi- PHSP;
    0.11111 K_S0 K- PHSP;
    0.11111 K+ K- pi- PHSP;
    0.11111 K+ pi- pi- pi0 PHSP;
    0.11111 pi+ pi- pi- PHSP;
    0.11111 K_S0 pi- pi0 PHSP;
    0.11111 K+ pi- pi- pi- pi+ PHSP;
    0.11111 K_S0 pi- pi- pi+ PHSP;
    Enddecay

    Decay D+
    1.0000 mu+ nu_mu PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Per-tag-mode signal MC cards: psi(3770) -> D+ D- with a single tag mode,
# used to determine the detection efficiency of D+ -> mu+ nu_mu in each
# tagged D- mode (the paper's epsilon = 0.6403 +- 0.0012 is the weighted
# average over the nine modes).
tag_mode_cards = {
  "KPiPi"        => ["D-", "K+ pi- pi- PHSP"],
  "KsPi"         => ["D-", "K_S0 pi- PHSP"],
  "KsK"          => ["D-", "K_S0 K- PHSP"],
  "KKPi"         => ["D-", "K+ K- pi- PHSP"],
  "KPiPiPi0"     => ["D-", "K+ pi- pi- pi0 PHSP"],
  "PiPiPi"       => ["D-", "pi+ pi- pi- PHSP"],
  "KsPiPi0"      => ["D-", "K_S0 pi- pi0 PHSP"],
  "KPiPiPiPi"    => ["D-", "K+ pi- pi- pi- pi+ PHSP"],
  "KsPiPiPi"     => ["D-", "K_S0 pi- pi- pi+ PHSP"]
}

exMC_signal = tag_mode_cards.map do |label, (mother, channel)|
  card = <<~DECAYCARD
      Decay psi(3770)
      1.0000 D+ D- PHSP;
      Enddecay

      Decay #{mother}
      1.0000 #{channel};
      Enddecay

      Decay D+
      1.0000 mu+ nu_mu PHSP;
      Enddecay

      Decay K_S0
      1.0000 pi+ pi- PHSP;
      Enddecay

      Decay pi0
      1.0000 gamma gamma PHSP;
      Enddecay

      End
  DECAYCARD

  DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_Dp_to_munu_tag_#{label}"
    config.related_dataset = psipp_data
    config.events          = 100_000
    config.decay_card      = card
    config.cross_section   = :default
  end
end

# ---------------------------------------------------------------------
### Tag analysis (BOSS) ###
# ---------------------------------------------------------------------
alg_name = "DmTagMuNu"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .set_alias({"std::vector<double>" => "Vdouble"})

# Tag side: a single tagged D- (anti-charm, charm = -1) reconstructed from
# the nine hadronic D- decay modes. The tag mBC and deltaE are stored
# unconditionally; the D- mass window (the paper's red dashed lines) and
# the |E_beam - E_mKnpi| < 2.5 sigma_E requirement are applied in ROOT.
alg.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,
          :DptoKsPi,
          :DptoKsK,
          :DptoKKPi,
          :DptoKPiPiPi0,
          :DptoPiPiPi,
          :DptoKsPiPi0,
          :DptoKPiPiPiPi,
          :DptoKsPiPiPi
  t.charm -1
end

# Signal side: exactly one leftover positively charged track, classified as
# a muon (the mu+ of D+ -> mu+ nu_mu), plus the undetected neutrino. The
# remaining D+ is fully accounted for by mu+ + nu_mu, so no other particles
# are permitted. The muon identification exploits the MUC penetration
# depth; the missing mass squared M_miss^2 is formed from the fitted
# missing 4-vector and the tagged-D momentum and is windowed in ROOT.
alg.signal_side do |s|
  s.charged(mup: 1)    # exactly one signal-side muon, no other charged tracks
  s.require_charge 1
  s.missing :nu_mu     # MASSLESS (resolved mass 0) -> p4 AddMissTrack overload
end

# 4C kinematic fit: tag + mu+ + nu_mu = the measured CMS 4-vector (the beam
# energy and boost are read per run from MeasuredEcmsSvc).
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# Procedures that have no DSL construct are recorded for the downstream
# (systematic / ROOT) stage.
alg.note(:background_veto,
         "the maximum energy E_gamma_max of any extra good photon in the EMC is required to be < 300 MeV; no DSL construct for an extra-shower energy ceiling")
alg.note(:background_veto,
         "the D+ -> K_L0 pi+ and D+ -> pi+ pi0 peaking backgrounds, the D+ -> tau+ nu_tau feed-down, and the e+e- -> gamma_ISR psi(3686) / gamma_ISR J/psi / q qbar / tau+ tau- / non-DDbar backgrounds are estimated from MC samples 10x the data size and subtracted; their peaking shapes in M_miss^2 are modelled after scaling the simulated M_miss^2 resolution by a factor 1.194 measured with D+ -> K_S0 pi+ control events")
alg.note(:efficiency_curve,
         "the M_miss^2 resolution difference between data and MC is corrected by the factor 1.194; the |M_miss^2| < 0.12 GeV^2/c^4 signal window (0.5% systematic) and the radiative correction (1.0%) are handled outside the BOSS selection")
alg.note(:track_selection,
         "all charged tracks except those from K_S0 decays are required to have a distance of closest approach to the average interaction point of < 1.0 cm in the plane perpendicular to the beam and < 15.0 cm along the beam; charged tracks are constrained to a common vertex")
alg.note(:pid_correction_method,
         "pion (kaon) identification requires CL_pi > CL_K (CL_K > CL_pi) for p < 0.75 GeV/c and CL_pi > 0.1% (CL_K > 0.1%) for p > 0.75 GeV/c; the pi0 -> gamma gamma mass window uses a 1C kinematic fit with chi2 < 100")

alg.with_decay_card(decay_card_signal).apply

# Execute on data, inclusive MC and the per-tag-mode signal MC samples
alg.execute_on([psipp_data, psipp_incMC] + exMC_signal)
