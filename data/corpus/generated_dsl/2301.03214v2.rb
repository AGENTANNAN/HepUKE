# frozen_string_literal: true

### Dataset description ###
# ψ(3770) real data and matching inclusive MC (2.93 fb⁻¹, sample 712_3773)
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ---------------------------------------------------------------------
# Decay cards (EvtGen format) — top mother is the ψ(3770) resonance
# ---------------------------------------------------------------------
# D0 channel: flavor tag anti-D0 → K+ π-, signal D0 → π+ π+ π- X
decay_card_D0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay D0
    1.0000 pi+ pi+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# D+ channel: hadronic tag D- → K+ π- π-, signal D+ → π+ π+ π- X
decay_card_Dp = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay D+
    1.0000 pi+ pi+ pi- pi0 PHSP;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------
# Exclusive MC samples — 1,000,000 events per channel
# ---------------------------------------------------------------------
exMC_D0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_D0_KPi_3piX"
  config.related_dataset = psi3770_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_D0
  config.cross_section   = :default
end

exMC_Dp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3770_Dp_KPiPi_3piX"
  config.related_dataset = psi3770_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_Dp
  config.cross_section   = :default
end

# ---------------------------------------------------------------------
# Tag-based event selection (TagAnalysis)
# ---------------------------------------------------------------------

# =============================== D0 channel ===============================
alg_D0 = TagAnalysis.new("D0TagPiPiPiX")
alg_D0.set_header(["D0TagPiPiPiXAlg/D0TagPiPiPiX.h"])
      .set_constant({ "ECMS" => [:double, 3.773] })
      .with_decay_card(decay_card_D0)
      # C = -1 ψ(3770) → D0 anti-D0 production is quantum-coherent; the tag/signal
      # efficiency and the M_3π migration matrix must use the coherent treatment.
      .note(:quantum_coherence_correction,
            "quantum-coherence correction applied for the C = -1 psi(3770) tag correlation " \
            "in the D0 channel (anti-D0 -> K+ pi- flavor tag); the coherent D0 anti-D0 " \
            "production modifies the effective tag/signal efficiency and the M_3pi " \
            "migration matrix relative to the incoherent assumption")
      # K_S0 veto — a "reject pairs inside a window" selection has no DSL form
      # (tag fit only supports invariant_mass_of(...).between, i.e. a window requirement).
      .note(:background_veto,
            "K_S0 contribution vetoed: events in which any signal-side pi+ pi- pair has " \
            "invariant mass in [0.485, 0.510] GeV/c2 are rejected")

# Tag side: anti-D0 → K+ pi- (flavor tag). One tag_side call ⇒ single tag (ST).
alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi   # D0 -> K- pi+ mode, charge-conjugated to the anti-D0
  t.charm -1         # pin the tagged side to the anti-D0
end

# Signal side: inclusive D0 → pi+ pi+ pi- X
alg_D0.signal_side do |s|
  s.charged(pip: 2, pim: 1, at_least: true)  # at least two pi+ and one pi-
  s.require_charge 1                         # net charge +1 of the pi+pi+pi- system
end

# Four-momentum-constrained fit, chi2 < 200
alg_D0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_D0.apply
alg_D0.execute_on([psi3770_data, psi3770_incMC, exMC_D0])

# =============================== D+ channel ===============================
alg_Dp = TagAnalysis.new("DpTagPiPiPiX")
alg_Dp.set_header(["DpTagPiPiPiXAlg/DpTagPiPiPiX.h"])
      .set_constant({ "ECMS" => [:double, 3.773] })
      .with_decay_card(decay_card_Dp)
      # No quantum-correlation correction in the D+ channel.
      .note(:background_veto,
            "K_S0 contribution vetoed: events in which any signal-side pi+ pi- pair has " \
            "invariant mass in [0.485, 0.510] GeV/c2 are rejected")

# Tag side: D- → K+ pi- pi- (hadronic tag). Single tag.
alg_Dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi  # D+ -> K- pi+ pi+ mode, charge-conjugated to the D-
  t.charm -1          # pin the tagged side to the D-
end

# Signal side: inclusive D+ → pi+ pi+ pi- X
alg_Dp.signal_side do |s|
  s.charged(pip: 2, pim: 1, at_least: true)  # at least two pi+ and one pi-
  s.require_charge 1                         # net charge +1 of the pi+pi+pi- system
end

# Four-momentum-constrained fit, chi2 < 200
alg_Dp.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_Dp.apply
alg_Dp.execute_on([psi3770_data, psi3770_incMC, exMC_Dp])