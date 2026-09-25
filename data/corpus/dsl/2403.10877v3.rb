# BESIII DSL for arXiv:2403.10877v3
# Test of Lepton Universality: D0 → K-π0 μ+ ν_μ
# Tag-based (ST/DT): ST anti-D0, signal D0 semileptonic with missing nu_mu
# √s = 3.773 GeV, ψ(3770) → D0 D0bar

# ============================================================
# Datasets
# ============================================================
psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# ============================================================
# Decay card for ψ(3770) → D0 D0bar
# ============================================================
decay_card_psip3770 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    End
DECAYCARD

exMC_tag_sl = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "d0tag_semileptonic"
  config.related_dataset = psi3770_data
  config.events = 500_000
  config.decay_card = decay_card_psip3770
  config.cross_section = :default
end

# ============================================================
# TagAnalysis: ST anti-D0 + signal D0 → K-π0 μ+ ν_μ
# Tag modes from Table I:
#   K+π-, K+π-π0, K+π+π-π-, K_S0π+π-, K+π-π-π0, K+π+π-π-π0
# Signal side: K-, π0, μ+, missing nu_mu
# ============================================================
alg_d0_sl = TagAnalysis.new("D0ToKPi0MuNu")
alg_d0_sl.set_header(["D0ToKPi0MuNuAlg/D0ToKPi0MuNu.h"])
         .set_constant({"ECMS" => [:double, 3.773]})

# Tag side: anti-D0 → hadronic modes (charm -1 = anti-D0 tag)
alg_d0_sl.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKsPiPi,
          :D0toKPiPiPi0, :D0toKPiPiPiPi0
  t.charm -1
end

# Signal side: D0 → K- π0 μ+ ν_μ
# K- (km), π0 → γγ (2 photons), μ+ (mup), missing ν_μ
alg_d0_sl.signal_side do |s|
  s.photons 2                    # pi0 → γγ
  s.charged(km: 1, mup: 1)       # K- and μ+
  s.require_charge 0              # -1 + 1 = 0
  s.missing :nu_mu                # massless neutrino
end

# 4C fit: tag + signal (K + pi0 + mu + nu) = ecms_lab
# pi0 mass constraint from the 2 photons
alg_d0_sl.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

# Apply the tag analysis
alg_d0_sl.apply

# ============================================================
# Execute
# ============================================================
alg_d0_sl.execute_on([psi3770_data, psi3770_incMC, exMC_tag_sl])