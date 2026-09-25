# Core DSL classes and dependencies will be loaded automatically at execution
# ---------------------------------------------------------------------------
# Tag-based analysis: e+e- -> psi(4260) -> D_s*+ D_s-,  D_s*+ -> D_s+ gamma
#   * tag side   : D_s-  (via K_S K- and K+ K- pi-)
#   * signal side: D_s+ fully reconstructed
#         Mode I : D_s+ -> omega pi+   (omega -> pi+pi-pi0 Dalitz, pi0 -> gamma gamma)
#         Mode II: D_s+ -> omega K+    (omega -> pi+pi-pi0 Dalitz, pi0 -> gamma gamma)
# Data: sqrt(s) = 4.178 GeV, 3.19 fb^-1  ->  sample "703_4180"
# ---------------------------------------------------------------------------

### Dataset description ###
data_4178  = DatasetManager.real_data.find("703_4180")      # 4.178 GeV real data (3.19 fb^-1)
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")   # corresponding inclusive MC

### Decay cards (signal exclusive MC, EvtGen format) ###
# The tag side D_s- -> K_S K- is included so that a complete e+e- -> D_s*+ D_s- event is generated.
decay_card_omega_pi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.0000 D_s+ gamma PHSP;
    Enddecay

    Decay D_s+
    1.0000 omega pi+ PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay D_s-
    1.0000 K_S0 K- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_omega_k = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.0000 D_s+ gamma PHSP;
    Enddecay

    Decay D_s+
    1.0000 omega K+ PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay D_s-
    1.0000 K_S0 K- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples (2,000,000 events each) ###
exMC_omega_pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_dstag_omega_pi"
  config.related_dataset = data_4178
  config.events          = 2_000_000
  config.decay_card      = decay_card_omega_pi
  config.cross_section   = :default
end

exMC_omega_k = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_dstag_omega_k"
  config.related_dataset = data_4178
  config.events          = 2_000_000
  config.decay_card      = decay_card_omega_k
  config.cross_section   = :default
end

### Tag analysis — Mode I: D_s+ -> omega pi+ ###
alg_omega_pi = TagAnalysis.new("DsTagOmegaPi")
alg_omega_pi.set_header(["DsTagOmegaPiAlg/DsTagOmegaPi.h"])
            .set_constant({"ECMS" => [:double, 4.178]})
            .with_decay_card(decay_card_omega_pi)

# Tag side: single D_s- tag through K_S K- and K+ K- pi- modes
alg_omega_pi.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi     # K_S K- and K+ K- pi-
  t.charm -1                      # pin the tagged side to D_s-
end

# Signal side: D_s+ -> omega pi+  (pi0 -> gamma gamma)
alg_omega_pi.signal_side do |s|
  s.photons 2                     # two photons from pi0 -> gamma gamma
  s.min_photon_angle 10.0         # minimum photon angle to charged tracks
  s.charged(pip: 2, pim: 1)       # omega -> pi+ pi- pi0  +  pi+ from D_s+ -> omega pi+
  s.require_charge 1              # net charge +1
end

# 4-momentum-constrained kinematic fit, chi2 < 200
alg_omega_pi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_omega_pi.apply
alg_omega_pi.execute_on([data_4178, incMC_4178, exMC_omega_pi])

### Tag analysis — Mode II: D_s+ -> omega K+ ###
alg_omega_k = TagAnalysis.new("DsTagOmegaK")
alg_omega_k.set_header(["DsTagOmegaKAlg/DsTagOmegaK.h"])
           .set_constant({"ECMS" => [:double, 4.178]})
           .note(:background_veto, "Suppress the D_s+ -> K_S0 K+ pi0 background in the omega K+ mode by vetoing K_S0 candidates: reject events with |M(pi+pi-) - m_K_S0| < 0.03 GeV/c^2 AND K_S0 flight-length significance L/sigma_L > 2. The flight-length-significance requirement has no DSL representation.")
           .with_decay_card(decay_card_omega_k)

# Tag side: single D_s- tag through K_S K- and K+ K- pi- modes
alg_omega_k.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi     # K_S K- and K+ K- pi-
  t.charm -1                      # pin the tagged side to D_s-
end

# Signal side: D_s+ -> omega K+  (pi0 -> gamma gamma)
alg_omega_k.signal_side do |s|
  s.photons 2                     # two photons from pi0 -> gamma gamma
  s.min_photon_angle 10.0         # minimum photon angle to charged tracks
  s.charged(kp: 1, pip: 1, pim: 1)  # omega -> pi+ pi- pi0  +  K+ from D_s+ -> omega K+
  s.require_charge 1              # net charge +1
end

# 4-momentum-constrained kinematic fit, chi2 < 200
alg_omega_k.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_omega_k.apply
alg_omega_k.execute_on([data_4178, incMC_4178, exMC_omega_k])