# ============================================================================
# Semileptonic D decays at sqrt(s) = 3.773 GeV (psi(3770))
# Double-tag analysis of the three signal modes
#   I  : D0 -> K- eta e+ nu_e          (eta -> gamma gamma)
#   II : D+ -> K_S0 eta e+ nu_e        (K_S0 -> pi+ pi-, eta -> gamma gamma)
#   III: D+ -> eta eta e+ nu_e         (eta -> gamma gamma)
# The tag D is reconstructed from DTagTool (hadronic modes) -> TagAnalysis.
# ============================================================================

### Datasets ###
data_3773  = DatasetManager.real_data.find("712_3773")       # 3.773 GeV real data (psi(3770))
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")    # matching inclusive MC

### Decay cards ###
# Mode I: D0 -> K- eta e+ nu_e (tag side shown as a representative hadronic mode)
decay_card_modeI = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0 K- eta e+ nu_e PHSP;
  Enddecay

  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay

  Decay anti-D0
  1.0 K+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: D+ -> K_S0 eta e+ nu_e
decay_card_modeII = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D- PHSP;
  Enddecay

  Decay D+
  1.0 K_S0 eta e+ nu_e PHSP;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay

  Decay D-
  1.0 K+ pi- pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode III: D+ -> eta eta e+ nu_e
decay_card_modeIII = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D- PHSP;
  Enddecay

  Decay D+
  1.0 eta eta e+ nu_e PHSP;
  Enddecay

  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay

  Decay D-
  1.0 K+ pi- pi- PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC (100k events for each signal mode) ###
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_D0_Km_eta_e_nu"
  config.related_dataset = data_3773
  config.events         = 100_000
  config.decay_card     = decay_card_modeI
  config.cross_section  = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_Dp_Ks_eta_e_nu"
  config.related_dataset = data_3773
  config.events         = 100_000
  config.decay_card     = decay_card_modeII
  config.cross_section  = :default
end

exMC_modeIII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_Dp_eta_eta_e_nu"
  config.related_dataset = data_3773
  config.events         = 100_000
  config.decay_card     = decay_card_modeIII
  config.cross_section  = :default
end

# ============================================================================
# Mode I : D0 -> K- eta e+ nu_e   (tag: anti-D0 from D0 hadronic modes)
# ============================================================================
alg_modeI = TagAnalysis.new("D0KmEtaENu")
alg_modeI.set_header(["D0KmEtaENuAlg/D0KmEtaENu.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_modeI)

# Tag side: opposite-charm D0 (i.e. anti-D0) reconstructed from hadronic modes
alg_modeI.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1                       # tag D carries opposite charm to the signal D0
end

# Signal side: two photons (eta -> gamma gamma) + K- + e+ , missing nu_e
alg_modeI.signal_side do |s|
  s.photons 2
  s.charged(km: 1, ep: 1)
  s.missing :nu_e                  # massless neutrino
  s.require_charge 0               # zero net charge (K- e+)
  s.min_photon_angle 10.0          # photon-track angle > 10 degrees
  s.min_photon_energy 0.025        # photon energy > 25 MeV
end

# Kinematic fit: 4-momentum conservation + gamma-gamma mass constrained to eta
alg_modeI.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma).between(0.505, 0.575)   # gamma-gamma mass window (modes I-II)
  f.chi2_cut 200
end

alg_modeI
  .note(:pid_correction_method,
        "signal-side kaon/pion separated by likelihood ratios; electron selected with " \
        "L_e>0.001, L_e/(L_e+L_pi+L_K)>0.8 and E/p in [0.8,1.1] (fixed v1 PID thresholds, not DSL-tunable)")
  .note(:photon_selection,
        "photon candidates: barrel E>25 MeV (|cos theta|<0.80) or endcap E>50 MeV (0.86<|cos theta|<0.92), " \
        "EMC TDC in [0,700] ns")
  .note(:background_veto,
        "ROOT-level vetoes: M(P,eta,e)<1.80 GeV, extra-photon energy <0.5 GeV, extra pi0 veto " \
        "(extra pi0 built from gamma-gamma in [0.115,0.150] GeV with 1C-fit chi2<50); " \
        "mode I additionally requires M(K,eta)>1.30 GeV")
  .apply

alg_modeI.execute_on([data_3773, incMC_3773, exMC_modeI])

# ============================================================================
# Mode II : D+ -> K_S0 eta e+ nu_e   (tag: D- from D+ hadronic modes)
# ============================================================================
alg_modeII = TagAnalysis.new("DpKsEtaENu")
alg_modeII.set_header(["DpKsEtaENuAlg/DpKsEtaENu.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .with_decay_card(decay_card_modeII)

# Tag side: opposite-charm D+ (i.e. D-) reconstructed from hadronic modes
alg_modeII.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1                       # tag D carries opposite charm to the signal D+
end

# Signal side: two photons (eta) + pi+ pi- (K_S0) + e+ , missing nu_e
alg_modeII.signal_side do |s|
  s.photons 2
  s.charged(pip: 1, pim: 1, ep: 1)
  s.missing :nu_e
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# Kinematic fit: 4-momentum + gamma-gamma mass constrained to eta
alg_modeII.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma).between(0.505, 0.575)   # gamma-gamma mass window (modes I-II)
  f.invariant_mass_of(:pip, :pim).between(0.487, 0.511)       # K_S0 mass window
  f.chi2_cut 200
end

alg_modeII
  .note(:ks0_reconstruction,
        "K_S0 built from two oppositely charged tracks with |Vz|<20, secondary vertex-fit chi2<100, " \
        "invariant mass in [0.487,0.511] GeV and decay length > 2 sigma (mass window enforced in the fit)")
  .note(:pid_correction_method,
        "signal-side kaon/pion separated by likelihood ratios; electron selected with " \
        "L_e>0.001, L_e/(L_e+L_pi+L_K)>0.8 and E/p in [0.8,1.1] (fixed v1 PID thresholds, not DSL-tunable)")
  .note(:photon_selection,
        "photon candidates: barrel E>25 MeV (|cos theta|<0.80) or endcap E>50 MeV (0.86<|cos theta|<0.92), " \
        "EMC TDC in [0,700] ns")
  .note(:signal_charge,
        "description asks for zero net charge on modes I-II; the mode-II signal side (pi+ pi- e+) carries " \
        "net charge +1 balanced by the D- tag, so no signal-side charge requirement is imposed")
  .note(:background_veto,
        "ROOT-level vetoes: M(P,eta,e)<1.80 GeV, extra-photon energy <0.5 GeV, extra pi0 veto " \
        "(extra pi0 built from gamma-gamma in [0.115,0.150] GeV with 1C-fit chi2<50)")
  .apply

alg_modeII.execute_on([data_3773, incMC_3773, exMC_modeII])

# ============================================================================
# Mode III : D+ -> eta eta e+ nu_e   (tag: D- from D+ hadronic modes)
# ============================================================================
alg_modeIII = TagAnalysis.new("DpEtaEtaENu")
alg_modeIII.set_header(["DpEtaEtaENuAlg/DpEtaEtaENu.h"])
           .set_constant({"ECMS" => [:double, 3.773]})
           .with_decay_card(decay_card_modeIII)

# Tag side: opposite-charm D+ (i.e. D-) reconstructed from hadronic modes
alg_modeIII.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1                       # tag D carries opposite charm to the signal D+
end

# Signal side: four photons (two eta) + e+ , missing nu_e
alg_modeIII.signal_side do |s|
  s.photons 4
  s.charged(ep: 1)
  s.missing :nu_e
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# Kinematic fit: 4-momentum + gamma-gamma mass constrained to eta
alg_modeIII.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_modeIII
  .note(:pid_correction_method,
        "signal-side kaon/pion separated by likelihood ratios; electron selected with " \
        "L_e>0.001, L_e/(L_e+L_pi+L_K)>0.8 and E/p in [0.8,1.1] (fixed v1 PID thresholds, not DSL-tunable)")
  .note(:photon_selection,
        "photon candidates: barrel E>25 MeV (|cos theta|<0.80) or endcap E>50 MeV (0.86<|cos theta|<0.92), " \
        "EMC TDC in [0,700] ns")
  .note(:background_veto,
        "ROOT-level vetoes: M(P,eta,e)<1.80 GeV, extra-photon energy <0.5 GeV, extra pi0 veto " \
        "(extra pi0 built from gamma-gamma in [0.115,0.150] GeV with 1C-fit chi2<50)")
  .apply

alg_modeIII.execute_on([data_3773, incMC_3773, exMC_modeIII])