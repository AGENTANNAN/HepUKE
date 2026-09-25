### ==================== Dataset preparation ====================
# ψ(3770) data and inclusive MC (BOSS convention: [BOSS version]_[CMS energy in MeV])
data_3773  = DatasetManager.real_data.find("712_3773")      # ψ(3770) real data, 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC sample

### Decay cards (EvtGen format) ###
# Signal channel: D+ -> gamma e+ nu_e (phase space); the tagged D- decays inclusively
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D-    PHSP;
    Enddecay

    Decay D+
    1.000 gamma e+ nu_e   PHSP;
    Enddecay

    End
DECAYCARD

# Validation channel: D+ -> pi0 e+ nu_e, pi0 -> gamma gamma (phase space)
decay_card_validation = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D-    PHSP;
    Enddecay

    Decay D+
    1.000 pi0 e+ nu_e   PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma   PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples (1M events each) ###
exmc_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToGammaEnu"
  config.related_dataset = data_3773          # anchor to the ψ(3770) real data (beam energy, run range)
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

exmc_validation = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToPi0Enu"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_validation
  config.cross_section   = :default
end

# All datasets the two algorithms are run on (real data, inclusive MC, signal MC, validation MC)
datasets = [data_3773, incMC_3773, exmc_signal, exmc_validation]

### ==================== Signal channel: tagged D- vs D+ -> gamma e+ nu_e ====================
alg_sig = TagAnalysis.new("DpToGammaEnu")
alg_sig.set_header(["DpToGammaEnuAlg/DpToGammaEnu.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })   # c.m. energy 3.773 GeV
       .with_decay_card(decay_card_signal)
       .note(:fsr_recovery, "FSR photons within 5 degrees of the positron track are added
         back into the positron four-momentum before the 4C kinematic fit")
       .note(:gamma_selection, "when more than one signal-side photon passes the shower
         selection, the most energetic photon is retained and used in the fit")

# Tag side: D- (charm = -1) in six hadronic modes, with the loose common |DeltaE| < 60 MeV window.
# Precise DeltaE / mBC windows are applied later in the ROOT analysis (store-not-cut default).
alg_sig.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,        # D- -> K+ pi- pi-
          :DptoKPiPiPi0,     # D- -> K+ pi- pi- pi0
          :DptoKsPi,         # D- -> K_S pi-
          :DptoKsPiPi0,      # D- -> K_S pi- pi0
          :DptoKsPiPiPi,     # D- -> K_S pi- pi- pi+
          :DptoKKPi          # D- -> K+ K- pi-
  t.charm -1
  t.window :deltaE, abs: 0.06   # loose common |DeltaE| < 60 MeV (only sanctioned tag-side cut)
end

# Signal side: exactly one photon, one positron (net charge +1), one massless missing nu_e
alg_sig.signal_side do |s|
  s.photons 1
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

# 4C kinematic fit: tag + gamma + e+ + nu_e constrained to the c.m. energy, chi2 < 200
alg_sig.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_sig.apply                                   # takes no Selection argument
alg_sig.execute_on(datasets)

### ==================== Validation channel: tagged D- vs D+ -> pi0 e+ nu_e ====================
alg_val = TagAnalysis.new("DpToPi0Enu")
alg_val.set_header(["DpToPi0EnuAlg/DpToPi0Enu.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .with_decay_card(decay_card_validation)

# Same D- tag: six hadronic modes, charm = -1, loose common |DeltaE| < 60 MeV window
alg_val.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,
          :DptoKPiPiPi0,
          :DptoKsPi,
          :DptoKsPiPi0,
          :DptoKsPiPiPi,
          :DptoKKPi
  t.charm -1
  t.window :deltaE, abs: 0.06
end

# Signal side: two photons and one positron (net charge +1), one massless missing nu_e
alg_val.signal_side do |s|
  s.photons 2
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

# 4C kinematic fit: M(gamma gamma) window before the fit, pi0 nominal-mass constraint in the fit,
# tag + pi0 + e+ + nu_e constrained to the c.m. energy, chi2 < 200
alg_val.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)                  # M(gamma gamma) in [115, 150] MeV/c^2
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)     # pi0 mass constraint
  f.chi2_cut 200
end

alg_val.apply
alg_val.execute_on(datasets)