# =============================================================================
# Dataset preparation + BOSS event selection for the absolute-BF measurement of
# Lambda_c+ (twelve Cabibbo-favoured modes) at sqrt(s) = 4.599 GeV, using the
# double-tag technique on the pre-stored tag collection.
#
# Tag-based analysis -> TagAnalysis (NOT Algorithm + Selection).
# The tag carries its own selected / PID'd tracks and showers; the signal side
# is what the tag did not use.  No select_track / select_photon / pid here.
# =============================================================================

### ---------------------------- Datasets ---------------------------- ###
data_4599  = DatasetManager.real_data.find("703_4600")      # 567 pb^-1 at sqrt(s) = 4.599 GeV
incMC_4599 = DatasetManager.inclusive_mc.find("703_4600")   # matching inclusive MC

### ------- Decay cards (EvtGen) — one per signal mode; the recoil ------ ###
### ------- anti-Lambda_c- decays generically (left unspecified) -------- ###

# Mode  : p K_S0   (Lambda_c+ -> p K_S0)
decay_card_pKs = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 p+ K_S0 PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode : p K- pi+   (Lambda_c+ -> p K- pi+)
decay_card_pKPi = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 p+ K- pi+ PHSP;
  Enddecay

  End
DECAYCARD

# Mode : p K_S0 pi0
decay_card_pKsPi0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 p+ K_S0 pi0 PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode : p K_S0 pi+ pi-
decay_card_pKsPiPi = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 p+ K_S0 pi+ pi- PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode : p K- pi+ pi0
decay_card_pKPiPi0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 p+ K- pi+ pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode : Lambda pi+
decay_card_LambdaPi = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 Lambda0 pi+ PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode : Lambda pi+ pi0
decay_card_LambdaPiPi0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 Lambda0 pi+ pi0 PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode : Lambda pi+ pi- pi+
decay_card_LambdaPiPiPi = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 Lambda0 pi+ pi- pi+ PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode : Sigma0 pi+
decay_card_Sigma0Pi = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 Sigma0 pi+ PHSP;
  Enddecay

  Decay Sigma0
  1.000 Lambda0 gamma PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode : Sigma+ pi0
decay_card_SigmaPi0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 Sigma+ pi0 PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode : Sigma+ pi+ pi-
decay_card_SigmaPiPi = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 Sigma+ pi+ pi- PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode : Sigma+ omega
decay_card_SigmaOmega = <<~DECAYCARD
  Decay psi(4260)
  1.000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.000 Sigma+ omega PHSP;
  Enddecay

  Decay Sigma+
  1.000 p+ pi0 PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

### ---- Exclusive signal MC: 200k events per mode, recoil generic ---- ###
signal_cards = {
  "exmc_4599_LambdacPtoKsP"        => decay_card_pKs,
  "exmc_4599_LambdacPtoKPiP"       => decay_card_pKPi,
  "exmc_4599_LambdacPtoKsPPi0"     => decay_card_pKsPi0,
  "exmc_4599_LambdacPtoKsPPiPi"    => decay_card_pKsPiPi,
  "exmc_4599_LambdacPtoKPiPPi0"    => decay_card_pKPiPi0,
  "exmc_4599_LambdacPtoLambdaPi"   => decay_card_LambdaPi,
  "exmc_4599_LambdacPtoLambdaPiPi0" => decay_card_LambdaPiPi0,
  "exmc_4599_LambdacPtoLambdaPiPiPi" => decay_card_LambdaPiPiPi,
  "exmc_4599_LambdacPtoSigma0Pi"   => decay_card_Sigma0Pi,
  "exmc_4599_LambdacPtoSigmaPi0"   => decay_card_SigmaPi0,
  "exmc_4599_LambdacPtoSigmaPiPi"  => decay_card_SigmaPiPi,
  "exmc_4599_LambdacPtoSigmaOmega" => decay_card_SigmaOmega
}

exMCs = signal_cards.map do |sample_name, card|
  DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = sample_name      # one exclusive MC per signal mode
    config.related_dataset = data_4599        # same 4.599 GeV data point
    config.events          = 200_000          # 200k events per mode
    config.decay_card      = card             # mode-specific decay card
    config.cross_section   = :default
  end
end

### ----------------------- Event selection (BOSS) ----------------------- ###
alg_name = "LambdacDoubleTag"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])          # C++ header include
   .set_constant({"ECMS" => [:double, 4.599]})             # sqrt(s) = 4.599 GeV

# Tag side: reconstruct ONE Lambda_c from the pre-stored tag collection over all
# twelve Cabibbo-favoured modes, with the charm (charge) fixed to the
# anti-Lambda_c- (charm = -1).  A single tag_side call => single-tag pattern.
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP,            # p K_S0
          :LambdacPtoKPiP,           # p K- pi+
          :LambdacPtoKsPPi0,         # p K_S0 pi0
          :LambdacPtoKsPPiPi,        # p K_S0 pi+ pi-
          :LambdacPtoKPiPPi0,        # p K- pi+ pi0
          :LambdacPtoLambdaPi,       # Lambda pi+
          :LambdacPtoLambdaPiPi0,    # Lambda pi+ pi0
          :LambdacPtoLambdaPiPiPi,   # Lambda pi+ pi- pi+
          :LambdacPtoSigma0Pi,       # Sigma0 pi+
          :LambdacPtoSigmaPi0,       # Sigma+ pi0
          :LambdacPtoSigmaPiPi,      # Sigma+ pi+ pi-
          :LambdacPtoSigmaOmega      # Sigma+ omega
  t.charm -1                         # tagged side is the anti-Lambda_c-
end

# Signal side: the recoil Lambda_c+ reconstructed from the tracks / showers the
# tag did not use.  Mode-specific charged multiplicity (relaxed to >=), net
# charge +1, and 0 / 1 / 2 (up to 4 for the two-pi0 modes) photons.
alg.signal_side do |s|
  s.charged(at_least: true)          # mode-specific charged multiplicity
  s.require_charge(1)                # recoil Lambda_c+ has net charge +1
  s.photons 0..4                     # 0/1/2 photons, up to 4 for the two-pi0 modes
  s.min_photon_angle  10.0
  s.min_photon_energy 0.025
end

# Kinematic fit: 4-momentum conservation + tag-side Lambda_c mass constraint +
# the pi0 resonance constraint; loose chi2 cut (tight cuts applied in ROOT).
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")     # tag-side Lambda_c mass
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)    # pi0 resonance
  f.chi2_cut 200
end

# ---- Inexpressible BOSS-side procedures (carried to the systematics layer) ----
alg
  .note(:signal_side_selection,
        "recoil Lambda_c+ built from the tag's unused tracks/showers: charged " \
        "tracks must satisfy |cos(theta)|<0.93, |Vz|<10 cm, Vr<1 cm (K_S0 and " \
        "Lambda daughters exempt); photons require E>25 MeV in the barrel " \
        "(|cos(theta)|<0.8) or E>50 MeV in the endcap (0.84<|cos(theta)|<0.92) " \
        "with EMC time in (0,700) ns; pi0 candidates require 0.115<M(gammagamma)" \
        "<0.150 GeV and are mass-constrained to the nominal pi0; K_S0 and Lambda " \
        "come from oppositely charged pairs with DCA along the beam within " \
        "+-20 cm, vertex chi2<100, a decay vertex at least twice its resolution " \
        "from the IP, and mass windows 0.487-0.511 GeV (K_S0) and 1.111-1.121 GeV " \
        "(Lambda); Sigma0, Sigma+ and omega use 1.179-1.203, 1.176-1.200 and " \
        "0.760-0.800 GeV")
  .note(:pid_correction_method,
        "proton PID uses L(p)>L(K) and L(p)>L(pi); K/pi separation uses " \
        "L(K)>L(pi) or L(pi)>L(K); applied to the signal-side tracks")
  .note(:background_veto,
        "Sigma0, Sigma+ and omega candidates vetoed against Lambda, Sigma+ and " \
        "K_S0 backgrounds; charge-conjugate modes are implicit with all charges " \
        "reversed")
  .note(:fit_resonance_constraints,
        "besides the tag-side Lambda_c mass and pi0 constraints, the relevant " \
        "K_S0, Lambda, Sigma0 and Sigma+ resonances are constrained per mode; " \
        "the two modes with two pi0 use only the four-momentum and tag-side " \
        "Lambda_c mass constraints")
  .note(:deltaE_mbc_window,
        "mode-dependent DeltaE windows of (-0.020,0.020), (-0.030,0.020) or " \
        "(-0.050,0.030) GeV and an M_BC signal region of 2.276-2.300 GeV are " \
        "applied at ROOT level; single- and double-tag yields are extracted " \
        "from M_BC fits for the global branching-fraction fit")

# Render the tag spec (apply takes NO Selection argument) and run on datasets.
alg.apply

root_files = alg.execute_on([data_4599, incMC_4599] + exMCs)