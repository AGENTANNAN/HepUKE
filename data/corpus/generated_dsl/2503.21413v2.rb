# =============================================================================
# BOSS-side specification for the search for
#   e+e- -> anti-Lambda_c- Lambda_c(2595)+ / Lambda_c(2625)+,
#   Lambda_c*+ -> Lambda_c+ pi0 pi0   (only the two pi0 are reconstructed,
#                                      the Lambda_c+ itself is NOT)
# Tag side : anti-Lambda_c- -> anti-p K+ pi-   (PHSP),
#                               anti-p K_S0     (K_S0 -> pi+ pi-),
#                               anti-Lambda pi- (Lambda -> p pi-)
# Both charge-conjugate channels are covered by scanning both tag charges.
# =============================================================================

### ------------------------- Dataset preparation ------------------------- ###
data_4918  = DatasetManager.real_data.find("707_4914")      # sqrt(s) = 4.918 GeV
incMC_4918 = DatasetManager.inclusive_mc.find("707_4914")   # inclusive MC @ 4.918 GeV
data_4951  = DatasetManager.real_data.find("707_4946")      # sqrt(s) = 4.951 GeV
incMC_4951 = DatasetManager.inclusive_mc.find("707_4946")   # inclusive MC @ 4.951 GeV

# Signal decay card.  The tag baryon decays through the phase-space mode
# anti-Lambda_c- -> anti-p K+ pi-; the signal side Lambda_c+ pi0 pi0 is kept
# generic so that the SAME card covers both the Lambda_c(2595)+ and the
# Lambda_c(2625)+ searches (the two states are separated in the recoil-mass
# distribution M_recoil(anti-Lambda_c-) in the ROOT analysis).
# psi(4260) is used as the top mother (BESIII convention for the KKMC generator).
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 anti-Lambda_c- Lambda_c+ pi0 pi0 PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the peaking background Lambda_c*+ -> Lambda_c+ pi+ pi-
# (subtracted from the data using this dedicated MC sample).
decay_card_bkg = <<~DECAYCARD
    Decay psi(4260)
    1.0000 anti-Lambda_c- Lambda_c+ pi+ pi- PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Signal exclusive MC : the same signal MC is run over both energy points.
exMCs_signal = DatasetManager.create_exclusive_mc_for([data_4918, data_4951]) do |config|
  config.sample_name   = "exmc_lambdacbar_lambdacstar_pi0pi0"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# Peaking-background exclusive MC, also over both energy points.
exMCs_bkg = DatasetManager.create_exclusive_mc_for([data_4918, data_4951]) do |config|
  config.sample_name   = "exmc_lambdacbar_lambdacstar_pipi"
  config.events        = 100_000
  config.decay_card    = decay_card_bkg
  config.cross_section = :default
end

### --------------------------- Event selection --------------------------- ###
# The tag based analysis is carried out separately at each energy point
# (the declared ECMS is used by the 4C fit and for the beam-energy spread).
# ============================== 4.918 GeV ============================== #
alg_name_4918 = "LambdaCStarTag4918"
alg_4918 = TagAnalysis.new(alg_name_4918)
alg_4918.set_header(["#{alg_name_4918}Alg/#{alg_name_4918}.h"])
        .set_constant({ "ECMS" => [:double, 4.918] })   # sqrt(s) = 4.918 GeV
        .note(:tag_mass_windows,
              "tag-side mass windows are stored and applied in ROOT (store-not-cut): " \
              "2.27 < M(anti-Lambda_c-) < 2.30 GeV/c^2, M(K_S0) in (0.487, 0.511) GeV/c^2, " \
              "M(Lambda) in (1.111, 1.121) GeV/c^2; the K_S0 and Lambda windows are " \
              "internal to the DTagAlg tag reconstruction and are not re-applied here.")
        .note(:soft_photon_energy_cut,
              "signal-side showers must satisfy E_gamma < 150 MeV (soft photons from the " \
              "low-Q-value Lambda_c(2595/2625)+ -> Lambda_c+ pi0 pi0 decay); a maximum " \
              "shower-energy requirement is not expressible in the DSL (only a minimum " \
              "shower energy is available), so it is applied in the ROOT analysis.")
        .note(:peaking_background_subtraction,
              "the peaking background from Lambda_c*+ -> Lambda_c+ pi+ pi- is subtracted " \
              "using the dedicated exclusive MC sample generated above.")
        .note(:tag_mode_efficiency,
              "tag-mode efficiencies approximately 48% (anti-p K+ pi-), 50% (anti-p K_S0) " \
              "and 38% (anti-Lambda pi-); input for the systematic uncertainty evaluation.")
        .with_decay_card(decay_card_signal)

# Tag side: one tag_side call -> single tag.  charm is left unpinned so that
# both tag charges (anti-Lambda_c- and Lambda_c+) are scanned, covering the
# two charge-conjugate channels of the analysis.
alg_4918.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP,      # anti-p K+ pi-
          :LambdacPtoKsP,       # anti-p K_S0
          :LambdacPtoLambdaPi   # anti-Lambda pi-
end

# Signal side: at least four photons forming the two pi0; the Lambda_c+ is not
# reconstructed and is declared as the missing massive particle.
alg_4918.signal_side do |s|
  s.photons 4
  s.missing :"Lambda_c+"
end

# 4C kinematic fit over the derived participants (tag + signal photons + missing Lambda_c+).
alg_4918.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_4918.apply   # takes NO Selection argument
alg_4918.execute_on([data_4918, incMC_4918] + exMCs_signal + exMCs_bkg)

# ============================== 4.951 GeV ============================== #
alg_name_4951 = "LambdaCStarTag4951"
alg_4951 = TagAnalysis.new(alg_name_4951)
alg_4951.set_header(["#{alg_name_4951}Alg/#{alg_name_4951}.h"])
        .set_constant({ "ECMS" => [:double, 4.951] })   # sqrt(s) = 4.951 GeV
        .note(:tag_mass_windows,
              "tag-side mass windows are stored and applied in ROOT (store-not-cut): " \
              "2.27 < M(anti-Lambda_c-) < 2.30 GeV/c^2, M(K_S0) in (0.487, 0.511) GeV/c^2, " \
              "M(Lambda) in (1.111, 1.121) GeV/c^2; the K_S0 and Lambda windows are " \
              "internal to the DTagAlg tag reconstruction and are not re-applied here.")
        .note(:soft_photon_energy_cut,
              "signal-side showers must satisfy E_gamma < 150 MeV (soft photons from the " \
              "low-Q-value Lambda_c(2595/2625)+ -> Lambda_c+ pi0 pi0 decay); a maximum " \
              "shower-energy requirement is not expressible in the DSL (only a minimum " \
              "shower energy is available), so it is applied in the ROOT analysis.")
        .note(:peaking_background_subtraction,
              "the peaking background from Lambda_c*+ -> Lambda_c+ pi+ pi- is subtracted " \
              "using the dedicated exclusive MC sample generated above.")
        .note(:tag_mode_efficiency,
              "tag-mode efficiencies approximately 48% (anti-p K+ pi-), 50% (anti-p K_S0) " \
              "and 38% (anti-Lambda pi-); input for the systematic uncertainty evaluation.")
        .with_decay_card(decay_card_signal)

alg_4951.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP,
          :LambdacPtoKsP,
          :LambdacPtoLambdaPi
end

alg_4951.signal_side do |s|
  s.photons 4
  s.missing :"Lambda_c+"
end

alg_4951.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_4951.apply
alg_4951.execute_on([data_4951, incMC_4951] + exMCs_signal + exMCs_bkg)