# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# sqrt(s) = 4.178 GeV -> BOSS sample name convention: [BOSS]_[CMS energy in MeV]
ds_data  = DatasetManager.real_data.find("703_4180")       # 7.33 fb^-1 real data at the Ds threshold (4.178 GeV)
ds_incMC = DatasetManager.inclusive_mc.find("703_4180")    # Corresponding inclusive MC sample

# Decay card for the signal process (EvtGen syntax, EvtGen particle names):
#   Ds*+ -> gamma Ds+ (transition photon), Ds+ -> pi+ pi0 pi0 eta (phase space),
#   with pi0 -> gamma gamma and eta -> gamma gamma.
decay_card_signal = <<~DECAYCARD
  Decay Ds*+
  1.0000 gamma Ds+
  PHSP;
  Enddecay

  Decay Ds+
  1.0000 pi+ pi0 pi0 eta
  PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma
  PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma
  PHSP;
  Enddecay

  End
DECAYCARD

# 500k exclusive signal MC events for Ds*+ -> gamma Ds+, Ds+ -> pi+ pi0 pi0 eta
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_dsstar_ds_pipipi0pi0eta"
  config.related_dataset = ds_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based analysis (double-tag: tag Ds- , signal Ds+ -> pi+ pi0 pi0 eta) ###
alg_name = "DsDoubleTagPiPi0Pi0Eta"
ds_alg = TagAnalysis.new(alg_name)                              # TagAnalysis < Algorithm
ds_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
      .set_constant({"ECMS" => [:double, 4.178]})               # CMS energy 4.178 GeV
      .with_decay_card(decay_card_signal)

# ---- Tag side: recoiling Ds- through nine hadronic decay modes (single tag reconstruction) ----
ds_alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK,          # K_S K
          :DstoKKPi,         # K K pi
          :DstoKsKPi0,       # K_S K pi0
          :DstoKKPiPi0,      # K K pi pi0
          :DstoPiPiPi,       # pi pi pi
          :DstoPiEta,        # pi eta
          :DstoPiPi0Eta,     # pi pi0 eta
          :DstoPiEtaPrime,   # pi eta'
          :DstoKPiPi         # K pi pi
  t.charm -1                 # pin the tagged side to Ds-
end

# ---- Signal side: Ds+ -> pi+ pi0 pi0 eta (pi0/eta -> gamma gamma) plus transition photon ----
ds_alg.signal_side do |s|
  s.charged(pip: 1)          # exactly one charged pi+ track
  s.require_charge(1)        # net charge +1
  s.photons 7                # seven photons: 4 (two pi0) + 2 (eta) + 1 transition photon
  s.min_photon_angle 10.0    # minimum photon angle to the charged tracks (degrees)
end

# ---- 9C kinematic fit (4C + 2 pi0 + eta + tag Ds- mass + Ds*+ mass) ----
ds_alg.fit do |f|
  f.constrain_four_momentum                                                  # 4C four-momentum conservation
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)     # first  pi0 -> gamma gamma
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)     # second pi0 -> gamma gamma
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)     # eta -> gamma gamma
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)               # tag Ds- mass
  f.invariant_mass_of(:tag1, :gamma).constrain_to_nominal_mass_of(:Ds_star)  # Ds*+ mass (Ds gamma)
  f.chi2_cut 200                                                             # accept smallest chi^2, chi^2 < 200
end

# ---- BOSS-side procedures that the tag DSL cannot express ----
ds_alg.note(:transition_photon_energy,
            "the transition photon from Ds*+ -> gamma Ds+ is required to have laboratory
             energy below 0.18 GeV; this upper bound on one signal-side photon is applied
             as an event-selection cut on the signal photon list (only a photon-energy
             floor, min_photon_energy, is expressible in the DSL)")
ds_alg.note(:recoil_mass_window,
            "recoil mass against the Ds*+ (signal Ds+ plus transition photon) required to be
             within [1.93, 1.98] GeV/c^2; stored and windowed at BOSS level, since a
             tag-group invariant-mass window (.between over :tag1) is not expressible in the
             tag-fit DSL")
ds_alg.note(:background_veto,
            "photon pairs with invariant mass in [0.10, 0.14] GeV/c^2 are vetoed to suppress
             fake pi0 candidates; applied on the signal-side photon list before the 9C fit")
ds_alg.note(:ds_star_mass_constraint,
            "the Ds*+ mass constraint in the 9C fit is imposed from either the Ds+ gamma or the
             Ds- gamma combination, keeping the combination with the smallest chi^2")

# Validate + render the tag spec (apply takes NO Selection argument)
ds_alg.apply

# Execute on real data, inclusive MC and the exclusive signal MC
root_files = ds_alg.execute_on([ds_data, ds_incMC, exMC_signal])