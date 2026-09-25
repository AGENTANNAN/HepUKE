# ============================================================================
# BOSS-side DSL for  e+e- -> Ds*+ Ds-  with  Ds*+ -> gamma Ds+,
# Ds+ -> pi+ pi+ pi- eta,  eta -> gamma gamma,  using a Ds- double-tag
# (eight hadronic tag modes).  Energy points:
# 4.180, 4.190, 4.200, 4.210, 4.220, 4.230 GeV.
# ============================================================================

### Dataset description ###
# Six real-data points (XYZ scan) with their matching inclusive MC
data_points = [
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230")
]

incMC_points = [
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230")
]

# Decay card for the signal process (EvtGen format)
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000 Ds*+ Ds- PHSP;
  Enddecay

  Decay Ds*+
  1.000 gamma Ds+ PHSP;
  Enddecay

  Decay Ds+
  1.000 pi+ pi+ pi- eta PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 500k-event exclusive signal MC at each of the six scan points
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_dsstar_ds_signal"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) — tag-based double-tag ###
alg = TagAnalysis.new("DsStarDsDoubleTag")
alg.set_header(["DsStarDsDoubleTagAlg/DsStarDsDoubleTag.h"])
   .set_constant({"ECMS" => [:double, 4.180]})  # nominal scan energy; the fit uses the per-run measured CMS
   .with_decay_card(decay_card_signal)

# Tag side: Ds- reconstructed in eight hadronic modes, opposite charm to the signal Ds+
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK,        # K_S K
          :DstoKKPi,       # K K pi
          :DstoKsKPi0,     # K_S K pi0
          :DstoKKPiPi0,    # K K pi pi0
          :DstoKsKPimPip,  # K_S K- pi+ pi-
          :DstoKsKpPimPim, # K_S K+ pi- pi-
          :DstoPiPiPiEta,  # pi+ pi+ pi- eta
          :DstoKPiPi       # K pi pi
  t.charm -1
end

# Signal side: Ds+ built from the tracks/showers left over by the tag
alg.signal_side do |s|
  s.photons 2                 # exactly two photons (the eta -> gamma gamma pair)
  s.min_photon_angle 10.0     # minimum photon angle to charged tracks (degrees)
  s.charged(pip: 2, pim: 1)   # two pi+ and one pi-
  s.require_charge(1)
end

# 7C kinematic fit: 4-momentum conservation + eta mass + tag Ds- mass
# (the Ds*+ mass constraint involving the transition photon is applied later in ROOT)
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.chi2_cut 200   # minimum-chi2 combination kept automatically when multiple combinations (~15%) arise
end

# Inexpressible BOSS-side procedures (stored, cut downstream)
alg.note(:background_veto,
         "signal-side background vetoes, not expressible in the tag DSL: " \
         "(1) K_S0 veto — reject secondary-vertex M(pi+pi-) in [0.487, 0.511] GeV/c^2 with L/sigma_L > 2; " \
         "(2) eta' veto — reject M(pi+ pi- eta) < 1 GeV/c^2; " \
         "(3) pi0 cross-feed veto — reject M(gamma_eta gamma) in [0.115, 0.150] GeV/c^2")
   .note(:tag_mode_veto,
         "tag mode K- pi+ pi- (DstoKPiPi): exclude the tag-side M(pi+pi-) in [0.487, 0.511] GeV/c^2 " \
         "to remove the Ds- -> K_S0 K- overlap")
   .note(:bdt_selection,
         "multivariate BDT (signal MC vs inclusive MC) to suppress combinatorial/continuum backgrounds " \
         "and reach a tag purity above 85%; applied as a ROOT-level selection")

alg.apply
alg.execute_on(data_points + incMC_points + exMCs_signal)