# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
data_4178  = DatasetManager.real_data.find("703_4180")     # 4.178 GeV real data (3.19 fb^-1)
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")  # matching inclusive MC

# Decay card for the signal process
# e+e- -> Ds*+ Ds- -> gamma Ds+ Ds-, with Ds+ -> K+K-pi+
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Ds+ Ds*- PHSP;
  Enddecay

  Decay Ds*-
  1.0000 gamma Ds- PHSP;
  Enddecay

  Decay Ds+
  1.0000 K+ K- pi+ PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC: 500k events (psi(4260) -> Ds+ Ds*-, Ds*- -> gamma Ds-, Ds+ -> K+K-pi+)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4178_DsDsstar_gamma_KKpi"
  config.related_dataset = data_4178
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag analysis (BOSS) ###
alg_name = "DsDsStarDTag"
ds_dtag = TagAnalysis.new(alg_name)
ds_dtag.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 4.178]})
       .with_decay_card(decay_card_signal)

# Tag side: Ds- reconstructed through eight hadronic modes
ds_dtag.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKsK, :DstoKsKPiPi, :DstoKPiPi,
          :DstoKsKKPi, :DstoPiPiPi, :DstoPiPiPiEta, :DstoKKPiPi0
  t.charm -1
end

# Signal side: Ds+ reconstructed in K+K-pi+, best candidate ranked by invariant mass
ds_dtag.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm 1
  t.rank_by :inv
end

# Exactly one photon from Ds* -> Ds gamma; minimum opening angle 10 degrees
ds_dtag.signal_side do |s|
  s.photons 1..1
  s.min_photon_angle 10.0
end

# 6C kinematic fit: 4-momentum conservation + both Ds invariant masses at nominal Ds mass
ds_dtag.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:Ds)
  f.chi2_cut 200
end

# Background suppression applied on stored (pre-fit) quantities
ds_dtag.note(:background_veto,
  "soft pi+/pi- and pi0 with momentum below 0.1 GeV/c vetoed; K_S (0.487-0.511 GeV/c2), " \
  "pi0 (0.115-0.150 GeV/c2, chi2<30) and eta (0.490-0.580 GeV/c2, chi2<30) mass windows applied; " \
  "secondary-vertex significance L/sigma_L > 2 required")

# Selection windows on stored Ds / recoil observables (windowed in ROOT, store-not-cut)
ds_dtag.note(:efficiency_curve,
  "tag Ds mass window mode-dependent in 1.940-1.996 GeV/c2; signal Ds window 1.950-1.986 GeV/c2; " \
  "recoil mass required within 2.051-2.180 GeV/c2; all applied on stored quantities")

ds_dtag.apply

root_files = ds_dtag.execute_on([data_4178, incMC_4178, exMC_signal])