### Dataset description ###
# 4.178 GeV data sample (BOSS 703): real data + inclusive MC
data_4180  = DatasetManager.real_data.find("703_4180")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")

# Decay card for the signal process (EvtGen syntax).
# e+e- -> Ds+ Ds*- (KKMC convention top mother: psi(4260)),
# Ds*- -> Ds- gamma (soft photon), Ds- -> K+ K- pi-, Ds+ -> gamma e+ nu_e.
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000  D_s+  D_s*-  PHSP;
  Enddecay

  Decay D_s*-
  1.000  D_s-  gamma  PHSP;
  Enddecay

  Decay D_s-
  1.000  K+  K-  pi-  PHSP;
  Enddecay

  Decay D_s+
  1.000  gamma  e+  nu_e  PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive (signal) MC: 100k events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4180_ds_gamma_e_nu"
  config.related_dataset = data_4180
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) — modified double-tag via TagAnalysis ###
alg_name = "DsToGammaENu"
alg = TagAnalysis.new(alg_name)                       # TagAnalysis < Algorithm; no Selection object
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])      # generate header for this algorithm
   .set_constant({"ECMS" => [:double, 4.178]})        # 4.178 GeV

# Tag side: single-tag a Ds- (charm = -1) in the 13 hadronic modes.
# eta/eta' and the Ds*- -> Ds- soft photon are reconstructed internally by the tag algorithm.
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,        # K+ K- pi-
          :DstoKKPiPi0,     # K+ K- pi- pi0
          :DstoKsK,         # K_S K
          :DstoPiEta,       # pi eta
          :DstoPiPiPi,      # pi pi pi
          :DstoKPiPi,       # K pi pi
          :DstoKsKmPiPi,    # K_S K- pi pi
          :DstoKsKpPiPi,    # K_S K+ pi pi
          :DstoKsKsPi,      # K_S K_S pi
          :DstoKsKPi0,      # K_S K pi0
          :DstoKsPi,        # K_S pi
          :DstoPiPiPiEta,   # pi pi pi eta
          :DstoPiPi0Eta     # pi pi0 eta
  t.charm -1
end

# Signal side: the tag's unused tracks / showers (Ds+ -> gamma e+ nu_e).
alg.signal_side do |s|
  s.photons 1                 # require (at least) one signal photon
  s.min_photon_energy 0.025   # photon energy > 25 MeV
  s.min_photon_angle 10.0     # photon angle to nearest charged track > 10 deg
  s.charged(ep: 1)            # exactly one e+
  s.require_charge 1          # net charge = +1
  s.missing :nu_e             # massless missing neutrino (semileptonic)
end

# Kinematic fit: four-momentum constraint, chi2 < 200.
# The four Ds+ production hypotheses are resolved later, at the ROOT level.
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.apply                                                  # no Selection argument for a TagAnalysis
alg.execute_on([data_4180, incMC_4180, exMC_signal])       # produce ROOT files