# ===============================================================
# Dataset preparation
# ===============================================================
# Six 703 energy-point real datasets spanning 4.178–4.226 GeV
data_4178 = DatasetManager.real_data.find("703_4180")
data_4188 = DatasetManager.real_data.find("703_4190")
data_4198 = DatasetManager.real_data.find("703_4200")
data_4209 = DatasetManager.real_data.find("703_4210")
data_4218 = DatasetManager.real_data.find("703_4220")
data_4226 = DatasetManager.real_data.find("703_4230")

# Matching inclusive MC samples for the same energy points
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4188 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4198 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4209 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4218 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")

data_points = [data_4178, data_4188, data_4198, data_4209, data_4218, data_4226]
incMCs      = [incMC_4178, incMC_4188, incMC_4198, incMC_4209, incMC_4218, incMC_4226]

# Decay card for e+e- -> Ds*- Ds+ (Ds*- -> gamma Ds-), EvtGen format
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Ds*- Ds+ PHSP;
    Enddecay

    Decay Ds*-
    1.0000 gamma Ds- PHSP;
    Enddecay

    Decay Ds-
    1.0000 K+ K- pi- PHSP;
    Enddecay

    Decay Ds+
    1.0000 pi+ pi0 eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 500k events of e+e- -> Ds*- Ds+, Ds- -> K+K-pi-, Ds+ -> pi+pi0eta'
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4226_DsstDs"
  config.related_dataset = data_4226
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# ===============================================================
# Tag-based event selection (BOSS)
# ===============================================================
tag_alg = TagAnalysis.new("DsTagDsStar")
tag_alg.set_header(["DsTagDsStarAlg/DsTagDsStar.h"])
       .set_constant({ "ECMS" => [:double, 4.226] })   # c.m. energy 4.226 GeV
       .with_decay_card(decay_card_signal)

# Tag side: Ds- reconstructed from pre-stored tag candidates in 12 hadronic
# modes; charm = -1 pins the tagged Ds- side.
tag_alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKsKPi0, :DstoKsKKPiPi, :DstoKsKPiPiPi0,
          :DstoPiPiPiPi0, :DstoPiPi0Eta, :DstoKKPiPi, :DstoKsKPiPi0,
          :DstoPiPiPi, :DstoPiPiPiPi0Eta, :DstoPiEta
  t.charm(-1)
end

# Signal side: everything the tag did not use — Ds+ -> pi+ pi0 eta'
# (eta' -> pi+ pi- eta, eta -> gamma gamma, pi0 -> gamma gamma).
# PID is implicit in the tag-mode hypotheses and the signal charged multiplicities.
tag_alg.signal_side do |s|
  s.photons 4                    # at least four good showers
  s.min_photon_energy 0.025      # E_gamma > 25 MeV
  s.min_photon_angle 10.0        # photon opening angle > 10 degrees
  s.charged(pip: 2, pim: 1)      # exactly two pi+ and one pi-
end

# Kinematic fit: four-momentum conservation, pi0 and eta nominal-mass
# constraints on the two gamma-gamma pairs, tag Ds mass constrained to nominal.
tag_alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.chi2_cut 200
end

tag_alg.apply
root_files = tag_alg.execute_on(data_points + incMCs + [exMC_signal])