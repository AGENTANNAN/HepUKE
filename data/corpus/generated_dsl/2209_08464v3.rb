# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# Seven threshold points of the e+e- -> Lambda_c+ Lambda_c- scan (4.600 - 4.699 GeV)
sample_names = %w[703_4600 706_4610 706_4620 706_4640 706_4660 706_4680 706_4700]

data_points  = sample_names.map { |n| DatasetManager.real_data.find(n) }      # matching real data
incmc_points = sample_names.map { |n| DatasetManager.inclusive_mc.find(n) }  # matching inclusive MC

# Decay card for the signal process (EvtGen format, EvtGen particle names),
# e+e- -> Lambda_c+ Lambda_c- with the signal-side cascade Lambda_c+ -> Lambda pi+ pi0
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 Lambda0 pi+ pi0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC (1M events per energy point, phase space), same selection chain for all points
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_LcLc_LcToLpipi0"
  config.events        = 1_000_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Tag-based event selection (BOSS) ###
alg_name = "LcSTTag"
my_Algorithm = TagAnalysis.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({ "ECMS" => [:double, 4.600] })   # beam energy (per-run value from DB at fit time)
            .with_decay_card(decay_card_signal)

# Tag side: single tag of a Lambda_c through the hadronic mode Lambda_cP -> K pi p, tagged charge pinned to -1
my_Algorithm.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP
  t.charm -1
end

# Signal side: Lambda_c+ -> Lambda pi+ pi0 with Lambda -> p pi- and pi0 -> gamma gamma
my_Algorithm.signal_side do |s|
  s.photons 2                                   # two photons from pi0 -> gamma gamma
  s.charged(prp: 1, pim: 1, pip: 1)             # one p (from Lambda), one pi- (from Lambda), one pi+
  s.require_charge 1                            # net charge of the signal side is +1
  s.min_photon_angle 10.0                       # min angle of signal photons to charged tracks
  s.min_photon_energy 0.025                     # shower energy floor (barrel threshold)
end

# Kinematic fit: 4C four-momentum conservation + Lambda mass + pi0 mass + tag Lambda_c+ mass
# (the paper's 3C fit emitted here as the equivalent 7C fit), stored chi2 with a loose cut
my_Algorithm.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:prp, :pim).constrain_to_nominal_mass_of(:Lambda)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

# BOSS-side procedures that have no DSL construct — preserved for the systematic-uncertainty stage
my_Algorithm
  .note(:track_quality, "signal-side charged tracks use the standard quality cuts |cos(theta)| < 0.93, |Vz| < 10 cm (20 cm for the Lambda daughters) and Vr < 1 cm")
  .note(:photon_selection, "signal photons require E > 25 MeV in the EMC barrel and E > 50 MeV in the endcap, with EMC time within 700 ns")
  .note(:mass_window, "signal-side pi0(gamma gamma) invariant mass required in [0.115, 0.150] GeV and Lambda(p pi-) in [1.111, 1.121] GeV before the kinematic fit")
  .note(:vertex_fit, "signal-side Lambda -> p pi- secondary vertex fit required with vertex-fit chi2 < 100 and decay length significance > 2 sigma")
  .note(:background_veto, "Sigma0 veto: events whose signal-side Lambda_c+ candidate also satisfies the Lambda_c+ -> Sigma0 pi+ (Sigma0 -> Lambda gamma) selection are removed")
  .note(:candidate_selection, "signal-side Lambda_c+ candidate chosen with minimal |DeltaE| restricted to [-0.03, 0.02] GeV; tag-side mBC and DeltaE are stored unconditionally and windowed in the ROOT stage (partial-wave analysis of Lambda_c+ -> Lambda pi+ pi0, with rho(770)+, Sigma(1385)+, Sigma(1385)0, also left to ROOT)")

# Render the tag spec (takes no Selection argument) and run on data, inclusive MC and signal MC
my_Algorithm.apply
root_files = my_Algorithm.execute_on(data_points + incmc_points + exMCs)