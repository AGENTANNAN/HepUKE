# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset preparation ###
# Seven D_s energy-scan points spanning sqrt(s) = 4.128-4.226 GeV (7.33 fb^-1 total)
ds_scan_data = [
  DatasetManager.real_data.find("705_4130"),   # 4128.5 MeV
  DatasetManager.real_data.find("705_4160"),   # 4157.4 MeV
  DatasetManager.real_data.find("703_4180"),   # 4178.0 MeV (reference point, ECMS = 4.178 GeV)
  DatasetManager.real_data.find("703_4190"),   # 4188.8 MeV
  DatasetManager.real_data.find("703_4200"),   # 4198.9 MeV
  DatasetManager.real_data.find("703_4210"),   # 4209.2 MeV
  DatasetManager.real_data.find("703_4230"),   # 4226.3 MeV
]

# Inclusive MC available at 4.130 and 4.180 GeV only
ds_incMC = [
  DatasetManager.inclusive_mc.find("705_4130"),
  DatasetManager.inclusive_mc.find("703_4180"),
]

# Decay card for the pure-phase-space signal decay
#   D_s+ -> pi+ pi+ pi- pi0 pi0,  pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
  Decay D_s+
  1.000 pi+ pi+ pi- pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 500k-event exclusive signal MC, one sample per scan energy point
exMCs_signal = DatasetManager.create_exclusive_mc_for(ds_scan_data) do |config|
  config.sample_name   = "DsToPiPiPiPiPi0Pi0"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Double-tag analysis (tag-based) ###
alg_name = "DsToPiPiPiPiPi0Pi0"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.178] })   # reference scan energy
   .with_decay_card(decay_card_signal)

# --- Tag side: hadronic D_s- (charm -1) reconstructed in any of the three modes ---
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKKPiPi0   # K_S0 K-, K+ K- pi-, K+ K- pi- pi0
  t.charm -1
end

# --- Signal side: D_s+ decay built from the tracks/showers left over by the tag ---
alg.signal_side do |s|
  s.charged(pip: 2, pim: 1)   # two pi+ and one pi-
  s.require_charge 1          # net charge +1
  s.photons 4                 # two pi0 -> four photons
  s.min_photon_energy 0.025   # each photon energy above 25 MeV
end

# --- Kinematic fit: four-momentum conservation + tag D_s mass constraint ---
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)  # tag D_s -> nominal D_s mass
  f.chi2_cut 200
end

# --- BOSS-side procedures that the tag-DSL vocabulary cannot express directly ---
alg.note(:extended_kinematic_fits,
         "the 4C fit (four-momentum conservation + tag D_s mass constraint, chi2 < 200, " \
         "minimum-chi2 candidate kept per event) is extended to a 6C fit that additionally " \
         "constrains the mass of (tag D_s- + D_s* transition gamma) to the nominal D_s* mass, " \
         "and to a 7C fit that additionally constrains the signal D_s+ mass to the known D_s " \
         "mass; the alternative chi2 values are stored for the ROOT-level hypothesis comparison")
     .note(:transition_photon,
         "D_s* -> gamma D_s transition photon required to have lab energy below 0.2 GeV and a " \
         "recoil mass against the tag D_s- inside [1.95, 2.00] GeV/c^2")
     .note(:background_veto,
         "veto eta -> pi+pi-pi0 (any pi+pi-pi0 invariant mass in 0.49-0.58 GeV/c^2); " \
         "veto K_S0 -> pi0pi0 (M(pi0pi0) outside 0.487-0.511 GeV/c^2); " \
         "veto K_S0 -> pi+pi- (secondary-vertex flight-length significance L/sigma_L > 2); " \
         "veto misidentified open-charm D0/D0bar " \
         "(|M(K-pi+pi0)-M_D0| < 30 MeV/c^2 and |M(K+pi+pi-pi-)-M_D0bar| < 30 MeV/c^2 " \
         "simultaneously)")
     .note(:competing_hypothesis_veto,
         "reject events in which the rho+ eta'(-> pi+pi-gamma) background-hypothesis chi2 is " \
         "smaller than the signal-hypothesis chi2")

alg.apply
alg.execute_on(ds_scan_data + ds_incMC + exMCs_signal)