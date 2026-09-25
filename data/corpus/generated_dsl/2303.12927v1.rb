# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Real data: eight c.m. energy points between 4.128 and 4.226 GeV (7.33 fb^-1 in total)
data_4128 = DatasetManager.real_data.find("705_4130")   # 4.128 GeV
data_4157 = DatasetManager.real_data.find("705_4160")   # 4.157 GeV
data_4178 = DatasetManager.real_data.find("703_4180")   # 4.178 GeV
data_4189 = DatasetManager.real_data.find("703_4190")   # 4.189 GeV
data_4199 = DatasetManager.real_data.find("703_4200")   # 4.199 GeV
data_4209 = DatasetManager.real_data.find("703_4210")   # 4.209 GeV
data_4219 = DatasetManager.real_data.find("703_4220")   # 4.219 GeV
data_4226 = DatasetManager.real_data.find("703_4230")   # 4.226 GeV
data_points = [data_4128, data_4157, data_4178, data_4189,
               data_4199, data_4209, data_4219, data_4226]

# Corresponding inclusive MC samples for each energy point
incMC_points = [
  DatasetManager.inclusive_mc.find("705_4130"),
  DatasetManager.inclusive_mc.find("705_4160"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4210"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230")
]

# Decay card for the signal process (EvtGen format, EvtGen particle names):
# e+e- -> D_s*+ D_s-, D_s*+ -> gamma D_s+, D_s+ -> f_0(980) e+ nu_e, f_0(980) -> pi+ pi-,
# the tag side D_s- -> K+ K- pi-
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000  D_s*+  D_s-     PHSP;
  Enddecay

  Decay D_s*+
  1.0000  gamma  D_s+     PHSP;
  Enddecay

  Decay D_s+
  1.0000  f_0  e+  nu_e   PHSP;
  Enddecay

  Decay f_0
  1.0000  pi+  pi-       PHSP;
  Enddecay

  Decay D_s-
  1.0000  K+  K-  pi-     PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive signal MC: 500k events, distributed over the eight energy points by
# luminosity weight (passing the Array of datasets selects the multi-round form)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_dsstar_ds_f0enu"
  config.related_dataset = data_points
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) -- tag-based analysis (D_s tag) ###
# The event is tagged by a hadronically reconstructed D_s- (from the pre-stored
# EvtRecDTag collection); the recoiling D_s*+ -> gamma D_s+ -> gamma f_0(980) e+ nu_e
# is the signal side.  Pattern: one tag side + a signal side carrying a massless
# missing nu_e  =>  single tag + missing particle.
alg_name = "DsTagF0ENu"
ds_analysis = TagAnalysis.new(alg_name)
ds_analysis.set_header(["#{alg_name}Alg/#{alg_name}.h"])
           .set_constant({"ECMS" => [:double, 4.178]})  # representative c.m. energy;
                                                        # the fit builds the per-run
                                                        # measured CMS four-vector itself
           .with_decay_card(decay_card_signal)
           # Inexpressible BOSS-side facts captured for the systematics layer:
           .note(:tag_mode_unavailable,
                 "the tag modes K_S0 K- pi0, K- pi+ pi- and K_S0 K- pi+ pi- are not " \
                 "available in the tagging algorithm; those tag modes are reconstructed " \
                 "and combined downstream in the ROOT analysis")
           .note(:tag_candidate_ranking,
                 "when several tag candidates survive in one event the one whose recoil " \
                 "mass M_rec is closest to the nominal D_s* mass is retained; this " \
                 "ranking is not expressible in the tag DSL and is applied together with " \
                 "the M_rec^2 / M_miss^2 windows")
           .note(:efficiency_curve,
                 "the f_0(980) lineshape is described by the Flatte formula (BESII " \
                 "parameters); the signal efficiency is therefore non-flat over the " \
                 "pi+pi- mass and a weighted efficiency of (35.44 +/- 0.07)% is used")

# --- Tag side: hadronic D_s- , charm = -1 --------------------------------
ds_analysis.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,        # K+ K- pi-
          :DstoKKPiPi0,     # K+ K- pi- pi0
          :DstoPiPiPi,      # pi+ pi- pi-
          :DstoKsK,         # K_S0 K-
          :DstoKsKPiPi,     # K_S0 K+ pi- pi-
          :DstoPiEta,       # pi- eta (eta -> gamma gamma)
          :DstoRhoEta,      # rho- eta (rho- -> pi- pi0)
          :DstoPiEtaPrime   # pi- eta' (eta' -> pi+ pi- eta and eta' -> gamma rho0)
  t.charm -1              # pin the tagged side to the D_s-
  # No window declared: tag mBC / deltaE and the missing mass are stored
  # unconditionally (store-not-cut) and windowed later in the ROOT analysis.
end

# --- Signal side: the D_s*+ -> gamma D_s+ -> gamma f_0(980) e+ nu_e ---------
ds_analysis.signal_side do |s|
  s.photons 1                        # the D_s*+ -> gamma D_s+ transition photon
  s.charged(pip: 1, pim: 1, ep: 1)   # f_0(980) -> pi+ pi- , electron via the lepton key
  s.missing :nu_e                    # massless missing neutrino
end

# --- Kinematic fit over the derived participants ---------------------------
ds_analysis.fit do |f|
  f.constrain_four_momentum                                             # 4C constraint
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:f_0)    # f_0(980) mass
  f.invariant_mass_of(:pip, :pim).between(0.6, 1.6)                     # M(pi+pi-) window
  f.chi2_cut 200                                                        # loose chi2 cut
end

# Validate, render and run.  apply takes NO Selection argument for a tag analysis.
ds_analysis.apply
root_files = ds_analysis.execute_on(data_points + incMC_points + exMC_signal)