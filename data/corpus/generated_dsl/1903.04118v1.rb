# Core DSL classes and dependencies will be loaded automatically at execution
# =====================================================================
# Double-tag analysis of D_s+ -> pi+ pi0 eta  (pi0->gg, eta->gg)
# sqrt(s) = 4.178 GeV ; tag side : D_s- (charm = -1) in seven ST modes
# =====================================================================

### Dataset preparation ###
data_4178   = DatasetManager.real_data.find("703_4180")      # 4.178 GeV real data set (3.19 fb^-1)
incMC_4178  = DatasetManager.inclusive_mc.find("703_4180")   # corresponding inclusive MC sample

# Decay card for the signal process, generated with the ConExc generator:
#   e+e- -> D_s*- D_s+ ,  D_s*- -> gamma D_s-  (direct photon, tag side)
#                         D_s+  -> pi+ pi0 eta (signal side), pi0/eta -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.0000 ConExc -1 vhdr;
    Enddecay

    Decay vhdr
    1.0000 D_s*- D_s+ PHSP;
    Enddecay

    Decay D_s*-
    1.0000 gamma D_s- PHSP;
    Enddecay

    Decay D_s+
    1.0000 pi+ pi0 eta PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal (single energy point -> one ExclusiveMC)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4178_ds_pipipi0eta_dt"
  config.related_dataset = data_4178
  config.events          = 100_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS, tag-based) ###
alg_name = "DsToPiPi0EtaDT"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.178]})   # sqrt(s) = 4.178 GeV
   .with_decay_card(decay_card_signal)

# ---- Tag side: D_s- reconstructed from the seven single-tag modes, charm = -1
# (one tag_side call + signal_side content -> the ST tag plus the DT signal side)
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK,        # K_S K-
          :DstoKKPi,       # K+ K- pi-
          :DstoKsKPi0,     # K_S K- pi0
          :DstoKKPiPi0,    # K+ K- pi- pi0
          :DstoKsKPiPi,    # K_S K+ pi- pi-
          :DstoPiEta,      # pi- eta
          :DstoPiEtaPrime  # pi- eta'
  t.charm(-1)              # tag side pinned to charm = -1 (D_s-)
end

# ---- Signal side: D_s+ -> pi+ pi0 eta with pi0 -> gamma gamma and eta -> gamma gamma
#      (four photons), exactly one pi+ and net charge +1
alg.signal_side do |s|
  s.photons 4            # four good showers from pi0 -> gg and eta -> gg
  s.charged(pip: 1)      # exactly one pi+
  s.require_charge 1     # net charge of the signal side = +1
end

# ---- 7C kinematic fit: 4-momentum conservation + m(pi0) + m(eta) + m(D_s) (tag side)
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # pi0 -> gamma gamma
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)  # eta -> gamma gamma
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)            # tagged D_s mass
  f.chi2_cut 1000
  f.store_fitted_momenta    # save the fitted four-momenta (incl. those of the nominal 7C fit)
end

# ---- BOSS-side procedures that the tag-fit surface cannot express
alg.note(:gamma_dir_selection,
         "the direct photon from D_s*- -> gamma D_s- (gamma_dir) is not one of the declared " \
         "signal-side photons; in the generated selection the remaining good shower of the event " \
         "is assigned to gamma_dir by scanning the photon combinations and keeping the one that " \
         "gives the smallest 7C chi2")
   .note(:background_veto,
         "7CA variant: a second, non-nominal 7C fit hypothesis (4-momentum conservation, m(pi0), " \
         "m(eta) plus the constraint of gamma_dir together with the tagged D_s+- to the D_s* mass) " \
         "is evaluated and its chi2 stored per event for the background suppression performed in ROOT")

alg.apply    # takes no Selection argument for a tag-based analysis
root_files = alg.execute_on([data_4178, incMC_4178, exMC_signal])