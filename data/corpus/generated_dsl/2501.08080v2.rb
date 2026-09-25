### Dataset preparation ###
# 3.097 GeV J/psi: real data and matching inclusive MC (BOSS 7.0.8 / 3097 MeV)
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Signal decay card: J/psi -> D0 mu+ mu-  (D0 -> K- pi+)
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 D0 mu+ mu- PHSP;
  Enddecay

  Decay D0
  1.0000 K- pi+ PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for the signal channel (100k events); only the K-pi+ D0 mode is generated
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_D0mumu_KPi"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
# The D0 is taken from the pre-stored DTag candidate collection -> TagAnalysis (not Algorithm + Selection)
alg_name = "JpsiD0MuMuTag"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })
   .with_decay_card(decay_card_signal)
   # Signal-side muon ID: custom EMC energy window (v1 lepton PID thresholds are fixed in the tag framework)
   .note(:muon_pid_selection,
         "signal-side muon identification uses an EMC energy window 0.11 < E < 0.25 GeV; " \
         "the tag framework's lepton PID thresholds are fixed in v1, so this custom EMC " \
         "window is applied downstream")
   # Event-level residual-momentum cut
   .note(:missing_momentum_cut,
         "event-level residual momentum |p_miss| < 0.05 GeV/c, computed from the tracks/showers " \
         "left unused by the tag and the signal side")
   # Competing-hypothesis / physics vetoes handled outside the tag framework
   .note(:background_veto,
         "events vetoed against K_S0, four-pion and recoil-mass backgrounds before the final " \
         "M(D0 mu+ mu-) fit")
   # Tag-side D0 quality selections (cannot be a tag_side window: only :deltaE/:mBC are supported)
   .note(:tag_side_d0_selection,
         "tag-side windows applied downstream: M(D0) in [1.84,1.89] GeV/c^2 for D0->K-pi+ and " \
         "D0->K-pi+pi+pi-, [1.80,1.91] GeV/c^2 for D0->K-pi+pi0; tag D0 vertex chi2 < 4.5, 3.1, " \
         "3.7 for the three modes respectively")

# Tag side: single D0 tag reconstructed in three modes from pre-stored candidates
# (both D0 flavours scanned; no tag-side window declared -> store-not-cut)
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
end

# Signal side: exactly one mu+ and one mu-, net charge zero
alg.signal_side do |s|
  s.charged(mup: 1, mum: 1)
  s.require_charge(0)
end

# 5C kinematic fit: 4-momentum conservation (4C) + D0 nominal-mass constraint (1C)
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg.apply
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])