# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# Seven energy points covering 4.600–4.699 GeV (4.5 fb^-1 in total)
data_points = [
  DatasetManager.real_data.find("703_4600"),   # 4.600 GeV (BOSS 7.0.3)
  DatasetManager.real_data.find("706_4610"),   # 4.610 GeV (BOSS 7.0.6)
  DatasetManager.real_data.find("706_4620"),   # 4.620 GeV
  DatasetManager.real_data.find("706_4640"),   # 4.640 GeV
  DatasetManager.real_data.find("706_4660"),   # 4.660 GeV
  DatasetManager.real_data.find("706_4680"),   # 4.680 GeV
  DatasetManager.real_data.find("706_4700")    # 4.699 GeV
]
# Corresponding inclusive MC samples at the same seven energy points
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700")
]

# Decay card for the signal process e+e- -> Lambda_c+ anti-Lambda_c-,
# Lambda_c+ -> n e+ nu_e (Cabibbo-suppressed beta decay)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 n0 e+ nu_e PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC for Lambda_c+ -> n e+ nu_e, 100k events per energy point
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "LcToNeNu_exclusive_mc"  # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) — tag-based double-tag analysis ###
alg_name = "LcToNeNuTag"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.600]})  # nominal ECMS; the per-run measured beam energy is used in the fit
   .note(:pid_correction_method, "positron identification uses likelihoods with L(e) > 0.001, "
         "L(e)/(L(e)+L(pi)+L(K)) > 0.8 and EMC energy over MDC momentum > 0.5; the tag-based "
         "signal-side lepton PID applies fixed SimplePIDSvc thresholds and cannot be re-tuned in DSL")
   .note(:neutron_pid, "the neutron is identified only through its EMC shower pattern by a "
         "ParticleNet graph-neural-network classifier trained and calibrated on J/psi and Lambda "
         "control samples; no DSL primitive exists for this classifier")
   .with_decay_card(decay_card_signal)

# Tag side: hadronic anti-Lambda_c- (charge -1) from the pre-stored DTag candidates
alg.tag_side(:Lambdac) do |t|
  t.mode_group(:hadronic)   # the 10 hadronic anti-Lambda_c- tag modes
  t.charm(-1)               # tag the anti-Lambda_c- (charm = -1)
end

# Signal side: everything the tag did not use
alg.signal_side do |s|
  s.charged(ep: 1)      # exactly one positron
  s.require_charge(1)   # signal-side net charge = +1
  s.photons 1           # at least one shower (EMC shower pattern of the neutron)
  s.missing :nu_e       # one missing massless nu_e (semileptonic beta decay)
end

# 4C kinematic fit: tag + e+ + nu_e constrained to the measured CMS four-momentum
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# Validate/render the tag specification (takes no Selection argument)
alg.apply

# Execute the same chain on all datasets (7 data points, their inclusive MC, and the signal MC)
root_files = alg.execute_on(data_points + incMC_points + exMC_signal)