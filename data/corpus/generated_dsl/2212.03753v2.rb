### Dataset preparation ###
# Seven XYZ energy points of real data (4.600 - 4.699 GeV), each with its matching inclusive MC.
data_4600 = DatasetManager.real_data.find("703_4600")   # 4.600 GeV
data_4612 = DatasetManager.real_data.find("706_4610")   # 4.612 GeV
data_4628 = DatasetManager.real_data.find("706_4620")   # 4.628 GeV
data_4640 = DatasetManager.real_data.find("706_4640")   # 4.640 GeV
data_4661 = DatasetManager.real_data.find("706_4660")   # 4.661 GeV
data_4682 = DatasetManager.real_data.find("706_4680")   # 4.682 GeV
data_4699 = DatasetManager.real_data.find("706_4700")   # 4.699 GeV

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4612 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4628 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4661 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4682 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4699 = DatasetManager.inclusive_mc.find("706_4700")

data_points  = [data_4600, data_4612, data_4628, data_4640, data_4661, data_4682, data_4699]
incMC_points = [incMC_4600, incMC_4612, incMC_4628, incMC_4640, incMC_4661, incMC_4682, incMC_4699]

# Decay card for the inclusive semileptonic signal Lambda_c+ -> X e+ nu_e (X inclusive, phase space).
# No explicit top mother in the description -> psi(4260) KKMC convention.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c-   PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 X e+ nu_e   PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive signal MC at each of the seven energy points
# (same decay card / cross section / event count, differing only by related dataset)
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Lambdac_Xenu_signal_exMC"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) — tag-and-missing (single tag) ###
alg_name = "LambdacXenuTag"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.640] })

# Tag side: one hadronic Lambda_c anti-charm tag taken from the pre-stored DTag collection.
# The paper's 12 hadronic tag modes are approximated by the hadronic mode group, with the
# p K pi and p K_S modes explicitly included.
alg.tag_side(:Lambdac) do |t|
  t.mode_group :hadronic
  t.modes :LambdacPtoKPiP, :LambdacPtoKsP
  t.charm -1
end

# Signal side: exactly one positron (charge +1) plus one missing massless nu_e.
alg.signal_side do |s|
  s.charged(ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

# 4C kinematic fit: tag Lambda_c- + e+ + nu_e constrained to the centre-of-mass energy.
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# BOSS-side electron PID recipe cannot be expressed through the tag signal_side declaration
# (the 'ep' key uses the fixed SimplePIDSvc thresholds).
alg.note(:pid_correction_method,
  "signal-side positron identified by E/p > 0.8 together with PID likelihoods " \
  "L'_e > 0.001 and L'_e/(L'_e + L'_pi + L'_K) > 0.8; the tag-layer 'ep' key " \
  "instead uses the fixed SimplePIDSvc thresholds, so this dedicated electron " \
  "selection must be applied in the generated BOSS code")

alg.apply
alg.execute_on(data_points + incMC_points + exMCs_signal)