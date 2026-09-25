### Dataset description ###
# Real data + inclusive MC at seven c.m. energies (4599.53 - 4698.82 MeV)
data_4600  = DatasetManager.real_data.find("703_4600")   # 4599.53 MeV (BOSS 703)
data_4610  = DatasetManager.real_data.find("706_4610")   # 4611.86 MeV (BOSS 706)
data_4620  = DatasetManager.real_data.find("706_4620")   # 4628.00 MeV
data_4640  = DatasetManager.real_data.find("706_4640")   # 4640.91 MeV
data_4660  = DatasetManager.real_data.find("706_4660")   # 4661.24 MeV
data_4680  = DatasetManager.real_data.find("706_4680")   # 4681.92 MeV
data_4700  = DatasetManager.real_data.find("706_4700")   # 4698.82 MeV

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")

data_points  = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680, data_4700]
incMC_points = [incMC_4600, incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

# Decay card for the double-tag signal: psi(4260) -> Lambda_c+ Lambda_c-
# (signal Lambda_c+ -> n K_S0 pi+ pi0, tag Lambda_c- -> pbar K_S0 pi+ pi0 as described for the MC),
# with K_S0 -> pi+ pi- and pi0 -> gamma gamma.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 n0 K_S0 pi+ pi0 PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K_S0 pi+ pi0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 1M exclusive MC events per energy point (same signal MC over the 7 scan points)
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_Lambdac_doubletag"   # auto-suffixed per energy point
  config.events        = 1_000_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Tag-based event selection (BOSS) ###
alg_name = "LambdacTagNKsPiPi0"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.640]})   # nominal c.m. energy of the scan
   .with_decay_card(decay_card_signal)

# Tag side: Lambda_c- reconstructed in eight single-tag modes.
# (Lambda_c double-tagging is not natively supported, so only a single tag side is declared.)
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP,            # p K_S0
          :LambdacPtoKPiP,           # p K+ pi-
          :LambdacPtoLambdaPi,       # Lambda pi-
          :LambdacPtoLambdaPiPi0,    # Lambda pi- pi0
          :LambdacPtoLambdaPiPiPi,   # Lambda pi- pi+ pi-
          :LambdacPtoSigmaPi0,       # Sigma- pi0
          :LambdacPtoSigma0Pi,       # Sigma0 pi-
          :LambdacPtoSigmaPiPi       # Sigma- pi+ pi-
  t.charm(-1)                        # pin the tagged side to Lambda_c-
end

# Signal side: Lambda_c+ -> n K_S0 pi+ pi0 -> (2 gamma)(missing n) [+ K_S0, pi+ handled outside DSL]
alg.signal_side do |s|
  s.photons 2                       # two photons from pi0
  s.min_photon_energy 0.025         # photon energy > 25 MeV
  s.min_photon_angle 10.0           # photon opening angle > 10 degrees
  s.missing :n, mass: 0.939565      # one missing neutron (massive form)
end

# Kinematic fit: four-momentum constraint against the measured CMS, chi2 < 200
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# BOSS-side procedures that cannot be expressed in the DSL
alg.note(:signal_ks_reconstruction,
    "signal-side K_S0 -> pi+ pi- is reconstructed with a secondary-vertex fit requiring " \
    "L/sigma_L > 2 together with a pi+pi- invariant-mass window; the tag-based signal side " \
    "cannot build a secondary-vertex K_S0, so this step is not encoded in BOSS")
   .note(:signal_pi0_reconstruction,
    "signal-side pi0 is selected with a gamma-gamma invariant-mass window 0.115-0.150 GeV/c^2 " \
    "and a 1C (gamma-gamma mass-constrained) Kalman fit with chi2 < 200; only the two-photon " \
    "requirement is expressible in the DSL signal side")
   .note(:direct_pion_pid,
    "the direct pi+ from Lambda_c+ is required to satisfy the kaon-veto PID L(pi) > L(K); this " \
    "per-track likelihood comparison is not expressible in the tag-based signal side")
   .note(:signal_vertex_cuts,
    "separate vertex cuts are applied to the K_S0 daughters (|Vz| < 20 cm) and to the direct pi+ " \
    "(|Vz| < 10 cm); not expressible in the tag-based signal side")
   .note(:background_veto,
    "post-fit peaking-background vetoes in the M(n pi+/-) and M(n pi0) distributions are applied " \
    "after the kinematic fit; handled at the ROOT-analysis level")
   .note(:tag_mode_limitation,
    "eight of the eleven paper tag modes have available DTagAlg channel symbols (pK_S0, p K+ pi-, " \
    "Lambda pi-, Lambda pi- pi0, Lambda pi- pi+ pi-, Sigma- pi0, Sigma0 pi-, Sigma- pi+ pi-); the " \
    "remaining three tag modes have no BOSS channel symbol and Lambda_c double-tagging is not " \
    "natively supported, so only the eight single-tag modes are covered")

alg.apply                                        # TagAnalysis: takes no Selection argument
alg.execute_on(data_points + incMC_points + exMCs_signal)